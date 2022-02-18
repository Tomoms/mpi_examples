#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include "mpi.h"
#include "kd_tree.h"

#define MMPI_TAG_METADATA	0
#define MMPI_TAG_DATA		1
#define MMPI_TAG_NODE		2

void perror_exit(char *msg);

void setup_node_type(MPI_Datatype *node_type)
{
	int lengths[4] = {1, 2, 1, 1};
	struct kdnode dummy;
	MPI_Aint displacements[4];
	MPI_Aint base_address;
	MPI_Get_address(&dummy, &base_address);
	MPI_Get_address(&dummy.axis, &displacements[0]);
	MPI_Get_address(&dummy.split, &displacements[1]);
	MPI_Get_address(&dummy.left, &displacements[2]);
	MPI_Get_address(&dummy.right, &displacements[3]);
	displacements[0] = MPI_Aint_diff(displacements[0], base_address);
	displacements[1] = MPI_Aint_diff(displacements[1], base_address);
	displacements[2] = MPI_Aint_diff(displacements[2], base_address);
	displacements[3] = MPI_Aint_diff(displacements[3], base_address);
	MPI_Datatype types[4] = {MPI_CHAR, MPI_DOUBLE, MPI_AINT, MPI_AINT};
	MPI_Type_create_struct(4, lengths, displacements, types, node_type);
	MPI_Type_commit(node_type);
}

void send_data(double *left, int left_len, int left_rank, double *right, int right_len, int right_rank)
{
	MPI_Request request_left, request_left_len, request_right, request_right_len;
	MPI_Status status_left, status_left_len, status_right, status_right_len;
	MPI_Isend(&left_len, 1, MPI_INT, left_rank, MMPI_TAG_METADATA, MPI_COMM_WORLD, &request_left_len);
	MPI_Wait(&request_left_len, &status_left_len);
	MPI_Isend(left, 2 * left_len, MPI_DOUBLE, left_rank, MMPI_TAG_DATA, MPI_COMM_WORLD, &request_left);
	MPI_Isend(&right_len, 1, MPI_INT, right_rank, MMPI_TAG_METADATA, MPI_COMM_WORLD, &request_right_len);
	MPI_Wait(&request_right_len, &status_right_len);
	MPI_Isend(right, 2 * right_len, MPI_DOUBLE, right_rank, MMPI_TAG_DATA, MPI_COMM_WORLD, &request_right);
	MPI_Wait(&request_left, &status_left);
	MPI_Wait(&request_right, &status_right);
}

void send_exit_signal(int left_rank, int right_rank)
{
	int dummy = -1;
	MPI_Request request_left_len, request_right_len;
	MPI_Isend(&dummy, 1, MPI_INT, left_rank, MMPI_TAG_METADATA, MPI_COMM_WORLD, &request_left_len);
	MPI_Isend(&dummy, 1, MPI_INT, right_rank, MMPI_TAG_METADATA, MPI_COMM_WORLD, &request_right_len);
}

int main(int argc, char *argv[])
{
	int nproc, rank, ndead;
	double *my_data;
	int my_data_len;
	struct kdnode *tree;
	struct kdnode my_node;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	MPI_Datatype node_type;
	setup_node_type(&node_type);

	int my_left_child = 2 * rank + 1;
	int my_right_child = my_left_child + 1;
	int my_parent = (int) ((rank - 1) / 2);

	if (!rank) {
		if (argc != 2) {
			printf("Exactly one command line argument is expected. Quitting.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
		my_data = load_points(argv[1]);
		my_data_len = DATASET_SIZE;
		tree = malloc(sizeof(struct kdnode) * nproc);
		if (!tree)
			perror_exit("malloc()");
	} else {
		MPI_Recv(&my_data_len, 1, MPI_INT, my_parent, MMPI_TAG_METADATA, MPI_COMM_WORLD, &status);
		if (my_data_len == -1) {
			my_node.axis = -2;
			MPI_Send(&my_node, 1, node_type, 0, MMPI_TAG_NODE, MPI_COMM_WORLD);
			MPI_Finalize();
			return 0;
		}
		my_data = malloc(sizeof(double) * my_data_len * 2);
		if (!my_data)
			perror_exit("malloc()");
		MPI_Recv(my_data, 2 * my_data_len, MPI_DOUBLE, my_parent, MMPI_TAG_DATA, MPI_COMM_WORLD, &status);
	}

	MPI_Bcast(&tree, 1, MPI_AINT, 0, MPI_COMM_WORLD);

	if (my_data_len > 1) {
		my_node = build_node(my_data, my_data_len);
		my_node.left = tree + my_left_child;
		my_node.right = tree + my_right_child;
		double *data_left, *data_right;
		int data_left_len, data_right_len;
		split_data(my_data, my_data_len, &data_left, &data_left_len, &data_right, &data_right_len, &my_node);
#ifdef VERBOSE
		for (int i = 0; i < data_left_len * 2; i++)
			printf("data_left[%d] = %lf\n", i, data_left[i]);
		for (int i = 0; i < data_right_len * 2; i++)
			printf("data_right[%d] = %lf\n", i, data_right[i]);
#endif
		send_data(data_left, data_left_len, my_left_child, data_right, data_right_len, my_right_child);
		free(data_left);
		free(data_right);
	} else {
		my_node.axis = -1;
		if (my_data_len == 1) {
			my_node.axis = 0;
			my_node.split[0] = my_data[0];
			my_node.split[1] = my_data[1];
#ifdef VERBOSE
			printf("process %d, terminal node: x = %lf, y = %lf\n", rank, my_data[0], my_data[1]);
#endif
			if (my_left_child < nproc && my_right_child < nproc)
				send_exit_signal(my_left_child, my_right_child);
		} else
#ifdef VERBOSE
			printf("process %d, terminal empty node\n", rank);
#endif
		my_node.left = NULL;
		my_node.right = NULL;
	}

	if (!rank) {
		*tree = my_node;
		for (int i = 1; i < nproc; i++)
			MPI_Recv(&tree[i], 1, node_type, i, MMPI_TAG_NODE, MPI_COMM_WORLD, &status);
		struct kdnode *n;
		for (int i = 0; i < nproc; i++) {
			n = &tree[i];
			if (n->axis != -2)
				printf("Node %d: address = %p, (%.1lf,%.1lf), left = %p, right = %p\n", i, n, n->split[0], n->split[1], n->left, n->right);
		}
	} else {
		MPI_Send(&my_node, 1, node_type, 0, MMPI_TAG_NODE, MPI_COMM_WORLD);
	}

	if (!rank)
		free(tree);
	free(my_data);
	MPI_Finalize();
	return 0;
}
