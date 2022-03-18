#include <stdbool.h>
#include <string.h>
#include <alloca.h>
#include <sys/param.h>
#include "kd_tree.h"

#define MMPI_TAG_METADATA	0
#define MMPI_TAG_DATA		1
#define MMPI_TAG_NODE		2
#define MMPI_TAG_EXIT		3

MPI_Datatype node_type, metadata_type;
MPI_Comm ring_comm;

struct task_metadata {
	int node_index;
	int data_size;
};

void setup_node_type(MPI_Datatype *node_type)
{
	int lengths[5] = {1, 2, 1, 1, 1};
	struct kdnode dummy;
	MPI_Aint displacements[5];
	MPI_Aint base_address;
	MPI_Get_address(&dummy, &base_address);
	MPI_Get_address(&dummy.axis, &displacements[0]);
	MPI_Get_address(&dummy.split, &displacements[1]);
	MPI_Get_address(&dummy.ordinal, &displacements[2]);
	MPI_Get_address(&dummy.left, &displacements[3]);
	MPI_Get_address(&dummy.right, &displacements[4]);
	displacements[0] = MPI_Aint_diff(displacements[0], base_address);
	displacements[1] = MPI_Aint_diff(displacements[1], base_address);
	displacements[2] = MPI_Aint_diff(displacements[2], base_address);
	displacements[3] = MPI_Aint_diff(displacements[3], base_address);
	displacements[4] = MPI_Aint_diff(displacements[4], base_address);
	MPI_Datatype types[5] = {MPI_CHAR, MPI_DOUBLE, MPI_INT, MPI_AINT, MPI_AINT};
	MPI_Type_create_struct(5, lengths, displacements, types, node_type);
	MPI_Type_commit(node_type);
}

void setup_metadata_type(MPI_Datatype *metadata_type)
{
	int lengths[2] = {1, 1};
	struct task_metadata dummy;
	MPI_Aint displacements[2];
	MPI_Aint base_address;
	MPI_Get_address(&dummy, &base_address);
	MPI_Get_address(&dummy.node_index, &displacements[0]);
	MPI_Get_address(&dummy.data_size, &displacements[1]);
	displacements[0] = MPI_Aint_diff(displacements[0], base_address);
	displacements[1] = MPI_Aint_diff(displacements[1], base_address);
	MPI_Datatype types[2] = {MPI_INT, MPI_INT};
	MPI_Type_create_struct(2, lengths, displacements, types, metadata_type);
	MPI_Type_commit(metadata_type);
}

void send_data(double *left, int left_len, int left_rank, double *right, int right_len, int right_rank)
{
	MPI_Request request_left, request_right;
	MPI_Status status_left, status_right;
	MPI_Isend(left, 2 * left_len, MPI_DOUBLE, left_rank, MMPI_TAG_DATA, ring_comm, &request_left);
	MPI_Isend(right, 2 * right_len, MPI_DOUBLE, right_rank, MMPI_TAG_DATA, ring_comm, &request_right);
	MPI_Wait(&request_left, &status_left);
	MPI_Wait(&request_right, &status_right);
}

void build_metadata(int my_node_index, int left_len, int right_len, struct task_metadata *meta_left, struct task_metadata *meta_right)
{
	meta_left->node_index = get_left_child(my_node_index);
	meta_left->data_size = left_len;
	meta_right->node_index = get_right_child(my_node_index);
	meta_right->data_size = right_len;
}

void send_metadata(struct task_metadata *meta_left, int left_rank, struct task_metadata *meta_right, int right_rank)
{
	MPI_Request request_left, request_right;
	MPI_Status status_left, status_right;
	MPI_Isend(meta_left, 1, metadata_type, left_rank, MMPI_TAG_METADATA, ring_comm, &request_left);
	MPI_Isend(meta_right, 1, metadata_type, right_rank, MMPI_TAG_METADATA, ring_comm, &request_right);
	MPI_Wait(&request_left, &status_left);
	MPI_Wait(&request_right, &status_right);
}

int main(int argc, char *argv[])
{
	int rank, nproc;
	const int one[1] = {1};
	int neighbors[2];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Status status;

	int ret = MPI_Cart_create(MPI_COMM_WORLD, 1, &nproc, one, 1, &ring_comm);
	if (ret) {
		perror("MPI_Cart_create()");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		return EXIT_FAILURE;
	}

	MPI_Cart_shift(ring_comm, 0, 1, &neighbors[0], &neighbors[1]);
	MPI_Comm_rank(ring_comm, &rank);
	// rank 0 is only responsible of collecting nodes
	(!neighbors[0] && (neighbors[0] = nproc - 1));
	(!neighbors[1] && (neighbors[1] = 1));

	setup_node_type(&node_type);
	setup_metadata_type(&metadata_type);

	double *my_data = NULL;
	double *data_left = NULL, *data_right = NULL;
	int my_data_len, data_left_len, data_right_len, total_nodes;
	size_t tree_size;
	struct kdnode my_node, *tree;
	int node_index = 0;
	if (!rank) {
		if (argc != 2) {
			printf("Exactly one command line argument is expected. Quitting.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
		my_data_len = count_points(argv[1]);
		my_data = load_points(argv[1], my_data_len);
		total_nodes = compute_total_nodes(my_data_len);
		tree_size = compute_tree_size(total_nodes);
		my_node = build_node(my_data, my_data_len, node_index);
		tree = malloc(tree_size * sizeof(struct kdnode));
		if (!tree)
			perror_exit("malloc()");
		bool *valid_indexes = calloc(tree_size, sizeof(bool));
		if (!valid_indexes)
			perror_exit("calloc()");
		tree[node_index] = my_node;
		valid_indexes[node_index] = 1;
		if (my_data_len > 1) {
			split_data(my_data, my_data_len, &data_left, &data_left_len, &data_right, &data_right_len, &my_node);
			struct task_metadata meta_left, meta_right;
			build_metadata(node_index, data_left_len, data_right_len, &meta_left, &meta_right);
			send_metadata(&meta_left, neighbors[0], &meta_right, neighbors[1]);
			send_data(data_left, data_left_len, neighbors[0], data_right, data_right_len, neighbors[1]);
			free(data_left);
			free(data_right);
			free(my_data);
		}
		int recv_nodes = 1; //  my own
		// suboptimal - this leaves holes where non-existing nodes would be, except for final non-existing nodes
		while (recv_nodes < total_nodes) {
			struct kdnode buf;
			MPI_Recv(&buf, 1, node_type, MPI_ANY_SOURCE, MMPI_TAG_NODE, ring_comm, &status);
			node_index = buf.ordinal;
			/*
			if (node_index + 1 > current_tree_size) {
				tree = extend_tree(tree, node_index + 1);
				valid_indexes = extend_indexes(valid_indexes, node_index + 1);
				current_tree_size = node_index + 1;
			}
			*/
			tree[node_index] = buf;
			valid_indexes[node_index] = 1;
			recv_nodes++;
		}
		for (int i = 1; i < nproc; i++)
			MPI_Send(NULL, 0, MPI_CHAR, i, MMPI_TAG_EXIT, ring_comm);

		fixup_pointers(tree, tree_size, valid_indexes);
		print_tree(tree, tree_size, valid_indexes);
		free(tree);
		free(valid_indexes);
	} else {
		bool should_exit = 0;
		while (!should_exit) {
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, ring_comm, &status);
			//printf("rank %d probed something: status.MPI_SOURCE = %d, Tag= %d\n", rank, status.MPI_SOURCE, status.MPI_TAG);
			if (status.MPI_SOURCE == 0 && status.MPI_TAG == MMPI_TAG_EXIT) {
				MPI_Recv(NULL, 0, MPI_CHAR, 0, MMPI_TAG_EXIT, ring_comm, &status);
				should_exit = 1;
				continue; // re-eval condition which should lead to breaking out of the loop
			} else {
				struct task_metadata my_metadata;
				MPI_Recv(&my_metadata, 1, metadata_type, MPI_ANY_SOURCE, MMPI_TAG_METADATA, ring_comm, &status);
				my_data_len = my_metadata.data_size;
				node_index = my_metadata.node_index;
				if (my_data_len) {
					my_data = malloc(my_data_len * sizeof(double) * 2);
					if (!my_data)
						perror_exit("malloc()");
					MPI_Recv(my_data, my_data_len * 2, MPI_DOUBLE, status.MPI_SOURCE, MMPI_TAG_DATA, ring_comm, &status);
					/*printf("rank %d has data for ordinal %d: ", rank, node_index);
					for(int i = 0;  i < my_data_len; i++) {
						printf("(%.1lf, %.1lf)\t", my_data[i], my_data[i + my_data_len]);
					}
					printf("\n");
					*/
					my_node = build_node(my_data, my_data_len, node_index);
					if (my_data_len > 1) {
						split_data(my_data, my_data_len, &data_left, &data_left_len, &data_right, &data_right_len, &my_node);
						/*for (int i = 0; i < data_left_len; i++)
							printf("left[%d]: (%.1lf, %.1lf)\t", i, data_left[i], data_left[i + data_left_len]);
						printf("\n");
						for (int i = 0; i < data_right_len; i++)
							printf("right[%d]: (%.1lf, %.1lf)\t", i, data_right[i], data_right[i + data_right_len]);
						printf("\n");
						*/
						struct task_metadata meta_left, meta_right;
						build_metadata(node_index, data_left_len, data_right_len, &meta_left, &meta_right);
						send_metadata(&meta_left, neighbors[0], &meta_right, neighbors[1]);
						//printf("rank %d build node (%.1lf, %.1lf) and is now sending data\n", rank, my_node.split[0], my_node.split[1]);
						send_data(data_left, data_left_len, neighbors[0], data_right, data_right_len, neighbors[1]);
						free(data_left);
						free(data_right);
						free(my_data);
					}
				} else { // empty terminal node
					//printf("rank %d has no data for ordinal %d\n", rank, node_index);
					MPI_Recv(NULL, 0, MPI_DOUBLE, status.MPI_SOURCE, MMPI_TAG_DATA, ring_comm, &status);
					my_node = build_node(NULL, my_data_len, node_index);
				}
				//printf("rank %d sent node\n", rank);
				MPI_Send(&my_node, 1, node_type, 0, MMPI_TAG_NODE, ring_comm);
			}
		}
	}

	//printf("rank %d exiting\n", rank);
	MPI_Finalize();
	return 0;
}
