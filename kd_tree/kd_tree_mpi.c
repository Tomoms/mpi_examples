#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include "mpi.h"

#define DATASET_SIZE	10
#define REAL_VARIANCE
#ifndef VERBOSE
#define VERBOSE			1
#endif

#define MMPI_TAG_METADATA	0
#define MMPI_TAG_DATA		1
#define MMPI_TAG_NODE		2

struct kdnode {
	char axis;
	double split[2];
	struct kdnode *left, *right;
};

void perror_exit(char *msg)
{
	perror(msg);
	MPI_Abort(MPI_COMM_WORLD, 1);
	exit(EXIT_FAILURE);
}

double *load_points(char *filename)
{
	FILE *data_file = fopen(filename, "r");
	if (data_file == NULL)
		perror_exit("fopen()");
	double *data = calloc(2 * DATASET_SIZE, sizeof(*data));
	if (!data)
		perror_exit("calloc()");
	double *x = data, *y = data + DATASET_SIZE;
	for (; x < data + DATASET_SIZE; x++, y++) {
		fscanf(data_file, "%lf", x);
		fscanf(data_file, "%lf", y);
	}
	if (fclose(data_file))
		perror_exit("fclose()");
#if VERBOSE
	for (double *datar = data; datar < data + DATASET_SIZE; datar++)
		printf("point x = %lf; y = %lf\n", *datar, *(datar + DATASET_SIZE));
#endif
	return data;
}

double get_mean(double *data, int len, char axis)
{
	double *startpos = data + axis * len;
	double *endpos = data + (axis + 1) * len;
	double mean = 0;
	for (double *scanner = startpos; scanner < endpos; scanner++)
		mean += *scanner;
	mean /= len;
	return mean;
}

double get_variance(double *data, int len, char axis, double mean)
{
	double *startpos = data + axis * len;
	double *endpos = data + (axis + 1) * len;
	double variance = 0;
	for (double *scanner = startpos; scanner < endpos; scanner++)
		variance += (*scanner - mean)
#ifdef REAL_VARIANCE
		* (*scanner - mean)
#endif
		;
#ifdef REAL_VARIANCE
	variance /= len;
#endif
	return variance;
}

int get_median_index(double *data, int len, char axis, double mean)
{
	int offset = axis * len;
	double delta = DBL_MAX, candidate_delta;
	int index;
	for (int i = 0; i < len; i++) {
		candidate_delta = fabs(data[i + offset] - mean);
		if (candidate_delta < delta) {
			delta = candidate_delta;
			index = i;
		}
	}
	return index;
}

struct kdnode build_node(double *my_data, int my_data_len)
{
	double mean[2];
	double variance[2];
	mean[0] = get_mean(my_data, my_data_len, 0);
	variance[0] = get_variance(my_data, my_data_len, 0, mean[0]);
	mean[1] = get_mean(my_data, my_data_len, 1);
	variance[1] = get_variance(my_data, my_data_len, 1, mean[1]);
	char max_variance_axis = variance[0] < variance[1];
	int median_idx = get_median_index(my_data, my_data_len, max_variance_axis, mean[max_variance_axis]);
#if VERBOSE
	printf("mean: x = %lf; y = %lf\n", mean[0], mean[1]);
	printf("variance: x = %lf; y = %lf\n", variance[0], variance[1]);
	printf("median: %lf\n", my_data[median_idx + max_variance_axis * my_data_len]);
#endif
	struct kdnode my_node;
	my_node.axis = max_variance_axis;
	my_node.split[0] = my_data[median_idx];
	my_node.split[1] = my_data[median_idx + my_data_len];
#if VERBOSE
	printf("node: axis = %d; point x = %lf, y = %lf\n", my_node.axis, my_node.split[0], my_node.split[1]);
#endif
	return my_node;
}

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

void split_data(double *data, int len, double *left, int left_len, double *right, int right_len)
{

}

int main(int argc, char *argv[])
{
	int nproc, rank;
	double *my_data;
	int my_data_len;
	struct kdnode *tree;

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
		/*
		unsigned int tree_levels = (unsigned int) ((log2(my_data_len) + 1e-10) + 1);
		unsigned long lower = 1 << (tree_levels - 1);
		unsigned long upper = 1 << tree_levels;
		tree_levels += lower > my_data_len;
		printf("lower = %ld; upper = %ld;\n", lower, upper);
		*/
		tree = malloc(sizeof(struct kdnode));
		if (!tree)
			perror_exit("malloc()");
	}

	if (!rank) {
		struct kdnode my_node = build_node(my_data, my_data_len);
		double *startpos = my_data + my_node.axis * my_data_len;
		double *endpos = my_data + (my_node.axis + 1) * my_data_len;
		int data_right_len = 0;
		for (double *scanner = startpos; scanner < endpos; scanner++) {
			if (*scanner > my_node.split[my_node.axis])
				data_right_len++;
		}
		int data_left_len = my_data_len - data_right_len - 1;
		double *data_left = malloc(sizeof(double) * data_left_len * 2);
		double *data_right = malloc(sizeof(double) * data_right_len * 2);
		if (!data_left || !data_right)
			perror_exit("malloc()");
		int rindex = 0, lindex = 0;
		for (int i = 0; i < my_data_len; i++) {
			if (my_data[i + my_node.axis * my_data_len] > my_node.split[my_node.axis]) {
				data_right[rindex] = my_data[i];
				data_right[rindex + data_right_len] = my_data[i + my_data_len];
				rindex++;
			} else {
				if (my_data[i] == my_node.split[0] && my_data[i + my_data_len] == my_node.split[1])
					continue;
				data_left[lindex] = my_data[i];
				data_left[lindex + data_left_len] = my_data[i + my_data_len];
				lindex++;
			}
		}
#if VERBOSE
		for (int i = 0; i < data_left_len * 2; i++)
			printf("data_left[%d] = %lf\n", i, data_left[i]);
		for (int i = 0; i < data_right_len * 2; i++)
			printf("data_right[%d] = %lf\n", i, data_right[i]);
#endif
		MPI_Request trash;
		MPI_Request req_left, req_right;
		MPI_Status status_left, status_right;
		MPI_Isend(&data_left_len, 1, MPI_INT, my_left_child, MMPI_TAG_METADATA, MPI_COMM_WORLD, &trash);
		MPI_Isend(data_left, 2 * data_left_len, MPI_DOUBLE, my_left_child, MMPI_TAG_DATA, MPI_COMM_WORLD, &req_left);
		MPI_Isend(&data_right_len, 1, MPI_INT, my_right_child, MMPI_TAG_METADATA, MPI_COMM_WORLD, &trash);
		MPI_Isend(data_right, 2 * data_right_len, MPI_DOUBLE, my_right_child, MMPI_TAG_DATA, MPI_COMM_WORLD, &req_right);
		MPI_Wait(&req_left, &status_left);
		MPI_Wait(&req_right, &status_right);
		free(data_left);
		free(data_right);
	} else {
		MPI_Recv(&my_data_len, 1, MPI_INT, my_parent, MMPI_TAG_METADATA, MPI_COMM_WORLD, &status);
		my_data = malloc(sizeof(double) * my_data_len * 2);
		if (!my_data)
			perror_exit("malloc()");
		MPI_Recv(my_data, 2 * my_data_len, MPI_DOUBLE, my_parent, MMPI_TAG_DATA, MPI_COMM_WORLD, &status);
		struct kdnode my_node = build_node(my_data, my_data_len);
	}

	free(my_data);
	MPI_Finalize();
	return 0;
}
