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
#ifndef VERBOSE
#define VERBOSE			1
#endif

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

// TODO remove useless division and squaring
double get_variance(double *data, int len, char axis, double mean)
{
	double *startpos = data + axis * len;
	double *endpos = data + (axis + 1) * len;
	double variance = 0;
	for (double *scanner = startpos; scanner < endpos; scanner++)
		variance += (*scanner - mean) * (*scanner - mean);
	variance /= len;
	return variance;
}

int get_median_index(double *data, int len, char axis, double mean)
{
	int offset = axis * len;
	double delta = DBL_MAX, candidate_delta;
	int ret;
	for (int i = 0; i < len; i++) {
		candidate_delta = fabs(data[i + offset] - mean);
		if (candidate_delta < delta) {
			delta = candidate_delta;
			ret = i;
		}
	}
	return ret;
}

struct kdnode *tree_shm_mmap(int rank)
{
	int oflags = rank ? O_RDWR : O_CREAT | O_RDWR | O_TRUNC;
	int fd = shm_open("/tree", oflags, 00666);
	if (fd == -1)
		perror_exit("shm_open()");
	struct kdnode *tree = (struct kdnode *) mmap(NULL, sizeof(struct kdnode) * nproc, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	close(fd);
	if (tree == MAP_FAILED)
		perror_exit("mmap()");
	return tree;
}

struct kdnode *build_node(double *my_data, int my_data_len)
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
	struct kdnode *my_node = (struct kdnode *) mmap(NULL, sizeof(struct kdnode), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
	my_node->axis = max_variance_axis;
	my_node->split[0] = my_data[median_idx];
	my_node->split[1] = my_data[median_idx + my_data_len];
#if VERBOSE
	printf("node: axis = %d; point x = %lf, y = %lf\n", my_node->axis, my_node->split[0], my_node->split[1]);
#endif
	//MPI_Send(&my_node, 1, MPI_LONG, 1, 0, MPI_COMM_WORLD);
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

	int my_first_child = 2 * rank + 1;
	int my_second_child = my_first_child + 1;
	int my_parent = (int) (rank / 2);

	if (!rank) {
		if (argc != 2) {
			printf("Exactly one command line argument is expected. Quitting.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
		my_data = load_points(argv[1]);
		my_data_len = DATASET_SIZE;

		tree = tree_shm_mmap(rank);
		MPI_Barrier();
	} else {
		MPI_Barrier();
		tree = tree_shm_mmap(rank);
	}

	if (!rank) {
		struct kdnode *my_node = build_node(my_data, my_data_len);
	} else {
		/*
		uintptr_t pointer;
		MPI_Recv(&pointer, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
		struct kdnode *father_node = (struct kdnode *) pointer;
		printf("node: axis = %d; point x = %lf, y = %lf\n", father_node->axis, father_node->split[0], father_node->split[1]);
		*/
	}

	MPI_Finalize();
	return 0;
}
