#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
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

double *load_points(char *filename)
{
	FILE *data_file = fopen(filename, "r");
	if (data_file == NULL) {
		perror("fopen()");
		MPI_Abort(MPI_COMM_WORLD, 1);
		exit(EXIT_FAILURE);
	}
	double *data = calloc(2 * DATASET_SIZE, sizeof(*data));
	double *x = data, *y = data + DATASET_SIZE;
	for (; x < data + DATASET_SIZE; x++, y++) {
		fscanf(data_file, "%lf", x);
		fscanf(data_file, "%lf", y);
	}
	if (fclose(data_file)) {
		perror("fclose()");
		MPI_Abort(MPI_COMM_WORLD, 1);
		exit(EXIT_FAILURE);
	}
#if VERBOSE
	for (double *datar = data; datar < data + DATASET_SIZE; datar++)
		printf("point x = %lf; y = %lf\n", *datar, *(datar + DATASET_SIZE));
#endif
	return data;
}

double get_mean(double *data, char axis)
{
	double *startpos = data + axis * DATASET_SIZE;
	double *endpos = data + (axis + 1) * DATASET_SIZE;
	double mean = 0;
	for (double *scanner = startpos; scanner < endpos; scanner++)
		mean += *scanner;
	mean /= DATASET_SIZE;
	return mean;
}

// TODO remove useless division and squaring
double get_variance(double *data, char axis, double mean)
{
	double *startpos = data + axis * DATASET_SIZE;
	double *endpos = data + (axis + 1) * DATASET_SIZE;
	double variance = 0;
	for (double *scanner = startpos; scanner < endpos; scanner++)
		variance += (*scanner - mean) * (*scanner - mean);
	variance /= DATASET_SIZE;
	return variance;
}

int get_median_index(double *data, char axis, double mean)
{
	int offset = axis * DATASET_SIZE;
	double delta = DBL_MAX, candidate_delta;
	int ret;
	for (int i = 0; i < DATASET_SIZE; i++) {
		if (candidate_delta = fabs(data[i + offset] - mean) < delta) {
			delta = candidate_delta;
			ret = i;
		}
	}
	return ret;
}

int main(int argc, char *argv[])
{
	int nproc, rank;
	int my_first_child = 2 * rank + 1;
	int my_second_child = my_first_child + 1;
	int my_parent = (int) (rank / 2);
	double *my_data;
	int my_data_len;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (!rank) {
		if (argc != 2) {
			printf("Exactly one command line argument is expected. Quitting.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
		my_data = load_points(argv[1]);
		my_data_len = DATASET_SIZE;
	}

	if (!rank) {
		double mean[2];
		double variance[2];
		mean[0] = get_mean(my_data, 0);
		variance[0] = get_variance(my_data, 0, mean[0]);
		mean[1] = get_mean(my_data, 1);
		variance[1] = get_variance(my_data, 1, mean[1]);
#if VERBOSE
		printf("mean: x = %lf; y = %lf\n", mean[0], mean[1]);
		printf("variance: x = %lf; y = %lf\n", variance[0], variance[1]);
#endif
		char max_variance_axis = variance[0] < variance[1];
		int median_idx = get_median_index(my_data, max_variance_axis, mean[max_variance_axis]);
#if VERBOSE
		printf("median: %lf\n", my_data[median_idx + max_variance_axis * DATASET_SIZE]);
#endif
		struct kdnode my_node;
		my_node.axis = max_variance_axis;
		my_node.split[0] = my_data[median_idx];
		my_node.split[1] = my_data[median_idx + DATASET_SIZE];
	}

	MPI_Finalize();
	return 0;
}
