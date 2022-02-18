#include <math.h>
#include <float.h>

#define DATASET_SIZE	10

struct kdnode {
	char axis;
	double split[2];
	struct kdnode *left, *right;
};

inline void perror_exit(char *msg)
{
	perror(msg);
#ifdef MPI
	MPI_Abort(MPI_COMM_WORLD, 1);
#endif
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
#ifdef VERBOSE
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
	for (double *scanner = startpos; scanner < endpos; scanner++) {
		double contribution = fabs(*scanner - mean);
		variance += contribution
#ifdef REAL_VARIANCE
		* contribution
#endif
		;
	}
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
#ifdef VERBOSE
	printf("mean: x = %lf; y = %lf\n", mean[0], mean[1]);
	printf("variance: x = %lf; y = %lf\n", variance[0], variance[1]);
	printf("median: %lf\n", my_data[median_idx + max_variance_axis * my_data_len]);
#endif
	struct kdnode my_node;
	my_node.axis = max_variance_axis;
	my_node.split[0] = my_data[median_idx];
	my_node.split[1] = my_data[median_idx + my_data_len];
#ifdef VERBOSE
	printf("node: axis = %d; point x = %lf, y = %lf\n", my_node.axis, my_node.split[0], my_node.split[1]);
#endif
	return my_node;
}
