#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

struct kdnode {
	char axis;
	double split[2];
	int ordinal;
	struct kdnode *left, *right;
};

static inline void perror_exit(char *msg)
{
	perror(msg);
#ifdef MPI
	MPI_Abort(MPI_COMM_WORLD, 1);
#endif
	exit(EXIT_FAILURE);
}

int count_points(char *filename)
{
	FILE *file = fopen(filename, "r");
	if (!file)
		perror_exit("fopen()");

	int lines = 0;
	while (EOF != (fscanf(file, "%*[^\n]"), fscanf(file, "%*c")))
		lines++;

	if (fclose(file))
		perror_exit("fclose()");

	return lines;
}

double *load_points(char *filename, int npoints)
{
	FILE *data_file = fopen(filename, "r");
	if (data_file == NULL)
		perror_exit("fopen()");
	double *data = calloc(2 * npoints, sizeof(*data));
	if (!data)
		perror_exit("calloc()");
	double *x = data, *y = data + npoints;
	for (; x < data + npoints; x++, y++) {
		fscanf(data_file, "%lf", x);
		fscanf(data_file, "%lf", y);
	}
	if (fclose(data_file))
		perror_exit("fclose()");
#ifdef VERBOSE
	for (double *datar = data; datar < data + npoints; datar++)
		printf("point x = %lf; y = %lf\n", *datar, *(datar + npoints));
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
	double contribution;
	for (double *scanner = startpos; scanner < endpos; scanner++) {
		contribution = fabs(*scanner - mean);
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

struct kdnode build_node(double *my_data, int my_data_len, int ordinal)
{
	struct kdnode my_node;

	if (!my_data_len) {
		my_node.axis = -1;
		my_node.split[0] = 0;
		my_node.split[1] = 0;
		my_node.ordinal = ordinal;
		my_node.left = NULL;
		my_node.right = NULL;
		goto print_exit;
	}

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
	my_node.axis = max_variance_axis;
	my_node.ordinal = ordinal;
	my_node.split[0] = my_data[median_idx];
	my_node.split[1] = my_data[median_idx + my_data_len];
print_exit:
#ifdef VERBOSE_RES
	printf("node: axis = %d; ordinal = %d; point: (%lf,%lf)\n", my_node.axis, my_node.ordinal, my_node.split[0], my_node.split[1]);
#endif
	return my_node;
}

void split_data(double *data, int len, double **left, int *left_len, double **right,
				int *right_len, struct kdnode *node)
{
	double *startpos = data + node->axis * len;
	double *endpos = data + (node->axis + 1) * len;
	int data_right_len = 0;
	for (double *scanner = startpos; scanner < endpos; scanner++) {
		if (*scanner > node->split[node->axis])
			data_right_len++;
	}
	*right_len = data_right_len;
	int data_left_len = len - data_right_len - 1;
	*left_len = data_left_len;
	*left = malloc(sizeof(double) * data_left_len * 2);
	*right = malloc(sizeof(double) * data_right_len * 2);
	if (!(*left) || !(*right))
		perror_exit("malloc()");
	int rindex = 0, lindex = 0;
	for (int i = 0; i < len; i++) {
		if (data[i + node->axis * len] > node->split[node->axis]) {
			(*right)[rindex] = data[i];
			(*right)[rindex + data_right_len] = data[i + len];
			rindex++;
		} else {
			if (data[i] == node->split[0] && data[i + len] == node->split[1])
				continue;
			(*left)[lindex] = data[i];
			(*left)[lindex + data_left_len] = data[i + len];
			lindex++;
		}
	}
}

int compute_total_nodes(int npoints)
{
	if (npoints <= 1)
		return 1;
	int left = (int) ((npoints - 1) / 2);
	int right = npoints - 1 - left;
	return 1 + compute_total_nodes(left) + compute_total_nodes(right);
}

static inline struct kdnode *extend_tree(struct kdnode *tree, size_t new_size)
{
	struct kdnode *new_tree = realloc(tree, new_size * sizeof(struct kdnode));
	if (!new_tree)
		perror_exit("realloc()");
	return new_tree;
}

static inline int get_left_child(int index)
{
	return 2 * index + 1;
}

static inline int get_right_child(int index)
{
	return (index + 1) * 2;
}

static inline int get_parent(int index)
{
	return (int) ((index - 1) / 2);
}

static inline void pointer_fixup(struct kdnode *tree, int node_index)
{
	int parent_index = get_parent(node_index);
	if (get_left_child(parent_index) == node_index)
		tree[parent_index].left = &tree[node_index];
	else
		tree[parent_index].right = &tree[node_index];
}
