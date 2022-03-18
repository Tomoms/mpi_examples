#include <omp.h>
#include "kd_tree.h"

struct kdnode *tree = NULL;
bool *valid_indexes = NULL;

void build_node_and_create_tasks(double *data, int len, int ordinal)
{
#ifdef VERBOSE
	int id = omp_get_thread_num();
#endif
	struct kdnode my_node;
	my_node = build_node(data, len, ordinal);
	if (len > 1) {
		double *data_left, *data_right;
		int data_left_len, data_right_len;
#ifdef VERBOSE
		printf("thread %d built node: x = %.2lf, y = %.2lf, left = %p, right = %p, and is about to write it at position %d of tree\n",
			   id, my_node.split[0], my_node.split[1], my_node.left, my_node.right, ordinal);
#endif
		split_data(data, len, &data_left, &data_left_len, &data_right, &data_right_len, &my_node);
		if (data_left_len >= 0) {
			#pragma omp task
			{
				build_node_and_create_tasks(data_left, data_left_len, get_left_child(ordinal));
			}
		}
		if (data_right_len >= 0) {
			#pragma omp task
			{
				build_node_and_create_tasks(data_right, data_right_len, get_right_child(ordinal));
			}
		}
	}
	tree[ordinal] = my_node;
	valid_indexes[ordinal] = 1;
	free(data);
}

int main(int argc, char *argv[])
{
	if (argc != 2) {
		printf("Exactly one command line argument is expected. Quitting.\n");
		return 1;
	}

	int npoints = count_points(argv[1]);
	double *data = load_points(argv[1], npoints);
	int total_nodes = compute_total_nodes(npoints);
	size_t tree_size = compute_tree_size(total_nodes);

	tree = malloc(tree_size * sizeof(struct kdnode));
	if (!tree)
		perror_exit("malloc()");
	valid_indexes = calloc(tree_size, sizeof(bool));
	if (!valid_indexes)
		perror_exit("calloc()");

	#pragma omp parallel
	{
		#pragma omp single
			build_node_and_create_tasks(data, npoints, 0);
	}

	fixup_pointers(tree, tree_size, valid_indexes);
	print_tree(tree, tree_size, valid_indexes);
	free(tree);
	free(valid_indexes);
	return 0;
}
