#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "kd_tree.h"

#define VERBOSE
#define NUM_THREADS	15

void perror_exit(char *msg);

struct kdnode *tree;

void build_node_and_create_tasks(double *data, int len, int node_ordinal)
{
	int id = omp_get_thread_num();
	struct kdnode my_node;
	if (len > 1) {
		my_node = build_node(data, len);
		double *data_left, *data_right;
		int data_left_len, data_right_len;
		int left_child_ordinal = 2 * node_ordinal + 1;
		int right_child_ordinal = left_child_ordinal + 1;
		my_node.left = tree + left_child_ordinal;
		my_node.right = tree + right_child_ordinal;
#ifdef VERBOSE
		printf("thread %d built node: x = %.2lf, y = %.2lf, left = %p, right = %p, and is about to write it at position %d of tree\n",
			   id, my_node.split[0], my_node.split[1], my_node.left, my_node.right, node_ordinal);
#endif
		split_data(data, len, &data_left, &data_left_len, &data_right, &data_right_len, &my_node);
		if (data_left_len >= 0) {
			#pragma omp task
			{
				build_node_and_create_tasks(data_left, data_left_len, left_child_ordinal);
			}
		}
		if (data_right_len >= 0) {
			#pragma omp task
			{
				build_node_and_create_tasks(data_right, data_right_len, right_child_ordinal);
			}
		}
	} else {
		my_node.axis = -1;
		if (len == 1) {
			my_node.axis = 0;
			my_node.split[0] = data[0];
			my_node.split[1] = data[1];
#ifdef VERBOSE
			printf("thread %d, terminal node: x = %lf, y = %lf, writing it ad %d\n", id, data[0], data[1], node_ordinal);
#endif
		} else {
#ifdef VERBOSE
			printf("thread %d, terminal empty node, writing it ad %d\n", id, node_ordinal);
#endif
		}
		my_node.left = NULL;
		my_node.right = NULL;
	}
	tree[node_ordinal] = my_node;
	free(data);
}

int main(int argc, char *argv[])
{
	if (argc != 2) {
		printf("Exactly one command line argument is expected. Quitting.\n");
		return 1;
	}

	double *data = load_points(argv[1]);
	int len = DATASET_SIZE;

	tree = malloc(sizeof(struct kdnode) * NUM_THREADS);
	if (!tree) {
		perror_exit("malloc()");
	}

	#pragma omp parallel
	{
		#pragma omp single
			build_node_and_create_tasks(data, len, 0);
	}
	free(tree);
	return 0;
}
