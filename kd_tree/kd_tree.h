#define DATASET_SIZE	10
#ifndef VERBOSE
#define VERBOSE
#endif

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
