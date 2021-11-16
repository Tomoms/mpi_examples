#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

static inline int finalize_exit(int code)
{
    MPI_Finalize();
    exit(code);
}

int main(int argc, char *argv[])
{
    int rank, nproc;
    int n;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!rank) {
        n = atoi(argv[1]);
        /*
         * In this program we assume n > 0, thus if
         * atoi() returns 0 it is surely an error.
         */
        if (n <= 0)
            finalize_exit(1);
        /* Un-buffer stdout
        int ret = setvbuf(stdout, 0, _IONBF, BUFSIZ);
        if (ret) {
            finalize_exit(ret);
        }
        */
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double base = 1.0 / (double) n;
    double x, y = 0;
    for (int i = rank; i < n; i += nproc) {
        x = base * i + base / 2;
        y += 4 / (1 + x * x);
    }
    double sum = y * base;
    double pi;
    
    MPI_Reduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (!rank)
        printf("The result for n = %d is pi = %.10f\n", n, pi);

    finalize_exit(0);
}
