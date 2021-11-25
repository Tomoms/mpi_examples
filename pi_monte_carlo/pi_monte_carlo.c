#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "mpi.h"

#define TAG_WORKER_READY    0
#define TAG_PAYLOAD         1
#define TAG_EXIT            2

static inline int finalize_exit(int code)
{
    MPI_Finalize();
    exit(code);
}

int main(int argc, char *argv[])
{
    int rank, nproc;
    int raw_coords[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    MPI_Comm workers;
    MPI_Status status;

    // process 0 is the server (rng)
    if (!rank) {
        int n = atoi(argv[1]);
        if (n <= 0)
            finalize_exit(1);

        /*
         * Dummy call to move the server into his own
         * communicator, misleadingly named 'workers'.
         * Need to check if this call is compulsory.
         */
        MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &workers);
        srand(2021);
        while (n) {
            MPI_Recv(MPI_BOTTOM, 0, MPI_INT, MPI_ANY_SOURCE, TAG_WORKER_READY, MPI_COMM_WORLD, &status);
            raw_coords[0] = rand();
            raw_coords[1] = rand();
            MPI_Send(&raw_coords[0], 2, MPI_INT, status.MPI_SOURCE, TAG_PAYLOAD, MPI_COMM_WORLD);
            n--;
        }

        for (int i = 1; i < nproc; i++) {
            MPI_Recv(MPI_BOTTOM, 0, MPI_INT, MPI_ANY_SOURCE, TAG_WORKER_READY, MPI_COMM_WORLD, &status);
            MPI_Send(MPI_BOTTOM, 0, MPI_INT, status.MPI_SOURCE, TAG_EXIT, MPI_COMM_WORLD);
        }
    } else {
        MPI_Comm_split(MPI_COMM_WORLD, 1, 0, &workers);
        // [0] is the no. of points inside the circle, [1] outside
        int results[2] = {0, 0};
        int global_results[2] = {0, 0};
        double coords[2];
        while (1) {
            MPI_Send(MPI_BOTTOM, 0, MPI_INT, 0, TAG_WORKER_READY, MPI_COMM_WORLD);
            MPI_Recv(&raw_coords[0], 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == TAG_EXIT)
                break;
            coords[0] = (double) raw_coords[0] / INT_MAX;
            coords[1] = (double) raw_coords[1] / INT_MAX;
            if ((coords[0] * coords[0] + coords[1] * coords[1]) <= 1)
                results[0]++;
            else
                results[1]++;
        }
        MPI_Allreduce(&results[0], &global_results[0], 2, MPI_INT, MPI_SUM, workers);
        double my_pi = ((double) results[0] / (results[1] + results[0])) * 4;
        double global_pi = ((double) global_results[0] / (global_results[1] + global_results[0])) * 4;
        printf("Process %d has estimated pi = %.15f based on local data, and pi = %.15f based on global data\n", rank, my_pi, global_pi);
    }

    finalize_exit(0);
}
