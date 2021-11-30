#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
    int rank, nproc;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status_add, status_subtract;

    const int right = (rank + 1) % nproc;
    const int left = !rank ? nproc - 1 : rank - 1;
    int payload_add = 0, payload_subtract = 0;
    int received = 0;

    if (!rank) {
        int former_add, former_subtract;
        for (int i = 0; i < nproc; i++) {
            former_add = payload_add;
            MPI_Recv(&payload_add, 1, MPI_INT, right, 10 * right, MPI_COMM_WORLD, &status_add);
            received++;
            MPI_Send(&former_add, 1, MPI_INT, left, 10 * rank, MPI_COMM_WORLD);

            former_subtract = payload_subtract;
            MPI_Recv(&payload_subtract, 1, MPI_INT, left, 10 * left, MPI_COMM_WORLD, &status_subtract);
            received++;
            MPI_Send(&former_subtract, 1, MPI_INT, right, 10 * rank, MPI_COMM_WORLD);
        }
    } else {
        for (int i = 0; i < nproc; i++) {
            payload_add += rank;
            MPI_Send(&payload_add, 1, MPI_INT, left, 10 * rank, MPI_COMM_WORLD);
            MPI_Recv(&payload_add, 1, MPI_INT, right, 10 * right, MPI_COMM_WORLD, &status_add);
            received++;

            payload_subtract -= rank;
            MPI_Send(&payload_subtract, 1, MPI_INT, right, 10 * rank, MPI_COMM_WORLD);
            MPI_Recv(&payload_subtract, 1, MPI_INT, left, 10 * left, MPI_COMM_WORLD, &status_subtract);
            received++;
        }
    }
    
    printf("I am process %d and I have received %d messages.", rank, received);
    printf("My final messages have got tag %d (from the right) and %d (from the left), and values %d (from the right) and %d (from the left)\n",
           status_add.MPI_TAG, status_subtract.MPI_TAG, payload_add, payload_subtract);

    MPI_Finalize();
    return 0;
}
