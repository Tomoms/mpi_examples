#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include "mpi.h"

#define SIZE    50

static inline int finalize_exit(int code)
{
    MPI_Finalize();
    exit(code);
}

int main(int argc, char *argv[])
{
    int rank, nproc;
    int vector[SIZE] = {0};

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    int limit = MIN(nproc - 1, SIZE);

    if (!rank) {
        srand(2021);
        int matrix[SIZE][SIZE] = {0};
        int *ptr = &matrix[0][0];
        int *stop = ptr + SIZE * SIZE;
        for (; ptr < stop; ptr++)
            *ptr = rand() % SIZE;
        
        printf("Manager: the matrix is:\n");
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                printf("%d\t", matrix[i][j]);
            }
            printf("\n");
        }
        printf("\n\n");

        ptr = &vector[0];
        stop = ptr + SIZE;
        for (; ptr < stop; ptr++)
            *ptr = rand() % SIZE;
        
        printf("Manager: the vector is:\n");
        for (int i = 0; i < SIZE; i++)
            printf("%d\t", vector[i]);
        printf("\n\n");

        MPI_Bcast(&vector[0], SIZE, MPI_INT, 0, MPI_COMM_WORLD);

        int sent = 0, received = 0, tmp;
        int result[SIZE], computers[SIZE];
        for (; sent < limit; sent++) {
            // sent is row number and also tag. row sent goes to process sent + 1
            MPI_Send(&matrix[sent], SIZE, MPI_INT, sent + 1, sent, MPI_COMM_WORLD);
        }
        while (received < SIZE) {
            MPI_Recv(&tmp, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            received++;
            // tag is row number
            int tag = status.MPI_TAG;
            result[tag] = tmp;
            computers[tag] = status.MPI_SOURCE;
            if (sent < SIZE) {
                MPI_Send(&matrix[sent], SIZE, MPI_INT, status.MPI_SOURCE, sent, MPI_COMM_WORLD);
                sent++;
            } else {
                // real tag is SIZE, which means "quit!"
                MPI_Send(MPI_BOTTOM, 0, MPI_INT, status.MPI_SOURCE, SIZE, MPI_COMM_WORLD);
            }
        }
        printf("Manager: the result is:\n");
        for (int i = 0; i < SIZE; i++)
            printf("%d\t", result[i]);
        printf("\n");
        printf("Manager: computed by:\n");
        for (int i = 0; i < SIZE; i++)
            printf("%d\t", computers[i]);
        printf("\n");
    } else {
        MPI_Bcast(&vector[0], SIZE, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank <= limit) {
            int row[SIZE];
            while (1) {
                MPI_Recv(&row[0], SIZE, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                if (status.MPI_TAG == SIZE)
                    finalize_exit(0);
                int dot_product = 0;
                for (int i = 0; i < SIZE; i++)
                    dot_product += row[i] * vector[i];
                MPI_Send(&dot_product, 1, MPI_INT, 0, status.MPI_TAG, MPI_COMM_WORLD);
            }
        }
    }

    finalize_exit(0);
}
