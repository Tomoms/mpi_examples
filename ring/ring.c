#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
	int rank, nproc;
	const int one[1] = {1};
	int neighbors[2];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm ring_comm;
	MPI_Status status_add, status_subtract;

	int ret = MPI_Cart_create(MPI_COMM_WORLD, 1, &nproc, one, 1, &ring_comm);
	if (ret) {
		perror("MPI_Cart_create()");
		MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	}

	MPI_Cart_shift(ring_comm, 0, 1, &neighbors[0], &neighbors[1]);
	MPI_Comm_rank(ring_comm, &rank);
	int payload_add = 0, payload_subtract = 0;
	int received = 0;

	if (!rank) {
		int former_add, former_subtract;
		for (int i = 0; i < nproc; i++) {
			former_add = payload_add;
			MPI_Recv(&payload_add, 1, MPI_INT, neighbors[1], 10 * neighbors[1], ring_comm, &status_add);
			received++;
			MPI_Send(&former_add, 1, MPI_INT, neighbors[0], 10 * rank, ring_comm);

			former_subtract = payload_subtract;
			MPI_Recv(&payload_subtract, 1, MPI_INT, neighbors[0], 10 * neighbors[0], ring_comm, &status_subtract);
			received++;
			MPI_Send(&former_subtract, 1, MPI_INT, neighbors[1], 10 * rank, ring_comm);
		}
	} else {
		for (int i = 0; i < nproc; i++) {
			payload_add += rank;
			MPI_Send(&payload_add, 1, MPI_INT, neighbors[0], 10 * rank, ring_comm);
			MPI_Recv(&payload_add, 1, MPI_INT, neighbors[1], 10 * neighbors[1], ring_comm, &status_add);
			received++;

			payload_subtract -= rank;
			MPI_Send(&payload_subtract, 1, MPI_INT, neighbors[1], 10 * rank, ring_comm);
			MPI_Recv(&payload_subtract, 1, MPI_INT, neighbors[0], 10 * neighbors[0], ring_comm, &status_subtract);
			received++;
		}
	}

	printf("I am process %d and I have received %d messages. ", rank, received);
	printf("My final messages have got tag %d (from the right) and %d (from the left), and values %d (from the right) and %d (from the left)\n",
		status_add.MPI_TAG, status_subtract.MPI_TAG, payload_add, payload_subtract);

	MPI_Finalize();
	return 0;
}
