CFLAGS := -O3 -march=native
LDFLAGS := -lm

all: kd_tree_mpi kd_tree_openmp

kd_tree_mpi: CC := mpicc
kd_tree_mpi: CFLAGS += -DMPI
kd_tree_mpi: kd_tree.h

kd_tree_openmp: CFLAGS += -fopenmp
kd_tree_openmp: kd_tree.h

clean:
	rm kd_tree_mpi kd_tree_openmp
