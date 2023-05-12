#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int count = 0, rank, size, tag = 123, prev, next;
    double x, y, z, pi;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int iterations_per_process = NUM_ITER / size;

    srand(SEED * rank);

    for (int iter = 0; iter < iterations_per_process; iter++)
    {
        x = (double)rand() / (double)RAND_MAX;
        y = (double)rand() / (double)RAND_MAX;
        z = sqrt((x * x) + (y * y));

        if (z <= 1.0)
        {
            count++;
        }
    }

    if (rank == 0)
    {
        double start = MPI_Wtime();
        for (int i = 1; i < size; i++)
        {
            MPI_Recv(&count, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            pi += (double)count / (double)iterations_per_process;
        }
        pi += (double)count / (double)iterations_per_process;
        pi *= 4.0;
        double end = MPI_Wtime();
        printf("The result is %f\n", pi);
        printf("Time elapsed is %f seconds\n", end - start);
    }
    else
    {
        int dest = 0;
        MPI_Send(&count, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
