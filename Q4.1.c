#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int count = 0;
    double x, y, z, pi;
    int rank, size, iter_per_proc;
    double start_time, end_time;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    iter_per_proc = NUM_ITER / size;
    
    srand(SEED * rank);
    
    start_time = MPI_Wtime();
    
    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < iter_per_proc; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        // Check if point is in unit circle
        if (z <= 1.0)
        {
            count++;
        }
    }
    
    int total_count;
    MPI_Reduce(&count, &total_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    end_time = MPI_Wtime();
    
    if (rank == 0) {
        // Estimate Pi and display the result
        pi = ((double)total_count / (double)NUM_ITER) * 4.0;
        printf("The result is %f\n", pi);
        
        printf("Execution time: %f seconds\n", end_time - start_time);
    }
    
    MPI_Finalize();
    return 0;
}
