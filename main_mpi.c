#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "helpers.h"
#include <time.h>
#include "mpi.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// mpicc -Wall -o main_mpi main_mpi.c helpers.c integrand.c -I/moto/home/ed2691/gsl/include -L/moto/home/ed2691/gsl/bin -std=c99 -Wall -lgsl -lgslcblas -lm

int main(int argc, char *argv[])
{
    double r,cpu_time_used,result,error,wall_time;
    double total = 0.0;
    double total_error = 0.0;
    clock_t start, end, start_wall,end_wall;

    start_wall = clock();

    double integration_error = 0;
    int max_num_points = 1000000;
    

    // Additional Variables needed for MPI
    int num_procs, rank, intervals;
    double segment_length, start_angle, end_angle,total_time;
    

    // Instance Definition
    int N = 50;
    int a[50] = {485, 326, 248, 421, 322, 795, 43, 845, 955, 252, 9, 901, 122, 94, 738, 574, 715, 882, 367, 984, 299, 433, 682, 72, 874, 138, 856, 145, 995, 529, 199, 277, 97, 719, 242, 107, 122, 70, 98, 600, 645, 267, 972, 895, 213, 748, 487, 923, 29, 674};
    int b = 995;

    r = optimal_radius(0,a,b,N);
    
    double alpha = 0.0;
    double beta = 2*M_PI;
    
    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    segment_length = (beta - alpha) / (double) num_procs;

    if (rank == 0)
        {
            printf("Number of processes (intervals) = %d\n", num_procs);
            printf("Interval Length: %f\n",segment_length);
            printf("Optimal Radius: %f\n",r);
        }
    
    MPI_Barrier(MPI_COMM_WORLD);

    printf("Process %d is starting...\n", rank);
    
  
    result = 0.0;
    error = 0.0;
    intervals = 0;

    start_angle = (double) rank * segment_length;
    end_angle = ((double) rank + 1.0) * segment_length;

    start = clock();
    
    integrate_path_circle(a,b,r,start_angle,end_angle,N,max_num_points,integration_error, &result, &error, &intervals);

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    MPI_Allreduce(&result, &total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&error, &total_error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&cpu_time_used, &total_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD);

    // Wait for all processes to finish before outputting Info
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Process: %d, Computation Time: %f\n", rank,cpu_time_used);
    printf("Process %d used %d intervals\n\n",rank,intervals);


    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("Total Integration Result: %f\n", total);
        printf("Total Integration Error: %f\n", total_error);
        printf("Maximum Computation Time: %f\n\n", total_time);
    }

    MPI_Finalize();

    end_wall = clock();

    if (rank == 0) 
    {
    wall_time = ((double) (end_wall - start_wall)) / CLOCKS_PER_SEC;
    printf("Wall Time: %f\n\n", wall_time);
    }

    return 0;
}
