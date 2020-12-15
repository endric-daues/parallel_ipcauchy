
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "helpers.h"
#include "integrand.h"
#include "shortest_path.h"
#include <time.h>

#include "mpi.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// mpicc -Wall -o sp_main_mpi sp_main_mpi.c shortest_path.c helpers.c integrand.c -I/moto/home/ed2691/gsl/include -L/moto/home/ed2691/gsl/bin -std=c99 -Wall -lgsl -lgslcblas -lm


int main(int argc, char *argv[])
{
    double cpu_time_used,wall_time;
    clock_t start, end,start_wall,end_wall;
    start_wall = clock();

    double alpha = 0.0;
    double beta = 2*M_PI;
    double integration_error = 0;
    int max_num_points = 1000000;

    int N = 50;
    int a[50] = {485, 326, 248, 421, 322, 795, 43, 845, 955, 252, 9, 901, 122, 94, 738, 574, 715, 882, 367, 984, 299, 433, 682, 72, 874, 138, 856, 145, 995, 529, 199, 277, 97, 719, 242, 107, 122, 70, 98, 600, 645, 267, 972, 895, 213, 748, 487, 923, 29, 674};
    int b = 995;

    // Additional Variables needed for MPI

    int num_procs, rank, num_radii_proc, num_angles_proc;
    double segment_length, start_angle, end_angle, total, total_error,total_time;
    

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double min_radius = optimal_radius(0,a,b,N);

    int num_radii = 100;
    int num_angles = 360;

    segment_length = (beta - alpha) / (double) num_procs;

    if (rank == 0)
        {
            printf("Number of processes (intervals) = %d\n", num_procs);
            printf("Interval Length: %f\n",segment_length);
            printf("Minimum Radius: %f\n",min_radius);
            printf("Total Number of nodes: %d\n",(num_radii+1)*(num_angles+1));
        }
    

    MPI_Barrier(MPI_COMM_WORLD);

    start_angle = (double) rank * segment_length;
    end_angle = ((double) rank + 1.0) * segment_length;

    num_radii_proc = num_radii;
    num_angles_proc = round(num_angles / num_procs) + 1;

    printf("Process %d is starting with %d nodes... start angle %f and end angle %f\n", rank,(num_radii_proc+1)*(num_angles_proc+1),start_angle,end_angle);

    double ** master_radius_intervals = (double**)malloc((num_procs+1) * sizeof(double *));
    double ** master_angle_intervals = (double**)malloc((num_procs+1) * sizeof(double *));

    master_radius_intervals[rank] = (double*)malloc(sizeof(double) * ((num_radii_proc + 1) * (num_angles_proc + 1)));
    master_angle_intervals[rank] = (double*)malloc(sizeof(double) * ((num_radii_proc + 1) * (num_angles_proc + 1)));

    start = clock();
    
    int length = shortest_path(a,b,N,start_angle,end_angle,min_radius,num_radii_proc,num_angles_proc,master_radius_intervals[rank],master_angle_intervals[rank]);

    double total_result_proc = 0.0;
    double total_error_proc = 0.0;

    for (int i = 0; i <= length; i++){
        double result = 0.0;
        double error = 0.0;
        int intervals = 0.0;
        //printf("Interval %d, from radius: %f at angle: %f to radius: %f at angle: %f\n",i,master_radius_intervals[rank][i],master_angle_intervals[rank][i],master_radius_intervals[rank][i+1],master_angle_intervals[rank][i+1]);

        if (master_radius_intervals[rank][i] == master_radius_intervals[rank][i+1]) {
            //printf("Interval %d, from angle: %f to angle: %f at radius: %f\n",i,angle_intervals[i],angle_intervals[i+1],radius_intervals[i]);
            integrate_path_circle(a,b,master_radius_intervals[rank][i],master_angle_intervals[rank][i],master_angle_intervals[rank][i+1],N,max_num_points,integration_error, &result, &error, &intervals);
            total_result_proc += result;
            total_error_proc += error;
        }

        else if (master_angle_intervals[rank][i] == master_angle_intervals[rank][i+1]) {
            //printf("Interval %d, from radius: %f to radius: %f at angle: %f\n",i,radius_intervals[i],radius_intervals[i+1],angle_intervals[i]);
            integrate_path_circle(a,b,master_radius_intervals[rank][i],master_angle_intervals[rank][i],master_angle_intervals[rank][i+1],N,max_num_points,integration_error, &result, &error, &intervals);
            total_result_proc += result;
            total_error_proc += error;
        }
    }

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    
    MPI_Allreduce(&total_result_proc, &total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&total_error_proc, &total_error, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&cpu_time_used, &total_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD);

    // Wait for all processes to finish before outputting Info
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Process: %d, Computation Time: %f, result: %f\n", rank,cpu_time_used, total_result_proc);

    free(master_radius_intervals[rank]);
    free(master_angle_intervals[rank]);
    free(master_radius_intervals);
    free(master_angle_intervals);

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