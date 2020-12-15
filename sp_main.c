
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "helpers.h"
#include "integrand.h"
#include "shortest_path.h"
#include <time.h>


int main(int argc, char *argv[])
{
    int N = 50;
    int a[50] = {485, 326, 248, 421, 322, 795, 43, 845, 955, 252, 9, 901, 122, 94, 738, 574, 715, 882, 367, 984, 299, 433, 682, 72, 874, 138, 856, 145, 995, 529, 199, 277, 97, 719, 242, 107, 122, 70, 98, 600, 645, 267, 972, 895, 213, 748, 487, 923, 29, 674};
    int b = 995;



    double total_result = 0.0;
    double total_error = 0.0;
    clock_t start, end;
    double cpu_time_used;

    int num_radii = 100;
    int num_angles = 360; 
    
    double integration_error = 0;
    int max_num_points = 1000000;

    double min_radius = optimal_radius(0,a,b,N);
    printf("Min Radius: %f\n",min_radius);

    double * radius_intervals;
    double * angle_intervals;

    radius_intervals = (double*)malloc(sizeof(double) * ((num_radii + 1) * (num_angles + 1)));
    angle_intervals = (double*)malloc(sizeof(double) * ((num_radii + 1) * (num_angles + 1)));
    int length;

    start = clock();


    length = shortest_path(a,b,N,0,2*M_PI,min_radius,num_radii,num_angles,radius_intervals,angle_intervals);


    for (int i = 0; i < length; i++){
        double result = 0.0;
        double error = 0.0;
        int intervals = 0.0;
        printf("Interval %d, from radius: %f at angle: %f to radius: %f at angle: %f\n",i,radius_intervals[i],angle_intervals[i],radius_intervals[i+1],angle_intervals[i+1]);

        if (radius_intervals[i] == radius_intervals[i+1]) {
            printf("Interval %d, from angle: %f to angle: %f at radius: %f\n",i,angle_intervals[i],angle_intervals[i+1],radius_intervals[i]);
            integrate_path_circle(a,b,radius_intervals[i],angle_intervals[i],angle_intervals[i+1],N,max_num_points,integration_error, &result, &error, &intervals);
            total_result += (result);
            total_error += error;
            
        }

        else if (angle_intervals[i] == angle_intervals[i+1]) {
            printf("Interval %d, from radius: %f to radius: %f at angle: %f\n",i,radius_intervals[i],radius_intervals[i+1],angle_intervals[i]);
            integrate_path_line(a,b,radius_intervals[i],radius_intervals[i+1],angle_intervals[i],N,max_num_points,integration_error, &result, &error, &intervals);
            total_result += (result);
            total_error += error;
        }
    }

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    
    free(radius_intervals);
    free(angle_intervals);

    printf("Total Integral Value: %f\n", total_result);
    printf("Total Error Value: %f\n", total_error);
    printf("Total Computation Time: %f\n",cpu_time_used);
    return 0;
}