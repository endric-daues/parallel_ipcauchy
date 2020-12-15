#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <complex.h>
#include "helpers.h"
#include "integrand.h"
#include <gsl/gsl_integration.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#ifndef INT64_MAX
	#define INT64_MAX 2147483647
#endif
 
/*Determines the optimal bypass point at any given angle, thus calculating
    the optimal radius for bypassing the ray at the given angle. Takes
    in the angle, the coefficient vector, the b value and the size of 
    the vector, N.*/
double optimal_radius(double angle,int a[],int b, int N)
{
    double r = 0.0;
    double current_val = 0.0;
    double lower_bound = 0.5;
    double upper_bound = 0.999999;
    int num_points = 1000;
    double current_min = INT64_MAX;
    double max_r = lower_bound;


    for (int i=0; i <= num_points; i++)
    {
        r = lower_bound + (double) i * ((upper_bound - lower_bound) / (double) num_points);
        current_val = fabs(integrand(a,b,angle,r,N));

        if (current_val < current_min)
        {
            current_min = current_val;
            max_r = r;
        }
        

    }

    return max_r;
}

/*Struct used to pass integral parameters to the integration function */
struct integrand_params {
    int b;
    double R;
    int N;
    int *a;
};

/*Struct used to pass integral parameters to the integration function for a straight line integral with a constant angle */
struct integrand_params_line {
    int b;
    double t;
    int N;
    int *a;
};

/*Integration caller function for implementing the path integral on a circular path on a given radius. Takes
    the coefficient vector, the b value, the radius and the angles between which to integrate, 
    as well as various variables to complete the integration such as the integration error,
    and pointers to the result and error variables.*/
int integrate_path_circle(int *a,int b,double r,double alpha, double beta,int N,int max_num_points,int integration_error, double *result, double *error, int *intervals)
{
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (max_num_points);

    struct integrand_params sigma;
    sigma.R = r;
    sigma.b = b;
    sigma.N = N;
    sigma.a = a;

    gsl_function F;
    F.function = &integrand_circle_gsl;
    F.params = &sigma;

    gsl_integration_qags(&F, alpha, beta, integration_error, 1e-8, max_num_points,w, result, error);

    int temp = (int) w->size;
    *intervals = temp;
    

    gsl_integration_workspace_free(w);

    return 0;
}

/*Integration caller function for implementing the path integral on a line with constant angle*/
int integrate_path_line(int *a,int b,double r1, double r2,double t,int N,int max_num_points,int integration_error, double *result, double *error, int *intervals)
{
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (max_num_points);

    struct integrand_params_line sigma;
    sigma.t = t;
    sigma.b = b;
    sigma.N = N;
    sigma.a = a;

    gsl_function F;
    F.function = &integrand_line_gsl;
    F.params = &sigma;

    gsl_integration_qags(&F, r1, r2, integration_error, 1e-8, max_num_points,w, result, error);

    int temp = (int) w->size;
    *intervals = temp;
    
    gsl_integration_workspace_free(w);

    return 0;
}

