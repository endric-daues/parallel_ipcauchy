#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "integrand.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


/*General gamma function*/
complex double gamma_(double t, double R)
{
    complex double val = (complex double) 0.0;
    val = (complex double) R * cexp(t * I);
    return val;
}

/*Returns generating function value. Takes in the a vector, the current position 
    on the parametrized variable t and the gamma function, as well as the length of the vector. */
complex double H(int *a,double t, double R, complex double (*gamma_)(double t, double R),int N)
{
    complex double val = (complex double) 1.0;

    for (int i = 0; i < N; i++)
    {
        val =  val * 1 / (complex double) (1 - (cpow((*gamma_)(t,R),*(a+i))));
        
    }
    
    return val;
}

/*Returns the integrand value for a circular path. Takes in the a vector, the b value, the current position 
    on the parametrized variable t and the gamma function, as well as the length of the vector. */
double integrand(int a[],int b,double t, double R,int N)
{
    complex double val = 0.0;
    complex double h;

    h = H(a,t,R,gamma_,N);
    val = (1/(2 * M_PI)) * (h / cpow(gamma_(t,R),b));
    return val;

}

/* Returns the absolute value of the integrand */
double integrand_abs(int a[],int b,double t, double R,int N)
{
    complex double val = 0.0;
    complex double h;

    h = H(a,t,R,gamma_,N);
    val = (1/(2 * M_PI)) * (h / cpow(gamma_(t,R),b));

    //printf("val: %f\n",cabs(val));
    return cabs(val);
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

/*Returns the integrand value for a circular path. Takes in the a vector, the b value, the current position 
    on the parametrized variable t and the gamma function, as well as the length of the vector. */
double integrand_circle_gsl(double t, void * p)
{
    struct integrand_params * params = (struct integrand_params *)p;
    double R = (params->R);
    int N = (params->N);
    int b = (params->b);
    int *a = (params->a);
    

    complex double val = 0.0;
    complex double h;

    h = H(a,t,R,gamma_,N);
    val = (1/(2 * M_PI)) * (h / cpow(gamma_(t,R),b));
    return val;
}

/*Returns the integrand value for a linear path. Takes in the a vector, the b value, the current position 
    on the parametrized variable t and the gamma function, as well as the length of the vector. */
double integrand_line_gsl(double R, void * p)
{
    struct integrand_params_line * params = (struct integrand_params_line *)p;
    double t = (params->t);
    int N = (params->N);
    int b = (params->b);
    int *a = (params->a);
    

    complex double val = 0.0;
    complex double h;

    h = H(a,t,R,gamma_,N);
    val = (1/(2 * M_PI)) * (h / cpow(gamma_(t,R),b));
    return val;
}