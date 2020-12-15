complex double gamma_(double t, double R);
complex double H(int *a,double t, double R, complex double (*gamma_)(double t, double R),int N);
double integrand(int a[],int b,double t, double R,int N);
double integrand_abs(int a[],int b,double t, double R,int N);
double integrand_circle_gsl(double t, void * p);
double integrand_line_gsl(double R, void * p);