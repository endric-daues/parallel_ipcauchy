double optimal_radius(double angle,int a[],int b,int N);
int integrate_path_circle(int *a,int b,double r,double alpha, double beta,int N,int max_num_points,int integration_error, double *result, double *error, int *intervals);
int integrate_path_line(int *a,int b,double r1, double r2,double t,int N,int max_num_points,int integration_error, double *result, double *error, int *intervals);
int integrate_path_circle_double_adaptive(int *a,int b,double r,double alpha, double beta,int N,int max_num_points,int integration_error, double result, double error);
int integrate_path_circle_romberg(int *a,int b,double r,double alpha, double beta,int N,int max_num_points,int integration_error, double result, double error);
int read_pisinger_instance(char file_name[],int a[],int *b);