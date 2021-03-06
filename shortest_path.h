double dist(complex double p1, complex double p2);
double edge_value(double R1,double R2,double t1, double t2, int* a, int b, int N);
void generate_graph(int num_radii, int num_angles,int nodes[num_radii+1][num_angles+1]);
void generate_cost_matrix(int num_radii, double radii[num_radii],int num_angles, double angles[num_angles],int nodes[num_radii+1][num_angles+1], int *a, int b, int N, int MAX,double C[MAX][MAX],double alpha,double beta);
int printSolution(int MAX, double dist[MAX]);
int minDistance(int MAX, double dist[MAX], int spSet[MAX]);
int dijkstra(int MAX, double C[MAX][MAX],int source_node,int des, int path[MAX]);
int shortest_path(int *a, int b, int N, double alpha, double beta, double min_radius,int num_radii, int num_angles,double *radius_intervals, double *angle_intervals);
double get_angle(int node, int num_radii, double *radii, int num_angles, double *angles,int nodes[num_radii+1][num_angles+1]);
double get_radius(int node, int num_radii, double *radii, int num_angles, double *angles,int nodes[num_radii+1][num_angles+1]);