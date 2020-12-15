#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "helpers.h"
#include "integrand.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/* Returns the Euclidean distance between two compelx poitns */
double dist(double r1, double t1, double r2, double t2)
{
    return sqrt(pow((r2*sin(t2) - r1*sin(t1)),2) + pow(r2*cos(t2) - r1*cos(t1),2));
}



/*Determines the cost of the edge between points p1 and p2 given a,b using the three point rule*/
double edge_value(double R1,double R2,double t1, double t2, int* a, int b, int N)
{   
    double t3;
    if (t1==t2){

        return 0.5 * fabs(integrand_abs(a,b,t1,R1,N) + integrand_abs(a,b,t1,R2,N)) * dist(R1,t1,R2,t2);
    }

    if (R1==R2) {

        t3 = 0.5*(t1+t2);
        return fabs(0.5*(integrand_abs(a,b,t1,R1,N)+integrand_abs(a,b,t3,R1,N)) * dist(R1,t1,R1,t3)) + fabs(0.5*(integrand_abs(a,b,t3,R1,N)+integrand_abs(a,b,t2,R1,N)) * dist(R1,t3,R1,t2));
    }    

    else {
        return INFINITY;
    }
}

/* Generates the graph to be used in the shortest path algorithm*/
void generate_graph(int num_radii, int num_angles,int nodes[num_radii+1][num_angles+1])
{   
    int counter = 0;

    for (int i = 0; i <= num_radii; i++){
        for (int j = 0; j <= num_angles; j++){
            nodes[i][j] = counter;
            counter += 1;
        }
    }
}


/* Generates the cost matrix for the shortest path problem */
void generate_cost_matrix(int num_radii, double radii[num_radii],int num_angles, double angles[num_angles],int nodes[num_radii+1][num_angles+1], int *a, int b, int N, int MAX,double C[MAX][MAX])
{
    double angle, radius, next_angle, next_radius, next_radius_in;
    int node, next_node_ahead,next_node_out, next_node_in;
    
    for (int i = 0; i < MAX; i++){
        for (int j = 0; j < MAX; j++){
            C[i][j] = INFINITY;
        }
    }

    for (int rad_index = 0; rad_index <= num_radii; rad_index++){
        for (int ang_index = 0; ang_index <= num_angles; ang_index++){

            radius = radii[rad_index];
            angle = angles[ang_index];
            node = nodes[rad_index][ang_index];
            // printf("rad_index: %d, radius: %f, ang_index: %d, angle: %f, node: %d\n",rad_index, radius,ang_index,angle,node);

            next_radius = radii[rad_index+1];
            next_angle = angles[ang_index+1];

            if (rad_index > 0){
                next_radius_in = radii[rad_index-1];
            }
            

            if ( radius < next_radius && (next_radius < 1) && (1e-12 < angle) && (ang_index < num_angles) ){
                
                next_node_out = nodes[rad_index+1][ang_index];
                // printf("Node: %d, Next node out: %d\n",node,next_node_out);
                C[node][next_node_out] = edge_value(radius,next_radius,angle,angle,a,b,N);
            
            }

            if (angle < next_angle && next_angle < 2 * M_PI && (ang_index < num_angles) ){
                
                next_node_ahead = nodes[rad_index][ang_index+1];
                // printf("Node: %d, Next node ahead: %d\n",node,next_node_ahead);
                C[node][next_node_ahead] = edge_value(radius,radius,angle,next_angle,a,b,N);
                
            }
            
            if ((rad_index > 0) && (1e-12 < angle) && (rad_index > 0) && (ang_index < num_angles)) {
                next_node_in = nodes[rad_index-1][ang_index];
                // printf("Node: %d, Next node in: %d\n",node,next_node_in);
                C[node][next_node_in] = edge_value(next_radius_in,radius,angle,angle,a,b,N); 
            }         
        }
    }
}

int printSolution(int MAX, double dist[MAX]) 
{ 
    printf("Vertex   Distance from Source\n"); 
    for (int i = 0; i < MAX; i++) 
        printf("%d        %f\n", i, dist[i]); 
    
    return 0;
} 

int minDistance(int MAX, double dist[MAX], int spSet[MAX]) 
{ 
    // Initialize min value 
    double min = INFINITY;
    int min_index; 
  
    for (int v = 0; v < MAX; v++) {
        if (spSet[v] == 0 && dist[v] < min) {
            min = dist[v];
            min_index = v; 
        }
    }  
    return min_index; 
} 

// Implementation of Dijkstras shortest path algorithm //
int dijkstra(int MAX, double C[MAX][MAX],int source_node,int des, int path[MAX])
{
    double key[MAX];
    int spSet[MAX];
    int parent[MAX];

    for (int i = 0; i < MAX; i++) {
        key[i] = INFINITY;
        spSet[i] = 0;
    }
    
    key[source_node] = 0;
    parent[source_node] =-1;
   

    for (int count = 0; count < MAX - 1; count++) {

        int u = minDistance(MAX, key, spSet);
        spSet[u] = 1; 

        for (int v = 0; v < MAX; v++){
            if (spSet[v]==0 && C[u][v] < INFINITY) {
                parent[v] = u;
                key[v] = C[u][v]; 
            }
        } 
    
    }
    int node = des;
    int count = 0;
    while (node > -1){
        node = parent[node];
        count+=1;
    }

    int temp[count];
    node = des;

    for(int i=0; i<count;i++){
        temp[i] = node;
        node = parent[node];
    }

    for(int i=0; i < count;i++){
        path[i] = temp[count-1 -i];
    
    }


   return count-1;
}

double get_radius(int node, int num_radii, double *radii, int num_angles, double *angles, int nodes[num_radii+1][num_angles+1]){
    for (int ang = 0; ang <= num_angles; ang++){
            for (int rad = 0; rad <= num_radii; rad++){            
                if (node == nodes[rad][ang]){
                    return radii[rad];
                }
            }
    }
    
    return 0.0;
}

double get_angle(int node, int num_radii, double *radii, int num_angles, double *angles,int nodes[num_radii+1][num_angles+1]){
    for (int ang = 0; ang <= num_angles; ang++){
            for (int rad = 0; rad <= num_radii; rad++){            
                if (node == nodes[rad][ang]){
                    return angles[ang];   
                }
            }
    }
    
    return 0.0;
}

int shortest_path(int *a, int b, int N, double alpha, double beta, double min_radius,int num_radii, int num_angles,double *radius_intervals, double *angle_intervals)
{
    int MAX = (num_radii + 1) * (num_angles + 2);

    
    int (*nodes)[num_radii+2];
    nodes = malloc(sizeof *nodes * (num_angles + 2));

    double *radii;
    double *angles;

    radii = (double*)malloc(sizeof(double) * (num_radii+2));
    angles = (double*)malloc(sizeof(double) * (num_angles+2));
    
    double (*C)[MAX];
    C = malloc(sizeof *C * MAX);

    generate_graph(num_radii, num_angles, nodes);

    double radius_interval = (0.999 - min_radius) / (double) num_radii;
    double angle_interval = (beta - 1e-12 - alpha ) / ((double) num_angles);

    for(int i = 0; i <= num_radii; i++){
        radii[i] = min_radius + (double) i * radius_interval;
    }
    
    for(int i = 0; i <=  num_angles; i++){
        angles[i] = alpha + (double) i * angle_interval;
    }

    generate_cost_matrix(num_radii,radii,num_angles,angles,nodes,a,b,N,MAX,C);

    double start_radius = optimal_radius(angles[0],a,b,N);

    int start_idx = 0;
    double min_dis = INFINITY;
    for (int i=0; i<num_radii;i++){
        if (min_dis > fabs(radii[i] - start_radius)){
            min_dis = fabs(radii[i] - start_radius);
            start_idx = i;
        }
    }

    int src = nodes[start_idx][0];
    double end_radius = optimal_radius(angles[num_angles+1],a,b,N);

    int end_idx = 0;
    min_dis = INFINITY;

    for (int i=0; i<num_radii;i++){
        if (min_dis > fabs(radii[i] - end_radius)){
            min_dis = fabs(radii[i] - end_radius);
            end_idx = i;
        }
    }


    int des = nodes[end_idx][num_angles];

    int path[MAX];
    int length;

    length = dijkstra(MAX,C,src,des,path);
    //printf("Length: %d\n", length);

    for (int i=0; i<= length ; i++){
        radius_intervals[i] = get_radius(path[i], num_radii, radii, num_angles,  angles, nodes);
        angle_intervals[i] = get_angle(path[i], num_radii, radii, num_angles,  angles, nodes);
        //printf("radius_interval[%d]: %f, angle_interval[%d]: %f\n",i,radius_intervals[i],i,angle_intervals[i]);
    }

    free(nodes);
    free(C);
    free(angles);
    free(radii);

return length;

}

