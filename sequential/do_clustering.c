#include "do_clustering.h"
#include "sequential_run.h"
#include <stdlib.h>

void run_clustering (int nclusters, int num_region, char *dist_path,
  int npass, char *output_path, double* error, int* ifound){
    int i, j;
    double div;
    char str[20];
    int clusterid[num_region];

    // construct the distance matrix
    //
    double **matrix = malloc(num_region*sizeof(double *)); 
            if(matrix == NULL){
                perror("malloc for div matrix");
                exit(-1);
            }
            matrix[0] = NULL;
            for(i = 1; i< num_region; i++){
                matrix[i] = malloc(i*sizeof(double));
                if(matrix[i] == NULL) break;
            }
            if(i< num_region){
                for(j = 0 ;j< i; j++)
                    free(matrix[j]);
                perror("malloc 2 for div matrix");
                exit(-2);
            }
     /*
    double** distance = (double **)malloc(sizeof(double *)*num_region);
    if(distance == NULL){
        perror("malloc");
        exit(-1);
    }
    for(i = 0; i< num_region; i++){
        distance[i] = (double *)malloc(sizeof(double)*num_region);
        if(distance[i] == NULL){
            perror("malloc");
            exit(-1);
        }
    }
    */
    
    //double distance[num_region][num_region];

    // open the divergence file
    FILE * p_dist_path = fopen(dist_path, "r");
    if(p_dist_path == NULL){
        perror("cannot access the dist_path file");
        exit(-1);
    }

    for(i = 1; i< num_region; i++){
        for(j = 0; j < i ; j++){
            // starting address of each region as input
            //div = get_divs( regions + i*d2*d3 , regions + j*d2*d3, region_length, k, div_func);
            //printf("\t divergence between region %d and %d is %.3f\n", i, j, div);
            fgets(str,60,p_dist_path);
            div = atof(str);
            //printf("\t divergence between region %d and %d is %.3lf\n", i, j, div);
            matrix[i][j] = div;
            //printf("read div between region No.%d and No.%d: %lf", i, j, div);
        }
    }
    fclose(p_dist_path);
    // reconstruct the distmatrix
    printf("***starting to run kmedoids algorithm\n");
    printf("nclusters = %d\nnum of elements %d\n npass %d\n. error = %.3lf\n", nclusters,num_region,npass,*error);
    
    kmedoids (nclusters, num_region, matrix, npass,clusterid,error,ifound);

    // free the divergence buffer
    for(i = 1; i< num_region; i++){
        if(matrix[i] != NULL) free(matrix[i]);
    }
    if(matrix != NULL) {
        free(matrix);
        printf("distance matrix freed\n");
    }
    
    // save cluster results into file
    FILE * f_clusterid = fopen(output_path, "w");
    if(f_clusterid == NULL){
        perror("file open error");
        exit(-1);
    }

    for(i = 0; i < num_region; i++){
        fprintf(f_clusterid, "%d\n",clusterid[i]);
    }
    
    fclose(f_clusterid);
}



int main(){
    // how many neighbors we need to look at in divs
    int k_npdiv = K_NPDIV;
    // what is the current timestep
    int timestep = 0;

    char dist_path[STRING_LENGTH]; 
    char output_path[STRING_LENGTH]; 
    //sprintf(dist_path, "data/sequential/all_dist_201_k_%d_t_%d.txt", k_npdiv, timestep);
    //sprintf(output_path, "data/sequential/clusterid_201_k_%d_t_%d.txt", k_npdiv, timestep);
    sprintf(dist_path, "data/sequential/all_dist_201_k_%d_t_%d_ragged.txt", k_npdiv, timestep);
    sprintf(output_path, "data/sequential/clusterid_201_k_%d_t_%d_ragged.txt", k_npdiv, timestep);
    /*
    char *dist_path = "data/all_dist_201_lg.txt";
    char *output_path = "data/clusterid_201_lg_1.txt";
    */
    int nclusters, num_region,npass, ifound;
    double error;

    // timer
    double time_start, time_current;

    nclusters = NCLUSTERS;
    num_region =NUM_REGION ; 
    npass =NPASS;

    //int clusterid[num_region];


    time_start = get_cur_time();

    //call the divergence function
    printf("read divergence from %s, num_regions is %d\n", dist_path, num_region );
    run_clustering (nclusters,num_region, dist_path,npass,output_path, &error,&ifound);

    time_current = get_cur_time();
    printf("clustering completed, result  is stored in %s\n", output_path);
    printf("%.4f time is used for clustering\n", time_current - time_start );
    printf("error is %.3lf, num of runs is %d\n", error, ifound);

    return 0;
}
