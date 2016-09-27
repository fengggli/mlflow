/*
 * Sep 10, by Feng Li, IUPUI
 * based on the divergences beteen the regions, do the clustering job
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../cluster/cluster.h"

/*
 *  kmedoid methods is implemented in cluster.h
 *
 */

void run_clustering (int nclusters, int num_region, char *dist_path,
  int npass, char *output_path, double* error, int* ifound){
    int i, j;
    double div;
    char str[20];
    int clusterid[num_region];

    // construct the distance matrix
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
    
    //double distance[num_region][num_region];

    // open the divergence file
    FILE * p_dist_path = fopen(dist_path, "r");
    if(p_dist_path == NULL){
        perror("cannot access the dist_path file");
        exit(-1);
    }

    for(i = 0; i< num_region; i++){
        for(j = i; j < num_region ; j++){
            // starting address of each region as input
            //div = get_divs( regions + i*d2*d3 , regions + j*d2*d3, region_length, k, div_func);
            //printf("\t divergence between region %d and %d is %.3f\n", i, j, div);
            fgets(str,60,p_dist_path);
            div = atof(str);
            //printf("\t divergence between region %d and %d is %.3lf\n", i, j, div);
            distance[i][j] = div;
            if(i != j){
                distance[j][i] = div;
            }
        }
    }
    fclose(p_dist_path);
    // reconstruct the distmatrix
    printf("***starting to run kmedoids algorithm\n");
    printf("nclusters = %d\nnum of elements %d\n npass %d\n. error = %.3lf\n", nclusters,num_region,npass,*error);
    
    kmedoids (nclusters, num_region, distance, npass,clusterid,error,ifound);
    
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

