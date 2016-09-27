#include <stdio.h>
#include <string.h>
#include "divide.h"


extern int read_data(const char* file_name, float **pressure, float **velocity, int *dim1, int *dim2, int *dim3);
extern int free_data(float *pressure, float *velocity);
extern float get_divs(float *A, float *B, int region_length, int k, int div_func);

int main(){
    
    printf(" read data to buffer\n");
    //const char *file_name = "data/test_1_2_3_4.h5";
    const char *file_name = "data/isotropic_201_201_1.h5";
    float *pressure, *velocity;
    int dim1, dim2, dim3;

    int i,j,ii, jj, p, q, ret;


    // divide configrations, length of 10 region will have 11*11 121 points
    int region_length = 10;

    // how many clusters 
    int ncluster = 3;
    int npass = 1;

    // clustering result
    /*
    int clusterid[];
    int centroid;
    float error;
    int ifound;
    */

    // read data to buffer
    read_data(file_name, &pressure, &velocity, &dim1, &dim2, &dim3);
    //print_p_data(pressure, dim1, dim2, dim3);
    //save visualized result
    
    printf("data is read, dimension is %d * %d * %d\n", dim1, dim2, dim3);
    


	// split the slice into regions
    // since we know each region's exact size, we only need record the starting address
    int num_region;

    // num_regions* points_in_region*3
    float  *regions;

    // dim1 = dim2 = 201, in a region, 11 points in each side, thus length 10, the whole block has length 201-1 = 200. So (200/10)^2 400 regions
    // d1 d2 d3 is the 3-dimension matrix 
    int d1 = dim1;
    int d2 = region_length*region_length;;
    int d3 = 3;
    divide(velocity, d1, region_length, &num_region, &regions);

    // distance matrix
    /*
    float **distance = (float **)malloc(sizeof(float *)* num_region);
    if(distance == NULL){
        perror("allocate error 1");
        exit(-1);
    }
       
    for(i = 0; i < num_region; i++){
        distance[i] = (float *)malloc(sizeof(float) * num_region);
        if(distance[i] == NULL){
            perror("allocation error 2");
            exit(-1);
        }
    }
    */
    float distance[num_region][num_region];
	
	// caculate the divergence between all regions
    // distance matrix can be very large
    
    float  div;
    // find k=5 nearest neighbours to determine the density
    int k = 5;
    // use L-2 divergence
    int div_func = 1;
    char *dist_path = "data/all_dist.txt";

    FILE * f_dist = fopen(dist_path, "w");
    if(p_dist_path == NULL){
        perror("cannot access the dist_path file");
        exit(-1);
    }


    for(i = 0; i< num_region; i++){
        for(j = i; j < num_region ; j++){
            // starting address of each region as input
            div = get_divs( regions + i*d2*d3 , regions + j*d2*d3, region_length, k, div_func);
            printf("\t divergence between region %d and %d is %.3f\n", i, j, div);
            fprintf(f_dist,"%f\n",div);
            /*
            distance[i][j] = div;
            if(i != j){
                distance[j][i] = div;
            }
            */
        }
    }

    printf("distance matrix is saved in %s\n", dist_path);
    fclose(f_dist);

    // do clustering
    //    kmedoid(nclusters, nelements, distance, npass, clusterid, &error, &ifound);

    // save the results(nclusters)?
    // or visulize the clustering results?
    
	// if regions are just pointers, no need to free 
	free(regions);
    // free buffer
    free_data(pressure, velocity);
    // also free pdata block
    
    return 1;
}
