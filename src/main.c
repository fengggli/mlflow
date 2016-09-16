#include <iostream>
#include "read_file.h"
#include "divide.h"
#include "get_div.h"


int main(){
	// read a data slice 
    char *file_name = "some turbulence datacut";

    int i,j,ret;


    // divide configrations
    int region_length = 11;

    // how many clusters 
    int ncluster = 3;
    int npass = 1;

    // clustering result
    int clusterid[];
    int centroid;
    double error;
    int ifound;

    // size of the turbulence cut
    int num_row, num_col;

    // flat representation of all velocity data
    double (*pdata)[3];

	ret = read_data(file_name, pdata,  double&num_row, &num_col);
    if(ret == 1){
        printf("velocity data is read from %s\n", file_name);
    }

	// split the slice into regions
    // since we know each region's exact size, we only need record the starting address
    double (* regions[])[3];
    int num_region;

    // distance matrix
    double **distance = (double **)malloc(sizeof(double *)* num_region);
    check_malloc(distance);
    for(i = 0; i < num_region; i++){
        distance[i] = (double *)malloc(sizeof(double) * num_region);
        check_malloc(distance);
    }

	num_region = divide(pdata, num_col,regions, region_length);
	
	// caculate the divergence between all regions
    // distance matrix can be very large
    
    double div;
    for(i = 0; i< num_region; i++){
        for(j = i+1; j < num_region ; j++){
            div = get_divs( regions[i], regions[j], region_length);

            distance[i][j] = div;
            dstance[j][i] = div;
        }

    // do clustering
    kmedoid(nclusters, nelements, distance, npass, clusterid, &error, &ifound);

    // save the results(nclusters)?
    // or visulize the clustering results?
    
	// if regions are just pointers, no need to free 
	free(regions);
    // also free pdata block

}

