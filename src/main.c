#include <iostream>
#include "read_file.h"
#include "divide.h"
#include "get_div.h"


int main(){
	// read a data slice 
    char *file_name = "some turbulence datacut";

    int i,j,ret;

    // how many clusters 
    int numcluster;
    the results
    int clusterid[];

    // size of the turbulence cut
    int num_row, num_col;

    // flat representation of all velocity data
    double (*pdata)[3];

    

	ret = readturbulence(file_name, pdata,  double&num_row, &num_col);
    if(ret == 1){
        printf("velocity data is read from %s\n", file_name);
    }


	// split the slice into regions
    char * regions[];
    int num_region;

    double **distance = (double **)malloc(sizeof(double *)* num_region);
    check_malloc(distance);
    for(i = 0; i < num_region; i++){
        distance[i] = (double *)malloc(sizeof(double) * num_region);
        check_malloc(distance);
    }

	num_region = divide(pdata, num_col,regions);
	
	// caculate the divergence between all regions
    double div;
    for(i = 0; i< num_region; i++){
        for(j = i+1; j < num_region ; j++){
            div = get_divs()

            distance[i][j]

    get_div
	

	// do the clustering job
void getclustermedoids(int nclusters, int nelements, double** distance,
  int clusterid[], int centroids[], double errors[]);

	// may visulaization
	
	// free in this iteration
	free(regions);

}

