#include <iostream>
#include "sequtial_run.h"

int main(){
	flann::Matrix<Scalar> *slice, *regions;
	// read a data slice 
	typeA * slice = readturbulence();

	// split the slice into regions
	regions = split_slice(slice);
	
	// caculate the divergence between all regions
	np_div()
	

	// do the clustering job
void getclustermedoids(int nclusters, int nelements, double** distance,
  int clusterid[], int centroids[], double errors[]);

	// may visulaization
	
	// free in this iteration
	free(regions);

}

