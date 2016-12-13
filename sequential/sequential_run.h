#ifndef SEQUENTIAL_RUN_H
#define SEQUENTIAL_RUN_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "region_def.h"
#include "divide.h"
#include "read_file.h"
#include "get_divs.h"
#include "common_utility.h"
#include "cluster.h"



// get current time
double get_cur_time();

/*
 * get all the divergence given all the regions
 * input:
 *  regions description:
 *      regions: start address of all the regions
 *      d2, d3: number of points in each region, number of data fields in each region
 *      region_length
 *  k_npdiv: the k in npdiv density estimation
 *  div_func: divergence functions
 *  
 * output:
 *  matrix: all divergence will saved into it
 *  divs_path: also divs will saved into this file(for debugging use)
 *
 * return: return 0 if success
 */
int get_all_divs(float *regions,int num_region, int d2, int d3, int region_length,double **matrix,char *divs_path, int k_npdiv, int div_func);



/*
 * run kmedoids method based on divergence matrix
 * input:
 *  ncluster, npass: see cluster/cluster.h
 *  num_region: number of regions
 *  matrix: all the divergence
 * output
 *  error, ifound: see cluster/cluster.h
 *  output_path: clusteringid
 */
int run_clustering (int nclusters, int num_region, double **matrix,
  int npass, char *output_path, double* error, int* ifound);
#endif
