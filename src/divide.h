#ifndef DIVIDE_H 
#define DIVIDE_H 

#include <stdio.h>
#include <stdlib.h>
/*
 * Sep 7, by Feng Li, IUPUI
 * Devide a slice of turbulence data into regions
 * before splitting, we have a block of a*b*1, default value of c is 1
 * each region is a block with a length of l, ovelapping length default value is 0
 */


/*
 * input
 *  pdata, pointer to data buffer
 *  dim: num of points in x and y dimension 
 *  l: region_length: for 11 points in each size, the length will be 10
 * output:
 *  p_num_region; number of regions
 *  p_regions: the tripple features(vx, vy, dc). there will be num_regions* points_in_each_region*3 float numbers,note that features are stored in flat format.
 */

void divide(float *pdata, int dim, int l, int *p_num_region, float **p_regions);

#endif

