#ifndef DIVIDE_H 
#define DIVIDE_H 
/*
 * Sep 7, by Feng Li, IUPUI
 * Devide a slice of turbulence data into regions
 * before splitting, we have a block of a*b*1, default value of c is 1
 * each region is a block with a length of l, ovelapping length default value is 0
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * generate sample regions
 * input:
 *  region_length, 10
 *  num_types: how many types of regions
 * output:
 *  sample_regions: buffer to save sample regions
 */
void gen_region_samples(float **sample_regions, int num_types,  int region_length);

/*
 * free sample regions
 * input:
 *  sample_regions: buffer to save sample regions
 *  num_types: number of types of regions
 */
void free_region_samples(float **sample_regions, int num_types);
    

/*
 * fill up regions with specific patters
 * input:
 */

void fill_region(float* to, float *from, size_t region_memory_size);

/*
 * Dec 22
 * generate synthetic regions 
 * input
 *  pdata, pointer to data buffer
 *  dim: num of points in x and y dimension 
 *  l: region_length: for 11 points in each size, the length will be 10
 * output:
 *  p_num_region; number of regions
 *  p_regions: the tripple features(vx, vy, dc). there will be num_regions* points_in_each_region*3 float numbers,note that features are stored in flat format.
 */

void divide_synthetic(int dim, int l, int *p_num_region, float **p_regions);

/*
 * given a large portion of multi-dimensional data, divide it into small regions
 * input
 *  pdata, pointer to data buffer
 *  dim: num of points in x and y dimension 
 *  l: region_length: for 11 points in each size, the length will be 10
 * output:
 *  p_num_region; number of regions
 *  p_regions: the tripple features(vx, vy, dc). there will be num_regions* points_in_each_region*3 float numbers,note that features are stored in flat format.
 */

void divide(float *pdata, unsigned int dims[3], int l, int *p_num_region, float *buffer_regions);

#endif
#ifdef __cplusplus
}
#endif

