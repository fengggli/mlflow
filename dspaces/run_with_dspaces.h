#ifdef  RUN_WITH_DSPACES_H
#define RUN_WITH_DSPACES_H

#include "divide.h"

// read a block cutout from hdf file
// divide into regions
// for each region get the density list
// for each pair of regions, get the divergence using their density lists
// output all the divergence into a txt file


extern int read_data(const char* file_name, float **pressure, float **velocity, int *dim1, int *dim2, int *dim3);
extern int free_data(float *pressure, float *velocity);
extern float get_divs(float *A, float *B, int region_length, int k, int div_func);

/* read a block cutout from hdf file
* divide into regions
* input: 
*   hdfpath:path of hdf source
*   region_length: length of 10 region will have 11*11 121 points
* output:
*   p_num_of_regions: number of regions
*   p_regions: buffer for all the regions
* return:
*   1 if succeed
*/


int generate_regions(char *hdfpath, int region_length, int *p_num_region, float **p_regions);


#endif
