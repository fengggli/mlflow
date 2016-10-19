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


#endif
