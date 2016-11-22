#ifndef READ_FILE_H
#define READ_FILE_H

#include "hdf5.h"
#include <stdlib.h>
//#include "read_file.h"
//#define FILE        "data/test_1_2_3_4.h5"
//#define FILE        "data/isotropic_255_255_5.h5"
#define DATASET_NAME
//#define DATASET_P "p00000" // dataset name for pressure
//#define DATASET_U "u00000" // dataset name for velocity 

#define RANK         4 // four dims: x y z and time 
#define p_malloc_error(I) {printf("malloc failuren at %d",I);exit(-1);} 


/*
 * * Get datatype and dataspace identifiers and then query
 * * read only one time slice
 * * dataset class, order, size, rank and dimensions.
 * */

int read_data(const char* file_name, float **pressure, float **velocity, int *dim1, int *dim2, int *dim3);

int free_data(float *pressure, float *velocity);
#endif

