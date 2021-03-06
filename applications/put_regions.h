#ifndef PUT_REGIONS_H
#define PUT_REGIONS_H
/*
 * it should  work like this:
 * note , currently only one process
 * 1. read data from hdf5
 * 2. divide into regions
 * 3. write regions in to dataspaces
 */

/* orignal comment from minmax_writer.c : Example 3: DataSpaces put 128 array
 * Bounding Box: 0,0,0 - 127,0,0 
 * 1 element at each space in the box
 * */

#include <stdint.h>
#include <unistd.h>
#include "run_with_dspaces.h"
#include "ds_adaptor.h"
// Size of array, if changing
// MUST also change in minmax_reader.c
#define ARRAY_SIZE 128

// read hdffile and generate all the regions
extern int generate_regions(char *hdfpath, int region_length, int *p_num_region, float **p_regions);

#endif


