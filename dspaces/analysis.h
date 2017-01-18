#ifndef PUT_REGIONS_H
#define PUT_REGIONS_H
/*
 * it should  work like this:
 * 1. only one process here, this process will get all the regions and put regions into dataspaces
 * 2. put the region data into daspaces
 * 3. initialize the 'divergence matrix' dataspace variable
 */

/* orignal comment from minmax_writer.c : Example 3: DataSpaces put 128 array
 * Bounding Box: 0,0,0 - 127,0,0 
 * 1 element at each space in the box
 * */

#include <stdint.h>
#include <unistd.h>
#include "dataspaces.h"
#include "run_with_dspaces.h"
#include "mpi.h"
#include "cluster.h"
#include "region_def.h"
#include "common_utility.h"
// Size of array, if changing
// MUST also change in minmax_reader.c
//#define ARRAY_SIZE 128

// read hdffile and generate all the regions
//extern int generate_regions(char *hdfpath, int region_length, int *p_num_region, float **p_regions);

#endif


