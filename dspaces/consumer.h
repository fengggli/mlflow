#ifndef CONSUMER_H
#define CONSUMER_H
/*
 * work like this
 *  1. read from dspaces(raw velocity data)
 *  2. divide into regions
 *  3. get all the divergences
 *  4. write divs into dspces
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
#include "region_def.h"
#include "common_utility.h"
#include "get_divs.h"
// Size of array, if changing
// MUST also change in minmax_reader.c
#define ARRAY_SIZE 128

// read hdffile and generate all the regions
extern int generate_regions(char *hdfpath, int region_length, int *p_num_region, float **p_regions);

#endif


