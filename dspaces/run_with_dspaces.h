#ifndef RUN_WITH_DSPACES_H
#define RUN_WITH_DSPACES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#include "divide.h"
#include "read_file.h"

// get current time
double get_cur_time() {
  struct timeval   tv;
  struct timezone  tz;
  double cur_time;

  gettimeofday(&tv, &tz);
  cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;

  return cur_time;
}

// read a block cutout from hdf file
// divide into regions
// for each region get the density list
// for each pair of regions, get the divergence using their density lists
// output all the divergence into a txt file

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
