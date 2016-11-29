#ifndef REGION_DEF_H
#define REGION_DEF_H
// define data block size, region size
// version number is from 1,2,3,4... MAX_VERSION
#define MAX_VERSION (10)

// side length(points num -1) in each region
#define REGION_LENGTH (10)

// points in each side  of this data block from hdf5 file
// 51 or 201
/*
#define POINTS_SIDE (51)
#define NUM_REGION (25)
*/

#define POINTS_SIDE (201)
#define NUM_REGION (400)
#define TIMING

// max length of a string
#define STRING_LENGTH (80)

// frequency when print divergence result
#define PER_FREQ (5)

#endif
