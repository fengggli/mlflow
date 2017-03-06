#ifndef REGION_DEF_H
#define REGION_DEF_H

// this file contains key configurations

#ifdef __cplusplus
extern "C" {
#endif

// define data block size, region size
// version number is from 1,2,3,4... MAX_VERSION

// use synthetic data to verify
// only tested in sequential implementation
//#define USE_SYNTHETIC 
//#define USE_CAVITY
#define INCLUDE_ML
#define USE_PARAL_CAVITY
#define SAMPLING
#define SAMPLE_SIZE (20)
 
// large scale experiment
#define USE_EXP

// whether or not using same lock eachtime, only affect the ds_adaptor
// only affect raw data, writer won't start new
//#define USE_SAME_LOCK

#include "stdlib.h"
#include "stdio.h"

#ifdef USE_EXP
    #define MAX_VERSION (100)
    //#define CASE_LENGTH (1024)
    //#define REGION_LENGTH (32)
    #define CASE_LENGTH (32)
    #define REGION_LENGTH (4)
    #ifndef PROCS_PER_DIM
        #error "need to define process per side"
    #endif
    // need to define from outside
    //#define PROCS_PER_DIM (2)
    // this will be 40*4+1 = 161 (161 points in each side for 4^2 = 16 procs)
    #define POINTS_SIDE (CASE_LENGTH*PROCS_PER_DIM) //2^13
    #define NUM_REGION (POINTS_SIDE/REGION_LENGTH)*(POINTS_SIDE*REGION_LENGTH) // (2^13/2^8)^2 = 2^10 regions

    #define K_NPDIV (5)
    #define NPASS (100)


#else
    #ifdef USE_CAVITY
        #define MAX_VERSION (100)
        // side length(points num -1) in each region
        #define REGION_LENGTH (4)
        #define POINTS_SIDE (41)
        #define NUM_REGION (100)

        // this is the k in density estimation
        #define K_NPDIV (5)

        // clustering
        #define NPASS (100)
    #else
        #ifdef USE_PARAL_CAVITY
            #define MAX_VERSION (100)
            #define CASE_LENGTH (40)
            #define REGION_LENGTH (10)
            #ifndef PROCS_PER_DIM
                #error "need to define process per side"
            #endif
            // need to define from outside
            //#define PROCS_PER_DIM (2)
            // this will be 40*4+1 = 161 (161 points in each side for 4^2 = 16 procs)
            #define POINTS_SIDE (CASE_LENGTH*PROCS_PER_DIM +1) //81
            #define NUM_REGION (((POINTS_SIDE -1.0)/REGION_LENGTH)*((POINTS_SIDE-1.0)/REGION_LENGTH))

            #define K_NPDIV (5)
            #define NPASS (100)
        #else
            #define MAX_VERSION (10)
            // side length(points num -1) in each region
            #define REGION_LENGTH (10)
            #define POINTS_SIDE (201)
            #define NUM_REGION (400)

            // this is the k in density estimation
            #define K_NPDIV (5)

            // clustering
            #define NPASS (100)
        #endif
    #endif
#endif


// points in each side  of this data block from hdf5 file
// 51 or 201
/*
#define POINTS_SIDE (51)
#define NUM_REGION (25)
*/


#define TIMING

// max length of a string
#define STRING_LENGTH (160)

// frequency when print divergence result
#define PER_FREQ (5)

// sequential run
//#define SEQ


/*
 * NPDIV related
 */

/*
 * clustering related
 */
#
#define NCLUSTERS (3)

typedef struct{
    int region_length;
    int side_num_region;
    int num_region;
    int region_num_cell;
    size_t region_memory_size;
}Region_Def;

//region definition
void fill_region_def(Region_Def *p_region_def);

// extract region definition
void extract_region_def(Region_Def *p_region_def, int *p_region_length,int * p_side_num_region, int *p_num_region,int * p_region_num_cell, size_t *p_region_memory_size);

#endif
#ifdef __cplusplus
}
#endif

