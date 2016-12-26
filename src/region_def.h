#ifndef REGION_DEF_H
#define REGION_DEF_H
// define data block size, region size
// version number is from 1,2,3,4... MAX_VERSION

// use synthetic data to verify
#define USE_SYNTHETIC 

#ifdef USE_SYNTHETIC
    #define MAX_VERSION (1)
    #define MAX_VELOCITY (0.5)
    #define REGION_LENGTH (2)

    #define POINTS_SIDE (41)
    #define NUM_REGION (400)

    #define K_NPDIV (4)

    // clustering
    //#define NPASS (1)
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


// points in each side  of this data block from hdf5 file
// 51 or 201
/*
#define POINTS_SIDE (51)
#define NUM_REGION (25)
*/


#define TIMING

// max length of a string
#define STRING_LENGTH (80)

// frequency when print divergence result
#define PER_FREQ (5)

// sequential run
#define SEQ


/*
 * NPDIV related
 */

/*
 * clustering related
 */
#
#define NCLUSTERS (3)


#endif
