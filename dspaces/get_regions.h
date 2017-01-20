#ifndef GET_REGIONS
#define GET_REGIONS
/* comment from original minmax_reader.c : Example 2: Min/Max/Average of Array using DataSpace
 * In this example, we will use a number of processes (specified by -np)
 * to compute the minimum and maximum element in an array and to compute
 * the average of all the values in the array.
 * You will see how DataSpaces accesses the values without reading from disk. 
*/


//notes
//it should run like this:
//1. assigned with pairs of regions
//2. get regions from dspaces
//3. calculate divergence
//4. put the divergence
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "dataspaces.h"
#include "mpi.h"
#include "cluster.h"
#include "get_divs.h"
#include "string.h"
#include "region_def.h"
#include "common_utility.h"
// Example using array size, 128. 
// If modifying, MUST change in minmax_writer.c as well.
#define ARRAY_SIZE 128


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

/*
 * we can use pair index to lookup this table for region index 
 * each entry like this:
 *  pairindex(key) index_a index_b
 *  (1,0),0 
 *  (2,0),1 (2,1),2
 *  (3,0),3 (3,1),4 (3,2) 5 
 *
 */
void generate_lookup_table(int num_region, int **p_table);


/*
 * free the lookuptable
 */
void free_lookup_table(int * table);


/* get the region index based on pair index
 * input:
 *  i: index of pair
 *  range: how many elements
 * output:
 *  a, b: index of two regions, a>b, start from (1,0)!
 *  linear time
 */
void get_pair_index(int *table, int index_pair,int *a, int *b);


/* 
 * modified on Jan 19
 * calculate all the divergence for this rank
 */
void cal_local_divs(float *buffer_all_regions,Region_Def * p_region_def, int k_npdiv, int div_func, int *table, int  pair_index_l,int  pair_index_h,  int rank, float** p_div_this_rank, double *time_used);

/* 
 * modified on Jan 19
 * get region info from dspaces
 */
void get_region_buffer(int timestep, Region_Def *p_region_def, int  rank, MPI_Comm * p_gcomm,  float ** p_buffer_all_regions, double * time_used);


/*
 * modified on Jan 19
 * put divergence into into dspaces
 */
void put_divs_buffer(int timestep,int pair_index_l, int pair_index_h , int num_tasks, int rank, MPI_Comm *p_gcomm, float ** p_divs_this_rank, double* time_used);



#endif
