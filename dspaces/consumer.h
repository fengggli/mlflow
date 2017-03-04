#ifndef CONSUMER_H
#define CONSUMER_H
/*
 * work like this, multiple processes
 *  1. read from dspaces(raw velocity data)
 *  2. divide into regions(currently allocate new space for region, this will waste space)
 *  3. get all the divergences of pairs of each rank
 *  4. write divs into dspces
 *
 * what should be defined
 *  1. region size
 *  2. how many ranks
 */



/* orignal comment from minmax_writer.c : Example 3: DataSpaces put 128 array
 * Bounding Box: 0,0,0 - 127,0,0 
 * 1 element at each space in the box
 * */

#include <stdint.h>
#include <unistd.h>
#include "ds_adaptor.h"
#include "run_with_dspaces.h"
#include "get_divs.h"
// Size of array, if changing
// MUST also change in minmax_reader.c
#define ARRAY_SIZE 128

// read hdffile and generate all the regions
extern int generate_regions(char *hdfpath, int region_length, int *p_num_region, float **p_regions);

/*
 * main routine to calculate divs which is assigned to each rank
 * INPUT:
 *      buffer_all_regions, all velocify info but orderd by region id
 * OUTPUT
 *      p_divs_this_rank, all divs of this rank
 */
void cal_local_divs(float *buffer_regions, int region_length, int k_npdiv, int div_func, int *table, int  pair_index_l,int  pair_index_h,  int rank, float*divs_this_rank, double *p_time_used);


/* sample s subset of regions
 * input:
 *      buffer_region, whole dataset
 *      elem_size_region, memory size of one region
 *      num_elems_region, how many regions in whole dataset
 *      num_elemns_sample, how many samples
 * output:
 *      buffer_sample, 
 */

static void prepare_sampled_buffer(float *buffer_region, float *buffer_sample, int num_elems_region, int num_elems_sample);


/*
 * assign cluster id to stripped regions using global samples and clusterid
 * input:
 *      k: number of clusters
 *      k_npdiv: how many nearest neighbours
 *      div_func: which div function
 */

static void assign_clusterid(float *buffer_region,int num_region, float * buffer_sample_all, int num_sample_all,float *buffer_medoids,int k,int region_length,int k_npdiv, int div_func,  float *clusterids);



#endif


