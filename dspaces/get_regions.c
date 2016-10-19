/* minmax_reader.c : Example 2: Min/Max/Average of Array using DataSpace
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
// Example using array size, 128. 
// If modifying, MUST change in minmax_writer.c as well.
#define ARRAY_SIZE 128

int main(int argc, char **argv)
{
    // region definition
    // those parameters are obtained after dividing 
    int region_length;
    int num_region;



	int err;
	int nprocs, rank, i;
	MPI_Comm gcomm;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Barrier(MPI_COMM_WORLD);
	gcomm = MPI_COMM_WORLD;

	// DataSpaces: Initalize and identify application
	// Usage: dspaces_init(num_peers, appid, Ptr to MPI comm, parameters)
	// Note: appid for get.c is 2 [for put.c, it was 1]
	dspaces_init(nprocs, 2, &gcomm, NULL);
	
 	// Name our data.
	char var_name[128];
	sprintf(var_name, "ex3_sample_data");

	// DataSpaces: Read-Lock Mechanism
	// Usage: Prevent other processies from changing the 
	// 	  data while we are working with it
	dspaces_lock_on_read("my_test_lock", &gcomm);

    // !!!feng: we should know how many regions in total so that we can asign the pairs
    // MPI_receive here?
    int num_region;
    int num_tasks = num_region*(num_region-1)/2;

	// Each process will need to compute its DataSpace index
	int tasks_per_proc = num_tasks/nprocs;
	int tasks_left_over = num_tasks%nprocs;
    int pair_index_l, pair_index_h;
	int ds_lb_index, ds_ub_index;
	
	if(rank<tasks_left_over){
		pair_index_l = rank*(tasks_per_proc+1);
		pair_index_h = pair_index_l + tasks_per_proc;
	}else{
		pair_index_l = rank*tasks_per_proc+tasks_left_over;
		pair_index_h = ds_lb_index + tasks_per_proc-1;
	}

    region_memory_size = (region_length+1)*(region_length+1)*3*sizeof(float);

    for(i = pair_index_l; i<= pair_index_h; i++){
        // the index of the two pairs
        int a, b;
        float div;

        // get the region index of that pair
        get_pair_index(i, a, b);

        float *buffer a = (float *)malloc(region_memory_size);
        float *buffer b = (float *)malloc(region_memory_size);
        // Define the dimensionality of the data to be received 
        int ndim = 1; 
            
        // Prepare LOWER and UPPER bound dimensions
        uint64_t lb[3] = {0}, ub[3] = {0};
        // get buffer for left region
        lb[0]= a;
        ub[0] =a;
        dspaces_get(var_name, 1, region_memory_size, ndim, lb, ub, buffer_a);

        // get buffer for right region
        lb[0]= a;
        ub[0] =a;
        dspaces_get(var_name, 1, region_memory_size, ndim, lb, ub, buffer_b);



        // we can get the divergence now!

        div = get_divs( regions + i*d2*d3 , regions + j*d2*d3, region_length, k, div_func);
        // put the divergence into divergence matrix
        char div_name[128];
        sprintf(div_name, "div_data");
		dspaces_lock_on_write("div_lock", &gcomm);

        ndim = 2
        lb[3] = {0}, ub[3] = {0};
        lb[0] = a, ub[0] = a;
        lb[1] = b, ub[1] = b;

		dspaces_put(var_name, timestep, sizeof(float), ndim, lb, ub, div);

        // how about the symmetric part?

		dspaces_unlock_on_write("div_lock", &gcomm);
    }
	

    // no need to reduce

	// Report data to user
	if(rank==0){
        // get the clustering done here
		printf("Max: %d, Min: %d, Average: %d\n",global_max,global_min,global_avg);		
	}

	// DataSpaces: Finalize and clean up DS process
	dspaces_finalize();

	MPI_Barrier(gcomm);
	MPI_Finalize();

	return 0;
}
