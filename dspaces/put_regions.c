#include "put_regions.h"

int main(int argc, char **argv)
{
	int err;
	int nprocs, rank;
	MPI_Comm gcomm;

    char * hdfpath = "data/isotropic_201_201_1.h5";
    int region_length = 10;
    int num_region;
    float *regions;
    size_t region_memory_size; 

    // MPI communicator
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Barrier(MPI_COMM_WORLD);
	gcomm = MPI_COMM_WORLD;

    // if rank == 0 get the data and divide into regions
    // other processes will wait here
    generate_regions(hdfpath, region_length, &num_region, &regions);
    printf("are regions are generated\n");

    // how large is one region?
    region_memory_size = (region_length+1)*(region_length+1)*3*sizeof(float);

	// Initalize DataSpaces
	// # of Peers, Application ID, ptr MPI comm, additional parameters
	// # Peers: Number of connecting clients to the DS server
	// Application ID: Unique idenitifier (integer) for application
	// Pointer to the MPI Communicator, allows DS Layer to use MPI barrier func
	// Addt'l parameters: Placeholder for future arguments, currently NULL.
	dspaces_init(1, 1, &gcomm, NULL);

	// Timestep notation left in to demonstrate how this can be adjusted
	int timestep=0;

    // we only need one timestamp
	while(timestep<1){
		timestep++;

		// DataSpaces: Lock Mechanism
		// Usage: Prevent other process from modifying 
		// 	  data at the same time as ours
		dspaces_lock_on_write("velocity_lock", &gcomm);

		//Name the Data that will be writen
		char var_name[128];
		sprintf(var_name, "velocity_data");

        // each "cell" is a region 
		// ndim: Dimensions for application data domain
		// In this case, our data array will be 1 dimensional
		int ndim = 1; 
		
		// Prepare LOWER and UPPER bound dimensions, init 0s
        // feng: assign correct bound based on rank and total number of regions
		uint64_t lb[3] = {0}, ub[3] = {0};

		// Bounds: 0,0,0 to 127,0,0
		ub[0] = num_region-1;

		// DataSpaces: Put data array into the space
		// 1 integer in each box, fill boxes 0,0,0 to 127,0,0
		dspaces_put(var_name, timestep, region_memory_size, ndim, lb, ub, regions);

		free(data);
		
		// DataSpaces: Release our lock on the data
		dspaces_unlock_on_write("velocity_lock", &gcomm);

	}

	// DataSpaces: Finalize and clean up DS process
	dspaces_finalize();

	MPI_Barrier(gcomm);
	MPI_Finalize();

	return 0;
}
