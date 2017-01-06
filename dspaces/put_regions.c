#include "put_regions.h"
//#define debug_1


int main(int argc, char **argv)
{
	int err;
	int nprocs, rank;
	MPI_Comm gcomm;

    
    // MPI communicator
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Barrier(MPI_COMM_WORLD);
	gcomm = MPI_COMM_WORLD;

    // Initalize DataSpaces
	// # of Peers, Application ID, ptr MPI comm, additional parameters
	// # Peers: Number of connecting clients to the DS server
	// Application ID: Unique idenitifier (integer) for application
	// Pointer to the MPI Communicator, allows DS Layer to use MPI barrier func
	// Addt'l parameters: Placeholder for future arguments, currently NULL.
    char msg[STRING_LENGTH];
    sprintf(msg, "try to init dataspaces");
    my_message(msg, rank, LOG_CRITICAL);

	dspaces_init(1, 1, &gcomm, NULL);


    sprintf(msg, "dataspaces init successfully");
    my_message(msg, rank, LOG_CRITICAL);


	// Timestep notation left in to demonstrate how this can be adjusted
	int timestep=0;

    //Name the Data that will be writen
    char var_name[STRING_LENGTH];
    sprintf(var_name, "region_data");

    char lock_name_regions[STRING_LENGTH];
    char lock_name_divs[STRING_LENGTH];

    /*
        snprintf(lock_name_regions, STRING_LENGTH, "region_lock_same");
        snprintf(lock_name_divs, STRING_LENGTH, "div_lock_same", timestep);
        */
    // we will receive each timestamp
    while(timestep < MAX_VERSION){
        
        snprintf(lock_name_regions, STRING_LENGTH, "region_lock_t_%d", timestep);
        snprintf(lock_name_divs, STRING_LENGTH, "div_lock_t_%d", timestep);

            if(rank == 0){
                sprintf(msg, "\n********************timestep %d now start!\n",timestep);
                my_message(msg, rank, LOG_WARNING);
            }

            //char * hdfpath = "data/isotropic_201_201_1.h5";
            char hdfpath[STRING_LENGTH];

            // 
            sprintf(hdfpath, "data/isotropic_%d_%d_1_t_%d.h5",POINTS_SIDE,POINTS_SIDE, timestep);

            //int region_length = 10;
            int num_region = -1;
            float *regions;
            size_t region_memory_size; 


            // dataspaces access return value
            int ret_put = -1;
            
            // how large is one region?
            region_memory_size = (REGION_LENGTH+1)*(REGION_LENGTH+1)*3*sizeof(float);

            // validate the path
            if( access( hdfpath, F_OK ) == -1){
#ifndef USE_SYNTHETIC
                sprintf(msg, "path %s does not exist", hdfpath);
                my_message(msg, rank, LOG_CRITICAL);
                exit(-1);
#endif
            }else{
                generate_regions(hdfpath, REGION_LENGTH, &num_region, &regions);
                sprintf(msg, "%d regions are generated, each region has size %ld bytes", num_region, region_memory_size);
                my_message(msg, rank, LOG_WARNING );
            }
            

            // each "cell" is a region 
            // ndim: Dimensions for application data domain
            // In this case, our data array will be 1 dimensional
            int ndim = 3; 
            
            // Prepare LOWER and UPPER bound dimensions, init 0s
            // feng: assign correct bound based on rank and total number of regions
            uint64_t lb[3] = {0}, ub[3] = {0};

            // Bounds: 0,0,0 to 127,0,0
            ub[0] = num_region-1;

            // DataSpaces: Lock Mechanism
            // Usage: Prevent other process from modifying 
            // 	  data at the same time as ours
            sprintf(msg, "try to acquired the region write lock %s", lock_name_regions);
            my_message(msg, rank, LOG_WARNING);

            dspaces_lock_on_write(lock_name_regions, &gcomm);
            sprintf(msg, "acquired the region write lock");
            my_message(msg, rank, LOG_WARNING);


            // DataSpaces: Put data array into the space
            // 1 integer in each box, fill boxes 0,0,0 to 127,0,0

            // check the variables 
#ifdef debug_1
            int ii;
            printf("contents of the points in regions 1 before dividing:\n");
            int num_cell = (REGION_LENGTH+1)*(REGION_LENGTH +1);
            int region_id = 1;
            float* A;
            A = regions + region_id*(num_cell*3);
            for(ii == 0 ; ii <  num_cell; ii++){
                printf("\t point %d: (%.3f %.3f %.3f)\n", ii, *(A + 3*ii+ 0), *(A +3*ii + 1),*(A + 3*ii +2)); 
            }
#endif


            ret_put = dspaces_put(var_name, timestep, region_memory_size, ndim, lb, ub, regions);

            // DataSpaces: Release our lock on the data
            dspaces_unlock_on_write(lock_name_regions, &gcomm);
            sprintf(msg, "released the region write lock");
            my_message(msg, rank, LOG_WARNING);


            if(ret_put == 0){
                sprintf(msg, "%d regions are written into Dspaces",num_region);
            }
            else{
                sprintf(msg, "ERROR when writing regions into dspacs");
            }
            my_message(msg, rank, LOG_WARNING);

            free(regions);
                    

        timestep++;
	}

    sprintf(msg, "now finalize the dspaces and exit");
    my_message(msg, rank, LOG_CRITICAL);

	// DataSpaces: Finalize and clean up DS process
	dspaces_finalize();

	MPI_Barrier(gcomm);
	MPI_Finalize();

	return 0;
}
