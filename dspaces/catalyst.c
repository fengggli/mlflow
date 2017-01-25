#include <stdint.h>
#include <unistd.h>
#include "dataspaces.h"
#include "mpi.h"
#include "region_def.h"
#include "common_utility.h"


int main(int argc, char **argv)
{
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
	dspaces_init(1, 4, &gcomm, NULL);

    /*
     * npdiv parameter
     */

    char msg[STRING_LENGTH];
    sprintf(msg, "ds init ok" );
    my_message(msg, rank, LOG_CRITICAL);

        /*
    sprintf(msg, "dataspaces init complete");
    my_message(msg, rank);
    */

	// Timestep notation left in to demonstrate how this can be adjusted
	int timestep=0;

    char lock_name_vel[STRING_LENGTH];

    int dims[3] = {1, POINTS_SIDE, POINTS_SIDE};
    int num_points = dims[0]*dims[1]*dims[2];


    // we will receive each timestamp
    while(timestep < MAX_VERSION){
        snprintf(lock_name_vel, STRING_LENGTH, "vel_lock");


            if(rank == 0){
                sprintf(msg, "\n********************timestep %d now start!\n",timestep);
                my_message(msg, rank, LOG_WARNING);
            }



            uint64_t lb[3] = {0}, ub[3] = {0};
            // save divergence in a 1-d array! this will save space
            int ndim = 3;
            char var_name_vel[STRING_LENGTH];
            sprintf(var_name_vel, "VEL");
            /*

            sprintf(msg, "now div variable has dimention %d",num_tasks);
            my_message(msg, rank, LOG_WARNING);
            */

            double t1, t2, t3;
            int i, j, ret_get;


            // now anlysis
            // get the clustering done here

            //reconstruct the divergence matrix, this is a ragged matrix, only le 
            //see the distancematrix function in cluster.c 
            int count = 0;
            int error_flag = 0;

            float* vel_data = (float*)malloc(num_points*sizeof(float)*3);

            lb[0] = 0; 
            ub[0] = num_points-1;

            // every rank need at least acquire the lock
            sprintf(msg, "try to acquired vel read lock %s", lock_name_vel);
            my_message(msg, rank, LOG_WARNING);

            dspaces_lock_on_read(lock_name_vel, &gcomm);
            sprintf(msg, "acquired vel read lock");
            my_message(msg, rank, LOG_WARNING);



            t1 = MPI_Wtime();
            ret_get = dspaces_get(var_name_vel, timestep, sizeof(float)*3 , ndim, lb, ub, vel_data);
            t2 = MPI_Wtime();

            dspaces_unlock_on_read(lock_name_vel, &gcomm);

            sprintf(msg, "vel read lock released ");
            my_message(msg, rank, LOG_WARNING);
            
            // reconstruct the matrix
            if(ret_get != 0){
                        error_flag  = 1;
                        exit(-1);
                    }
            sprintf(msg, "all the vel is read from dataspaces");
            my_message(msg, rank, LOG_WARNING);

            my_message(msg, rank, LOG_WARNING);

            if(error_flag == 1){
                sprintf(msg, "ERROR when read vel from Dspaces");
                my_message(msg, rank, LOG_CRITICAL);
            }


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
