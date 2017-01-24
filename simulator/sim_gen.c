//#include "FEDataStructures.h"
#include <mpi.h>
#include <stdlib.h>

//#include "FEAdaptor.h"
#include "region_def.h"
#include "common_utility.h"
#include "dataspaces.h"

void update_attributes(int timestep, int *dims,float *vel_data, float*p_data){
    int i, j, k;

    // writers
    float *tmp_vel = vel_data;
    float *tmp_p= p_data;

    for(i = 0; i< dims[0]; i++){
        for(j = 0; j < dims[1]; j++){
            for(k = 0; k < dims[2]; k++){
                // vel data
                tmp_vel[0] = j*timestep;
                tmp_vel[1] = 0;
                tmp_vel[2] = 0;

                // pressure
                tmp_p[0] = 0;

                tmp_vel +=3;
                tmp_p +=1;
            }
        }
    }
}

void put_vel_buffer(int timestep, Region_Def *p_region_def, int rank, MPI_Comm * p_gcomm, float **p_buffer_vel, double *p_time_used){
    char msg[STRING_LENGTH];
    double t1, t2;
    int ret_put = -1;
    
    if(p_region_def != NULL){
        printf("error, no region info required here\n");
        exit(-1);
    }


    // data layout
    int dims[3] = {1, POINTS_SIDE, POINTS_SIDE};
    int num_points = dims[0]*dims[1]*dims[2];

    size_t elem_size_vel = sizeof(float)*3;
    
    // prepare to read regions from dataspaces
    uint64_t lb[3] = {0}, ub[3] = {0};
    lb[0] = 0;
    ub[0] = num_points - 1;

    // Define the dimensionality of the data to be received 
    int ndim = 3;

    char var_name_vel[STRING_LENGTH];
    sprintf(var_name_vel, "VEL");

    uint64_t gdim_vel[3] = {num_points,1,1};
    dspaces_define_gdim(var_name_vel, 3, gdim_vel);

    char lock_name_vel[STRING_LENGTH];
    snprintf(lock_name_vel, STRING_LENGTH, "vel_lock_t_%d", timestep);

    
    // prepare space
    float * vel_data = (float *)malloc(num_points* sizeof(float)*3);
    if(vel_data == NULL){
          perror("vel data allocated error");
          exit(-1);
      }


    sprintf(msg, "try to acquired the vel write lock %s", lock_name_vel );
    my_message(msg, rank, LOG_WARNING);
    dspaces_lock_on_write(lock_name_vel, p_gcomm);

    sprintf(msg, "get the  the vel write lock");
    my_message(msg, rank, LOG_WARNING);

    // read all regions in once
    t1 = MPI_Wtime();

    ret_put = dspaces_put(var_name_vel, timestep, elem_size_vel, ndim, lb, ub, vel_data);

    t2 = MPI_Wtime();

    // now we can release region lock
    dspaces_unlock_on_write(lock_name_vel, p_gcomm);
    sprintf(msg, "release the vel write lock");
    my_message(msg, rank, LOG_WARNING);

    if(ret_put != 0){
        perror("write all vel error, now exit");
        exit(-1);
    }else{
        sprintf(msg, "write %d vel to dspaces, each has %ld bytes", num_points, elem_size_vel);
        my_message(msg, rank, LOG_WARNING);
    }

    *p_buffer_vel = vel_data;
    *p_time_used = t2-t1;
}


// Example of a simulation code
int main(int argc, char* argv[])
{

    // dataspaces preparation
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
    dspaces_init(1, 1, &gcomm, NULL);

    sprintf(msg, "dataspaces init successfully");
    my_message(msg, rank, LOG_CRITICAL);

    // data layout
    int dims[3] = {1, POINTS_SIDE, POINTS_SIDE};
    int num_points = dims[0]*dims[1]*dims[2];

    // time 
    double time_comm_vel;

    // prepare space
    float * vel_data = (float *)malloc(num_points* sizeof(float)*3);
    if(vel_data == NULL){
        perror("vel data allocated error");
        exit(-1);
    }
    float * pres_data = (float *)malloc(num_points* sizeof(float));

    if(pres_data == NULL){
        perror("pressure data allocated error");
        exit(-1);
    }


    for(unsigned int timestep=0;timestep<MAX_VERSION;timestep++)
    {

        sprintf(msg, "\n*******************simulation timestep %d now start!\n",timestep);
        my_message(msg, rank, LOG_WARNING);


        // this will simulation process risides
        update_attributes(timestep, dims, vel_data, pres_data);

        put_vel_buffer(timestep, NULL,rank, &gcomm, &vel_data, &time_comm_vel);
    }

    // free
    if(vel_data != NULL){
        free(vel_data);
    }
    if(pres_data != NULL){
        free(pres_data);
    }
    return 0;
}

