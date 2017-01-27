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





// Example of a simulation code
int main(int argc, char* argv[])
{

    // dataspaces preparation
    int nprocs, rank, ret;
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
    ret = dspaces_init(1, 1, &gcomm, NULL);

    if(ret == 0){
        sprintf(msg, "dataspaces init successfully");
        my_message(msg, rank, LOG_CRITICAL);
    }else{
        sprintf(msg, "dataspaces init error");
        my_message(msg, rank, LOG_CRITICAL);
        exit(-1);

    }

    // data layout
    int dims[3] = {1, POINTS_SIDE, POINTS_SIDE};
    int num_points = dims[0]*dims[1]*dims[2];

    // time 
    double time_comm_vel;

    
    char var_name_vel[STRING_LENGTH];
    char var_name_pres[STRING_LENGTH];
    sprintf(var_name_vel, "VEL");
    sprintf(var_name_pres, "PRES");
    /*
    uint64_t gdim_vel[3] = {num_points,1,1};
    dspaces_define_gdim(var_name_vel, 3, gdim_vel);

    uint64_t gdim_pres[3] = {num_points,1,1};
    dspaces_define_gdim(var_name_pres, 3, gdim_pres);
    */
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


        //put_vel_buffer(timestep, NULL,rank, &gcomm, &vel_data, &time_comm_vel);
        put_raw_buffer(timestep, NULL,rank, &gcomm, var_name_vel,  &vel_data,var_name_pres, &pres_data, &time_comm_vel);
    }

    // free
    if(vel_data != NULL){
        free(vel_data);
    }
    if(pres_data != NULL){
        free(pres_data);
    }

    dspaces_finalize();

    MPI_Barrier(gcomm);
    MPI_Finalize();
    return 0;
}

