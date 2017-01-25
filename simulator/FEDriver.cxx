#include "FEDataStructures.h"
#include <mpi.h>

#ifdef USE_CATALYST
#include "FEAdaptor.h"
#endif

// added by feng 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "dataspaces.h"
#include "common_utility.h"
#include "region_def.h"
#include "string.h"

using namespace std;

// this will get all vel and pres data
void get_raw_buffer(int timestep, void *extra_info, int rank, MPI_Comm * p_gcomm,char * var_name_vel, float **p_buffer_vel, char * var_name_pres, float **p_buffer_pres,  double *p_time_used){
    char msg[STRING_LENGTH];
    double t1, t2;
    int ret_get = -1;

    float *vel_data = *p_buffer_vel;

    if(extra_info != NULL){
        printf("no extra info required\n");
        exit(-1);
    }

    // data layout
    int dims[3] = {1, POINTS_SIDE, POINTS_SIDE};
    int num_points = dims[0]*dims[1]*dims[2];

    size_t elem_size_vel = sizeof(float)*3;
    size_t elem_size_pres = sizeof(float);
    
    // prepare to read regions from dataspaces
    uint64_t lb[3] = {0}, ub[3] = {0};
    lb[0] = 0;
    ub[0] = num_points - 1;

    // Define the dimensionality of the data to be received 
    int ndim = 3;

    

    char lock_name_vel[STRING_LENGTH];
    //snprintf(lock_name_vel, STRING_LENGTH, "vel_lock_t_%d", timestep);
    snprintf(lock_name_vel, STRING_LENGTH, "vel_lock");

    char lock_name_pres[STRING_LENGTH];
    //snprintf(lock_name_pres, STRING_LENGTH, "pres_lock_t_%d", timestep);
    snprintf(lock_name_pres, STRING_LENGTH, "pres_lock");

    

    sprintf(msg, "try to acquired the vel read lock %s", lock_name_vel );
    my_message(msg, rank, LOG_WARNING);
    dspaces_lock_on_read(lock_name_vel, p_gcomm);

    sprintf(msg, "get the  the vel read lock");
    my_message(msg, rank, LOG_WARNING);

    // read all regions in once
    t1 = MPI_Wtime();

    printf("var name is %s, timstep: %d, elem_size_vel = %d, ndim =%d. lb=[%d, %d, %d], hb=[%d, %d, %d], one data %3.4f", var_name_vel, timestep, elem_size_vel, ndim, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2], vel_data[num_points-1]);
    ret_get = dspaces_get(var_name_vel, timestep, elem_size_vel, ndim, lb, ub, vel_data);

    t2 = MPI_Wtime();

    // now we can release region lock
    dspaces_unlock_on_read(lock_name_vel, p_gcomm);
    sprintf(msg, "release the vel read lock");
    my_message(msg, rank, LOG_WARNING);

    if(ret_get != 0){
        perror("get all vel error, now exit");
        printf("error number %d \n", ret_get);
        exit(-1);
    }else{
        sprintf(msg, "read %d vel from dspaces, each has %zu bytes", num_points, elem_size_vel);
        my_message(msg, rank, LOG_WARNING);
    }


    *p_time_used = t2-t1;
    
    // do the same for pres data
    /*
    sprintf(msg, "try to acquired the pres read lock %s", lock_name_pres );
    my_message(msg, rank, LOG_WARNING);
    dspaces_lock_on_read(lock_name_pres, p_gcomm);

    sprintf(msg, "get the  the pres read lock");
    my_message(msg, rank, LOG_WARNING);

    // read all regions in once
    t1 = MPI_Wtime();

    ret_get = dspaces_get(var_name_pres, timestep, elem_size_pres, ndim, lb, ub, pres_data);

    t2 = MPI_Wtime();

    // now we can release region lock
    dspaces_unlock_on_read(lock_name_pres, p_gcomm);
    sprintf(msg, "release the pres read lock");
    my_message(msg, rank, LOG_WARNING);

    if(ret_get != 0){
        perror("get all pres error, now exit");
        printf("error number %d \n", ret_get);
        exit(-1);
    }else{
        sprintf(msg, "read %d pres from dspaces, each has %zu bytes", num_points, elem_size_pres);
        my_message(msg, rank, LOG_WARNING);
    }

    *p_buffer_pres = pres_data;
    *p_time_used += t2-t1;
    */
}
// Example of a C++ adaptor for a simulation code
// where the simulation code has a fixed topology
// grid. We treat the grid as an unstructured
// grid even though in the example provided it
// would be best described as a vtkImageData.
// Also, the points are stored in an inconsistent
// manner with respect to the velocity vector.
// This is purposefully done to demonstrate
// the different approaches for getting data
// into Catalyst. Note that through configuration
// that the driver can be run without linking
// to Catalyst.

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  /*
  Grid grid;
  unsigned int numPoints[3] = {201, 201, 1};
  double spacing[3] = {1, 1, 0 };
  grid.Initialize(numPoints, spacing);
  Attributes attributes;
  attributes.Initialize(&grid);
  */


   /*
   * start of my definition
   */
    // initialize dataspaces
    int nprocs, rank, ret;
    MPI_Comm gcomm;

    char result_path[STRING_LENGTH]="";

    // output results into a folded specified with a slurm jobid
    if(argc == 2){
        strcpy(result_path, argv[1]);
    }



    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
    gcomm = MPI_COMM_WORLD;

    //char msg[STRING_LENGTH];

// DataSpaces: Initalize and identify application
// Usage: dspaces_init(num_peers, appid, Ptr to MPI comm, parameters)
// Note: appid for get.c is 2 [for put.c, it was 1]
    ret = dspaces_init(1, 2, &gcomm, NULL);

    if(ret == 0){
        printf("dataspaces init successfully");
    }else{
        printf("dataspaces init error");
        exit(-1);

    }


  // vel and pressure buffer
  double time_used;

    char var_name_vel[STRING_LENGTH];
    char var_name_pres[STRING_LENGTH];
    sprintf(var_name_vel, "VEL");
    sprintf(var_name_pres, "PRES");


    int dims[3] = {1, POINTS_SIDE, POINTS_SIDE};
    uint64_t num_points = dims[0]*dims[1]*dims[2];
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

    // prepare space for pres
    float * pres_data = (float *)malloc(num_points* sizeof(float));
    if(pres_data == NULL){
          perror("pres data allocated error");
          exit(-1);
      }





  /*
   * end of my definition
   */



#ifdef USE_CATALYST
  //FEAdaptor::Initialize(argc, argv);
#endif
  //unsigned int numberOfTimeSteps = 100;
  int timestep = 0;
  for(timestep=0;timestep<MAX_VERSION;timestep++)
    {
    // use a time step length of 0.1
    double time = timestep * 0.1;

    // read data from dataspces
    // this will get blocked until new data available
    get_raw_buffer(timestep, NULL ,rank, &gcomm, var_name_vel, &vel_data, var_name_pres,  &pres_data, &time_used);

    //update using vel and pres info, if there is more ranks I need to partition first
   // attributes.UpdateFields(buffer_vel, buffer_pres);

    

    // also let the analysis part done.
    // add one pipeline 
    // show together
#ifdef USE_CATALYST
    //FEAdaptor::CoProcess(grid, attributes, time, timestep, timestep == numberOfTimeSteps-1);
#endif
    }
    // free buffer
    if(vel_data != NULL){
        free(vel_data);
        cout << "vel is freed" << endl;
    }
    if(pres_data != NULL){
        free(pres_data);
        cout << "pres is freed" << endl;
    }

#ifdef USE_CATALYST
  //FEAdaptor::Finalize();
#endif

    dspaces_finalize();

    MPI_Barrier(gcomm);
    MPI_Finalize();

  return 0;
}

