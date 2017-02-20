// this is calller to catalyst routine
// it wil read from dspaces and send to catalyst adaptor

#include "FEDataStructures.h"
#include <mpi.h>

#ifdef USE_CATALYST
#include "FEAdaptor.h"
#endif

// added by feng 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "ds_adaptor.h"

using namespace std;

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

  Grid grid;

  unsigned int dims[3] = {POINTS_SIDE, POINTS_SIDE, 1};
  uint64_t num_points = dims[0]*dims[1]*dims[2];

  //double spacing[3] = {1, 1, 0 };
  double spacing[3] = {0.1, 0.1, 0 };
  grid.Initialize(dims, spacing);
  Attributes attributes;
  attributes.Initialize(&grid);

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
        printf("dataspaces init successfully\n");
    }else{
        printf("dataspaces init error\n");
        exit(-1);

    }


    // vel and pressure buffer
    double time_used;

    char var_name_vel[STRING_LENGTH];
    char var_name_pres[STRING_LENGTH];
    sprintf(var_name_vel, "VEL");
    sprintf(var_name_pres, "PRES");


    /*

    uint64_t num_points = dims[0]*dims[1]*dims[2];

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
  FEAdaptor::Initialize(argc, argv);
#endif
  //unsigned int numberOfTimeSteps = 100;
  int timestep = 0;
  for(timestep=0;timestep<MAX_VERSION;timestep++)
    {
    // use a time step length of 0.1
    double time = timestep * 0.1;

    if(rank == 0){
        cout <<"-----current timestep" << timestep << endl;
    }

    // read data from dataspces
    // this will get blocked until new data available
    get_raw_buffer(timestep, NULL ,rank, &gcomm, var_name_vel, &vel_data, var_name_pres,  &pres_data, &time_used);

    //update using vel and pres info, if there is more ranks I need to partition first
    attributes.UpdateFields(vel_data, pres_data);

    // also let the analysis part done.
    // add one pipeline 
    // show together
#ifdef USE_CATALYST
    FEAdaptor::CoProcess(grid, attributes, time, timestep, timestep == MAX_VERSION-1);
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
  FEAdaptor::Finalize();
#endif

    dspaces_finalize();

    MPI_Barrier(gcomm);
    MPI_Finalize();

  return 0;
}

