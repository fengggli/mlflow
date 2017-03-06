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

static void map_regions(float *buffer_cluster, float *cluster_data, const unsigned int dims[3], const unsigned int region_length){
    int p, q, ii, jj, index_x, index_y, linear_index;

    int l = region_length;
    printf("dims== %d %d %d , region_length = %d\n", dims[0],dims[1], dims[2], region_length);

    float clusterid;

    if((dims[0])%l != 0|| dims[1]%l != 0){
        printf("not pefect division, try different region size of datacut size\n");
        exit(-1);
    }

    int num_region_row = dims[0]/l;
    int num_region_col =  dims[1]/l;

    for(p = 0; p < num_region_row; p++){
        for( q = 0; q < num_region_col; q++){
            clusterid = *(buffer_cluster + num_region_col*p + q);
            // for each region
            // this region will have corner:
            //  (pl, ql) , (pl, (q+1)l)
            //  ((p+1)l,ql), ((p+1)l,(q+1)l)
            //  also the region all have the center (pl+l/2, ql+l/2)
            for(ii = 0; ii < l; ii++){
                for(jj = 0; jj < l; jj++){
                    // for each point inside the region
                    index_x = p*l+ii;
                    index_y = q*l+jj;

                    // mapped the logic address into linear address
                    // will be overlap but it's okay
                    linear_index = index_x*dims[1] + index_y;
                    *(cluster_data +linear_index) = clusterid;
                }
            }
            //printf("    region %d-%d is tupled\n", p, q);
        }
    }

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

  Grid grid;

  unsigned int dims[3] = {POINTS_SIDE, POINTS_SIDE, 1};
  unsigned int num_points = dims[0]*dims[1]*dims[2];
  unsigned int region_length = REGION_LENGTH;
    


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
    printf("trying init dspaces for %d process\n", nprocs);
    ret = dspaces_init(nprocs, 4, &gcomm, NULL);

    if(ret == 0){
        printf("dataspaces init successfully\n");
    }else{
        printf("dataspaces init error\n");
        exit(-1);

    }


    // vel and pressure buffer
    //double time_used, time_used_cluster;
    double t1, t2;
    double time_comm_raw, time_comm_cluster, time_comm_vel, time_comm_pres;
    time_comm_raw = 0;
    time_comm_cluster =0;
    time_comm_vel=0;
    time_comm_pres = 0;
   
    double time_comp=0;

    // if we define strict max_reader in dspaces configurations, we can read the same variable
    
    /*
     * vel and pressure data
     */
    char var_name_vel[STRING_LENGTH];
    char var_name_pres[STRING_LENGTH];
    sprintf(var_name_vel, "VEL");
    sprintf(var_name_pres, "PRES");

    // data layout
    uint64_t gdims_raw[3] = {POINTS_SIDE, POINTS_SIDE,1};
    dspaces_define_gdim(var_name_vel, 3,gdims_raw);
    dspaces_define_gdim(var_name_pres, 3,gdims_raw);

    size_t elem_size_vel = sizeof(float)*3;
    size_t elem_size_pres = sizeof(float);


    //x_min,y_min,z_min,x_max_y_max_z_max
    int bounds[6]={0};
    bounds[3] = dims[0]-1;
    bounds[4] = dims[1]-1;


    // prepare space
    float * vel_data = (float *)malloc(num_points* elem_size_vel);
    if(vel_data == NULL){
          perror("vel data allocated error");
          exit(-1);
      }

    // prepare space for pres
    float * pres_data = (float *)malloc(num_points*elem_size_pres);
    if(pres_data == NULL){
          perror("pres data allocated error");
          exit(-1);
      }
#ifdef INCLUDE_ML

    int num_region = NUM_REGION;

    char var_name_cluster[STRING_LENGTH];
    sprintf(var_name_cluster, "cluster");

    int bounds_cluster[6]={0};
    // x_min
    bounds_cluster[0] = 0;

    // x_max
    bounds_cluster[3] = num_region;
    printf("    number of regions %d\n", num_region);

    int num_elems_cluster = bounds_cluster[3] - bounds_cluster[0]+1;
    size_t elem_size_cluster = sizeof(float);

    float *buffer_cluster = (float *)malloc(num_elems_cluster*elem_size_cluster);
    if(buffer_cluster== NULL){
        perror("    allocate space cluster");
        exit(-1);
    }


    //cluster info will be small and 
    float * cluster_data = (float *)malloc(num_points* sizeof(float));
#endif

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

    // 1. read simulation  data from dataspces
    // this will get blocked until new data available

    get_common_buffer(timestep,bounds,rank, &gcomm, var_name_vel, (void **)&vel_data, elem_size_vel, &time_comm_vel);
    get_common_buffer(timestep,bounds,rank, &gcomm, var_name_pres, (void **)&pres_data, elem_size_pres, &time_comm_pres);
    time_comm_raw = time_comm_vel+ time_comm_pres;


    //update using vel and pres info, if there is more ranks I need to partition first
#ifndef INCLUDE_ML
    attributes.UpdateFields(vel_data, pres_data);
#else
    // updated on March 2
    // 1. get medoids from 

    // also get cluster buffer here
    //get_cluster_buffer(timestep, &num_region ,rank, &gcomm, var_name_cluster, &cluster_raw , &time_comm_cluster);


    // 2. read cluster  data from dataspces
    get_common_buffer(timestep, bounds_cluster, rank, &gcomm, var_name_cluster,(void **) &buffer_cluster, elem_size_cluster, &time_comm_cluster);

    
    t1 = MPI_Wtime();
    map_regions(buffer_cluster, cluster_data, dims, region_length);
    attributes.UpdateFields(vel_data, pres_data, cluster_data);
#endif

    // add one pipeline 
    // show together
#ifdef USE_CATALYST
    //FEAdaptor::CoProcess(grid, attributes, time, timestep, timestep == MAX_VERSION-1);
    FEAdaptor::CoProcess(grid, attributes, time, timestep, timestep == MAX_VERSION-1);
#endif

    t2 = MPI_Wtime();
    time_comp = t2-t1;
    // timer
        double global_time_comm_cluster;
        double global_time_comm_raw;
        double global_time_comp;

        MPI_Reduce(&time_comm_cluster, &global_time_comm_cluster, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comm_raw, &global_time_comm_raw, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comp, &global_time_comp, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);

        // Print the result
        if (rank == 0) {
          printf("%d Computation Total %lf avg %lf\n",timestep,  global_time_comp , global_time_comp/ (nprocs));
          printf("%d cluster Total %lf avg %lf\n",timestep,  global_time_comm_cluster , global_time_comm_cluster/ (nprocs));
          printf("%d divs Total %lf avg %lf\n",timestep,  global_time_comm_raw , global_time_comm_raw/ (nprocs));
        }
    }

    MPI_Barrier(gcomm);
    // reduce all comm_time
    
    // free buffer
    if(vel_data != NULL){
        free(vel_data);
        cout << "vel is freed" << endl;
    }
    if(pres_data != NULL){
        free(pres_data);
        cout << "pres is freed" << endl;
    }
#ifdef INCLUDE_ML
    if(cluster_data != NULL){
        free(cluster_data);
        cout << "cluster is freed" << endl;
    }
    if(buffer_cluster!= NULL){
        free(buffer_cluster);
        cout << "cluster is freed" << endl;
    }
#endif

#ifdef USE_CATALYST
  FEAdaptor::Finalize();
#endif

    dspaces_finalize();

    MPI_Barrier(gcomm);
    MPI_Finalize();

  return 0;
}

