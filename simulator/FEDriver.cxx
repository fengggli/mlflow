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

// this will map cluster id buffer into the same size of pressure buffer
// this function is similar to divide function
// dim is the how many points in each side// 201
void map_regions(float *cluster_raw, float *cluster_data, const unsigned int dim, const unsigned int region_length){
    int p, q, ii, jj, index_x, index_y, linear_index;

    int l = region_length;
    printf("dim == %d, region_length = %d\n", dim, region_length);

    float clusterid;

    if((dim - 1)%l != 0){
        printf("not pefect division, try different region size of datacut size\n");
        exit(-1);
    }
    int side_num_region = (dim -1)/l; // this will be (201-1)/10 = 20
    for(p = 0; p < side_num_region; p++){
        for( q = 0; q < side_num_region; q++){
            clusterid = *(cluster_raw + side_num_region*p + q);
            // for each region
            // this region will have corner:
            //  (pl, ql) , (pl, (q+1)l)
            //  ((p+1)l,ql), ((p+1)l,(q+1)l)
            //  also the region all have the center (pl+l/2, ql+l/2)
            for(ii = 0; ii < l+1; ii++){
                for(jj = 0; jj < l+1; jj++){
                    // for each point inside the region
                    index_x = p*l+ii;
                    index_y = q*l+jj;

                    // mapped the logic address into linear address
                    // will be overlap but it's okay
                    linear_index = index_x*dim + index_y;
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
  unsigned int num_regions = NUM_REGION;
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
    double time_comm_raw, time_comm_cluster;
    time_comm_raw = 0;
    time_comm_cluster =0;
    double time_comp=0;

    // if we define strict max_reader in dspaces configurations, we can read the same variable
    char var_name_vel[STRING_LENGTH];
    char var_name_pres[STRING_LENGTH];
    char var_name_cluster[STRING_LENGTH];
    sprintf(var_name_vel, "VEL");
    sprintf(var_name_pres, "PRES");
    sprintf(var_name_cluster, "CLUSTER");
    // data layout
    uint64_t gdims_raw[3] = {POINTS_SIDE, POINTS_SIDE,1};
    dspaces_define_gdim(var_name_vel, 3,gdims_raw);
    dspaces_define_gdim(var_name_pres, 3,gdims_raw);


    //x_min,y_min,z_min,x_max_y_max_z_max
    int bounds[6]={0};
    bounds[3] = dims[0]-1;
    bounds[4] = dims[1]-1;


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
#ifdef INCLUDE_ML

    // this may be space efficiency but will work 
    float * cluster_raw = (float *)malloc(num_regions* sizeof(float));
    if(cluster_raw == NULL){
          perror("cluster data allocated error");
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

    // read data from dataspces
    // this will get blocked until new data available
    get_raw_buffer(timestep,bounds, NULL ,rank, &gcomm, var_name_vel, &vel_data, var_name_pres,  &pres_data, &time_comm_raw);


    //update using vel and pres info, if there is more ranks I need to partition first
#ifndef INCLUDE_ML
    attributes.UpdateFields(vel_data, pres_data);
#else
    // updated on March 2
    // 1. get medoids from 
    int num_region = NUM_REGION;

    // also get cluster buffer here
    get_cluster_buffer(timestep, &num_region ,rank, &gcomm, var_name_cluster, &cluster_raw , &time_comm_cluster);

    
    t1 = MPI_Wtime();
    map_regions(cluster_raw, cluster_data, dims[0], region_length);
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
    if(cluster_raw != NULL){
        free(cluster_raw);
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

