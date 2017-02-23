
#ifdef __cplusplus
extern "C" {
#endif

#include "dataspaces.h"
#include "common_utility.h"
#include "region_def.h"
#include "string.h"
#include "mpi.h"

//#include <mpi.h>

/*
 * read raw velocity and pres buffer from dataspaces
 * INPUT:
 *  timestep:
 *      int
 *  extra_info
 *      dimension info
 *  rank
 *  p_gcomm, var_name_vel,_var_name_pres
 *      required by dsput
 *
 * OUTPUT
 *  p_buffer
 *      receiving buffer
 */

// this will get all vel and pres data
void get_raw_buffer(int timestep, int bounds[6], void *extra_info, int rank, MPI_Comm * p_gcomm,char * var_name_vel, float **p_buffer_vel, char * var_name_pres, float **p_buffer_pres,  double *p_time_used);

 

/*
 * put raw velocity and pres buffer into dataspaces
 * INPUT:
 *  timestep:
 *      int
 *  extra_info
 *      dimension info
 *  rank
 *  p_gcomm, var_name_vel,_var_name_pres
 *      required by dsput
 *
 *  p_buffer
 *      send buffer
 */
void put_raw_buffer(int timestep,int bounds[6], void *extra_info, int rank, MPI_Comm * p_gcomm, char *var_name_vel, float **p_buffer_vel, char *var_name_pres, float **p_buffer_pres,  double *p_time_used);

void get_cluster_buffer(int timestep, void *extra_info, int rank, MPI_Comm * p_gcomm,char * var_name_cluster, float **p_buffer_cluster,  double *p_time_used);

void put_cluster_buffer(int timestep, void * extra_info, int rank, MPI_Comm * p_gcomm, char *var_name_cluster, float **p_buffer_cluster,  double *p_time_used);

#ifdef __cplusplus
}
#endif
