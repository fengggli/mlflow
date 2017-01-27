
#ifdef __cplusplus
extern "C" {
#endif

#include "dataspaces.h"
#include "common_utility.h"
#include "region_def.h"
#include "string.h"
#include <mpi.h>


// this will get all vel and pres data
void get_raw_buffer(int timestep, void *extra_info, int rank, MPI_Comm * p_gcomm,char * var_name_vel, float **p_buffer_vel, char * var_name_pres, float **p_buffer_pres,  double *p_time_used);

void put_raw_buffer(int timestep, Region_Def *p_region_def, int rank, MPI_Comm * p_gcomm, char *var_name_vel, float **p_buffer_vel, char *var_name_pres, float **p_buffer_pres,  double *p_time_used);
 

#ifdef __cplusplus
}
#endif
