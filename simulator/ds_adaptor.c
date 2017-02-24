#include "ds_adaptor.h"

#define debug_1
// this will get all vel and pres data
void get_raw_buffer(int timestep, int bounds[6], void *extra_info, int rank, MPI_Comm * p_gcomm,char * var_name_vel, float **p_buffer_vel, char * var_name_pres, float **p_buffer_pres,  double *p_time_used){
    // how many number of elements are actually written
    int num_elems;
    char msg[STRING_LENGTH];
    double t1, t2;
    int ret_get = -1;

    float *vel_data = *p_buffer_vel;
    float *pres_data = *p_buffer_pres;

    int num_points;

    // user provided the area size
    if(extra_info != NULL){
        //printf("no extra info required\n");
        //exit(-1);
        num_points = *(int*)extra_info;
    }
    else{

        // data layout
        int dims[3] = {POINTS_SIDE, POINTS_SIDE,1};
        num_points = dims[0]*dims[1]*dims[2];
    }

    size_t elem_size_vel = sizeof(float)*3;
    size_t elem_size_pres = sizeof(float);
    
    uint64_t lb[3] = {0}, ub[3] = {0};
    /*
    lb[0] = 0;
    ub[0] = num_points - 1;
    */
    lb[0] = bounds[1];
    lb[1] = bounds[0];
    //y
    ub[0] = bounds[4];
    //x
    ub[1] = bounds[3];

    // Define the dimensionality of the data to be received 
    int ndim = 3;

    char lock_name_vel[STRING_LENGTH];
#ifdef USE_SAME_LOCK
    snprintf(lock_name_vel, STRING_LENGTH, "%s_lock", var_name_vel);
#else
    snprintf(lock_name_vel, STRING_LENGTH, "%s_lock_t_%d",var_name_vel, timestep);
#endif

    char lock_name_pres[STRING_LENGTH];
#ifdef USE_SAME_LOCK
    snprintf(lock_name_pres, STRING_LENGTH, "%s_lock", var_name_pres);
#else
    snprintf(lock_name_pres, STRING_LENGTH, "%s_lock_t_%d", var_name_pres, timestep);
#endif
    sprintf(msg, "try to acquired the vel read lock %s", lock_name_vel );
    my_message(msg, rank, LOG_WARNING);
    dspaces_lock_on_read(lock_name_vel, p_gcomm);

    sprintf(msg, "get the  the vel read lock");
    my_message(msg, rank, LOG_WARNING);

    // read all regions in once
    t1 = MPI_Wtime();

    ret_get = dspaces_get(var_name_vel, timestep, elem_size_vel, ndim, lb, ub, vel_data);

#ifdef debug_1

    printf("matrix content:");
   // print_matrix(vel_data, 81, );

    snprintf(msg, STRING_LENGTH, "%s, var name is %s, timstep: %d, elem_size_vel = %d, ndim =%d. lb=[%d, %d, %d], hb=[%d, %d, %d] \n", __func__, var_name_vel, timestep, elem_size_vel, ndim, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2]);
    my_message(msg, rank, LOG_WARNING);

    snprintf(msg, STRING_LENGTH,"num_elem %d, first data %f %f %f, last data %f %f %f\n", num_points,vel_data[0],vel_data[1],vel_data[2],vel_data[3*num_points-3], vel_data[3*num_points-2],vel_data[3*num_points-1]);
    my_message(msg, rank, LOG_WARNING);
#endif

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
}
// this will get all vel and pres data
// bounds is [x,y,z] 
// but  lb and ub will {fasted, lowest}, for 2d,{y,x,0}
void put_raw_buffer(int timestep,int bounds[6], void * extra_info, int rank, MPI_Comm * p_gcomm, char *var_name_vel, float **p_buffer_vel, char *var_name_pres, float **p_buffer_pres,  double *p_time_used){
    char msg[STRING_LENGTH];


    double t1, t2;
    int ret_put = -1;

    int num_points;

    if(extra_info != NULL){
        //printf("no extra info required\n");
        //exit(-1);
        num_points = *(int*)extra_info;
    }
    else{
        printf("dim info required\n");
        exit(-1);
    }

    size_t elem_size_vel = sizeof(float)*3;
    size_t elem_size_pres = sizeof(float);
    
    // prepare to write regions to dataspaces
    uint64_t lb[3] = {0}, ub[3] = {0};
    /*
    lb[0] = 0;
    ub[0] = num_points - 1;
    */
    lb[0] = bounds[1];
    lb[1] = bounds[0];
    ub[0] = bounds[4];
    ub[1] = bounds[3];

    // Define the dimensionality of the data to be received 
    int ndim = 3;

    
    char lock_name_vel[STRING_LENGTH];
#ifdef USE_SAME_LOCK
    snprintf(lock_name_vel, STRING_LENGTH, "%s_lock",var_name_vel );
#else
    snprintf(lock_name_vel, STRING_LENGTH, "%s_lock_t_%d",var_name_vel,  timestep);
#endif
    //snprintf(lock_name_vel, STRING_LENGTH, "vel_lock");

    char lock_name_pres[STRING_LENGTH];
#ifdef USE_SAME_LOCK
    snprintf(lock_name_pres, STRING_LENGTH, "%s_lock", var_name_pres);
#else
    snprintf(lock_name_pres, STRING_LENGTH, "%s_lock_t_%d", var_name_pres, timestep);
#endif

    
    sprintf(msg, "try to acquired the vel write lock %s", lock_name_vel );
    my_message(msg, rank, LOG_WARNING);
    dspaces_lock_on_write(lock_name_vel, p_gcomm);

    sprintf(msg, "get the  the vel write lock");
    my_message(msg, rank, LOG_WARNING);

    // write all regions in once
    t1 = MPI_Wtime();

    ret_put = dspaces_put(var_name_vel, timestep, elem_size_vel, ndim, lb, ub, *p_buffer_vel);

    dspaces_put_sync();

    t2 = MPI_Wtime();

#ifdef debug_1


    float * vel_data = *p_buffer_vel;

    printf("matrix content:");
    //print_matrix(vel_data, 40, 40);

    snprintf(msg, STRING_LENGTH, "%s, var name is %s, timstep: %d, elem_size_vel = %d, ndim =%d. lb=[%d, %d, %d], hb=[%d, %d, %d] \n", __func__, var_name_vel, timestep, elem_size_vel, ndim, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2]);
    my_message(msg, rank, LOG_WARNING);

    snprintf(msg, STRING_LENGTH,"num_elem %d, first data %f %f %f, last data %f %f %f\n", num_points,vel_data[0],vel_data[1],vel_data[2],vel_data[3*num_points-3], vel_data[3*num_points-2],vel_data[3*num_points-1]);
    my_message(msg, rank, LOG_WARNING);
#endif


    // now we can release region lock
    dspaces_unlock_on_write(lock_name_vel, p_gcomm);
    sprintf(msg, "release the vel write lock");
    my_message(msg, rank, LOG_WARNING);

    if(ret_put != 0){
        perror("put all vel error, now exit");
        printf("error number %d \n", ret_put);
        exit(-1);
    }else{
        sprintf(msg, "write %d vel to dspaces, each has %ld bytes", num_points, elem_size_vel);
        my_message(msg, rank, LOG_WARNING);
    }


    *p_time_used = t2-t1;
    
    // do the same for pres data
    sprintf(msg, "try to acquired the pres write lock %s", lock_name_pres );
    my_message(msg, rank, LOG_WARNING);
    dspaces_lock_on_write(lock_name_pres, p_gcomm);

    sprintf(msg, "get the  the pres write lock");
    my_message(msg, rank, LOG_WARNING);

    // write all regions in once
    t1 = MPI_Wtime();

    ret_put = dspaces_put(var_name_pres, timestep, elem_size_pres, ndim, lb, ub, *p_buffer_pres);

    dspaces_put_sync();

    t2 = MPI_Wtime();

    // now we can release region lock
    dspaces_unlock_on_write(lock_name_pres, p_gcomm);
    sprintf(msg, "release the pres write lock");
    my_message(msg, rank, LOG_WARNING);

    if(ret_put != 0){
        perror("put all pres error, now exit");
        printf("error number %d \n", ret_put);
        exit(-1);
    }else{
        sprintf(msg, "write %d pres to dspaces, each has %ld bytes", num_points, elem_size_pres);
        my_message(msg, rank, LOG_WARNING);
    }

    *p_time_used += t2-t1;
}

void get_cluster_buffer(int timestep, void *extra_info, int rank, MPI_Comm * p_gcomm,char * var_name_cluster, float **p_buffer_cluster,  double *p_time_used){
    char msg[STRING_LENGTH];
    double t1, t2;
    int ret_get = -1;

    float *cluster_data = *p_buffer_cluster;

    if(extra_info == NULL){
        printf("region size required\n");
        exit(-1);
    }

    int    num_regions = *(int*)extra_info;


    size_t elem_size_cluster = sizeof(float);
    
    // prepare to read regions from dataspaces
    uint64_t lb[3] = {0}, ub[3] = {0};
    lb[0] = 0;
    ub[0] = num_regions - 1;

    // Define the dimensionality of the data to be received 
    int ndim = 3;

    char lock_name_cluster[STRING_LENGTH];

    // assumer that ML won't be too fast
    //snprintf(lock_name_cluster, STRING_LENGTH, "cluster_lock");
    snprintf(lock_name_cluster, STRING_LENGTH, "cluster_lock_t_%d", timestep);
    // do the same for cluster data
    sprintf(msg, "try to acquired the cluster read lock %s", lock_name_cluster );
    my_message(msg, rank, LOG_WARNING);
    dspaces_lock_on_read(lock_name_cluster, p_gcomm);

    sprintf(msg, "get the  the cluster read lock");
    my_message(msg, rank, LOG_WARNING);

    // read all regions in once
    t1 = MPI_Wtime();

    ret_get = dspaces_get(var_name_cluster, timestep, elem_size_cluster, ndim, lb, ub, cluster_data);

    t2 = MPI_Wtime();

    // now we can release region lock
    dspaces_unlock_on_read(lock_name_cluster, p_gcomm);
    sprintf(msg, "release the cluster read lock");
    my_message(msg, rank, LOG_WARNING);

    if(ret_get != 0){
        perror("get all cluster error, now exit");
        printf("error number %d \n", ret_get);
        exit(-1);
    }else{
        sprintf(msg, "read %d cluster from dspaces, each has %zu bytes", num_regions, elem_size_cluster);
        my_message(msg, rank, LOG_WARNING);
    }

    *p_buffer_cluster = cluster_data;
    *p_time_used += t2-t1;
}

// this will put cluster id info into dspaces
void put_cluster_buffer(int timestep, void * extra_info, int rank, MPI_Comm * p_gcomm, char *var_name_cluster, float **p_buffer_cluster,  double *p_time_used){
    char msg[STRING_LENGTH];


    double t1, t2;
    int ret_put = -1;

    
    if(extra_info == NULL){
        printf("region size required\n");
        exit(-1);
    }

    int    num_regions = *(int*)extra_info;


    size_t elem_size_cluster = sizeof(float);
    
    // prepare to write regions to dataspaces
    uint64_t lb[3] = {0}, ub[3] = {0};
    lb[0] = 0;
    ub[0] = num_regions - 1;

    // Define the dimensionality of the data to be received 
    int ndim = 3;

    char lock_name_cluster[STRING_LENGTH];
    //snprintf(lock_name_cluster, STRING_LENGTH, "cluster_lock");
    //assumer that ML won't be faster
    snprintf(lock_name_cluster, STRING_LENGTH, "cluster_lock_t_%d", timestep);

    
    sprintf(msg, "try to acquired the cluster write lock %s", lock_name_cluster );
    my_message(msg, rank, LOG_WARNING);
    dspaces_lock_on_write(lock_name_cluster, p_gcomm);

    sprintf(msg, "get the  the cluster write lock");
    my_message(msg, rank, LOG_WARNING);

    // write all regions in once
    t1 = MPI_Wtime();

    ret_put = dspaces_put(var_name_cluster, timestep, elem_size_cluster, ndim, lb, ub, *p_buffer_cluster);

    dspaces_put_sync();

    t2 = MPI_Wtime();

    // now we can release region lock
    dspaces_unlock_on_write(lock_name_cluster, p_gcomm);
    sprintf(msg, "release the cluster write lock");
    my_message(msg, rank, LOG_WARNING);

    if(ret_put != 0){
        perror("put all cluster error, now exit");
        printf("error number %d \n", ret_put);
        exit(-1);
    }else{
        sprintf(msg, "write %d cluster to dspaces, each has %ld bytes", num_regions, elem_size_cluster);
        my_message(msg, rank, LOG_WARNING);
    }

    *p_time_used += t2-t1;
}
