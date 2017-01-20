#include "get_regions.h"
//#define debug_1

void fill_region_def(Region_Def *p_region_def){
    // some pre-definition
    p_region_def->region_length = REGION_LENGTH;
    p_region_def->side_num_region = (POINTS_SIDE - 1)/REGION_LENGTH;
    p_region_def-> num_region =p_region_def->side_num_region*p_region_def->side_num_region ;

    p_region_def->region_num_cell = (p_region_def->region_length+1)*(p_region_def->region_length+1);
    // attention here: this is used only when malloc
    p_region_def->region_memory_size = p_region_def->region_num_cell*3*sizeof(float);
}

void extract_region_def(Region_Def *p_region_def, int *p_region_length,int * p_side_num_region, int *p_num_region,int * p_region_num_cell, size_t *p_region_memory_size){
    *p_region_length = p_region_def->region_length;
    *p_side_num_region = p_region_def->side_num_region;
    *p_num_region = p_region_def->num_region ;
    *p_region_num_cell = p_region_def->region_num_cell;
    // attention here: this is used only when malloc
    *p_region_memory_size = p_region_def->region_memory_size;
}



void generate_lookup_table(int num_region, int **p_table){
    int i, j, count;

    int num_pair = num_region*(num_region)/2;
    int *table = (int *)malloc(num_pair*sizeof(int)*2);
    if(table == NULL){
        perror("malloc lookupdable");
        exit(-1);
    }

    count = 0;// pair index
    for(i = 1; i< num_region ; i++){
        for(j = 0; j < i; j++){
            table[2*count] = i;
            table[2*count + 1] = j;
            count+=1;
        }
    }

    *p_table = table;
}

void free_lookup_table(int *table){
    free(table);
}

void get_pair_index(int *table, int index_pair,int *a, int *b){
    *a = table[2*index_pair];
    *b = table[2*index_pair +1];
}

void cal_local_divs(float *buffer_all_regions, Region_Def * p_region_def, int k_npdiv, int div_func, int *table, int  pair_index_l,int  pair_index_h,  int rank, float** p_divs_this_rank, double *time_used){
    int i;
    double t2, t3;
    // the index of the two pairs
    int a, b;
    float div = -1;
    float *buffer_a, *buffer_b;
    char msg[STRING_LENGTH];

    if(time_used != 0){
        time_used =0;
    }

    // some pre-definition
    int region_length, side_num_region,num_region,region_num_cell;
    size_t region_memory_size;

    // extract those 
    extract_region_def(p_region_def, &region_length, &side_num_region,&num_region,&region_num_cell, &region_memory_size);

    // timer
    float* divs_this_rank;

    // prepare buffer for divs
    int size_div = (pair_index_h - pair_index_l+1)*sizeof(float);
    snprintf(msg, STRING_LENGTH,  "div buffer has size %d", size_div);
    my_message(msg, rank, LOG_WARNING);

    // freed in get_regions.c main function
    divs_this_rank = (float *)malloc(size_div);

    if(divs_this_rank == NULL){
        perror("malloc error for div buffer, now exit");
        exit(-1);
    }


    for(i = pair_index_l; i<= pair_index_h; i++){
        // get the region index of that pair

        snprintf(msg, STRING_LENGTH,"i = %d, index_l = %d, index_h = %d\n", i, pair_index_l, pair_index_h);
        my_message(msg, rank, LOG_WARNING);

        get_pair_index(table, i, &a, &b);

        snprintf(msg, STRING_LENGTH,"try to access No.%d/%d pair, region %d and %d",i - pair_index_l, pair_index_h - pair_index_l +1, a, b);
        my_message(msg, rank, LOG_WARNING);

        buffer_a = buffer_all_regions + a*region_num_cell*3;
        buffer_b = buffer_all_regions + b*region_num_cell*3;

        t2 = MPI_Wtime();

#ifdef debug
        int aa, bb;
        // validate the two regions
        aa =  validate_regions(buffer_a,region_memory_size);
        bb = validate_regions(buffer_b, region_memory_size);

        if(aa == 1 && bb == 1){
            sprintf(msg, "No.%d/%d pair, region %d(%p) and %d(%p) is validate",i - pair_index_l, pair_index_h - pair_index_l +1, a,(void*)buffer_a, b,(void *)buffer_b);
            my_message(msg, rank, LOG_VERB);
        }

        // we can get the divergence now!
        snprintf(msg, STRING_LENGTH,"buffer a:%p b:%p, region_length %d, k_npdiv %d, div_func %d", buffer_a, buffer_b, region_length, k_npdiv, div_func);
        my_message(msg, rank, LOG_WARNING);

#endif
        div = get_divs( buffer_a , buffer_b, region_length, k_npdiv, div_func);

        //sprintf(msg, "No.%d/%d pair, region %d and %d: %.3f",i - pair_index_l, pair_index_h - pair_index_l +1, a, b, div);
        //my_message(msg, rank);

        t3 = MPI_Wtime();


#ifdef debug_1
        snprintf(msg, STRING_LENGTH,"got div %.6f ", div);
        my_message(msg, rank, LOG_WARNING);
#endif
        // save it into buffer first
        divs_this_rank[i - pair_index_l] = div;

        snprintf(msg, STRING_LENGTH,"div written");
        my_message(msg, rank, LOG_WARNING);

        // record the time for communication and calculation
        *time_used += t3 -t2;
    }
    *p_divs_this_rank = divs_this_rank;
}

void get_region_buffer(int timestep, Region_Def *p_region_def, int  rank, MPI_Comm * p_gcomm,  float ** p_buffer_all_regions, double * time_used){
    char msg[STRING_LENGTH];
    double t1, t2;
    int ret_get = -1;

    // some pre-definition
    int region_length, side_num_region,num_region,region_num_cell;
    size_t region_memory_size;

    // extract those 
    extract_region_def(p_region_def, &region_length, &side_num_region,&num_region,&region_num_cell, &region_memory_size);

    // prepare to read regions from dataspaces
    uint64_t lb[3] = {0}, ub[3] = {0};
    lb[0] = 0;
    ub[0] = num_region - 1;

    // Define the dimensionality of the data to be received 
    int ndim = 3;

    char var_name_regions[STRING_LENGTH];
    sprintf(var_name_regions, "region_data");


    char lock_name_regions[STRING_LENGTH];
    snprintf(lock_name_regions, STRING_LENGTH, "region_lock_t_%d", timestep);


    // freed in get_regions.c main function
    float * buffer_all_regions = (float*)malloc(region_num_cell*sizeof(float)*3*num_region);
    if(buffer_all_regions == NULL){
        perror("malloc error for regions buffer, now exit");
        exit(-1);
    }

    sprintf(msg, "try to acquired the region read lock %s", lock_name_regions );
    my_message(msg, rank, LOG_WARNING);
    dspaces_lock_on_read(lock_name_regions, p_gcomm);

    sprintf(msg, "get the  the region read lock");
    my_message(msg, rank, LOG_WARNING);

    // read all regions in once
    t1 = MPI_Wtime();

    ret_get = dspaces_get(var_name_regions, timestep, region_memory_size, ndim, lb, ub, buffer_all_regions);

    t2 = MPI_Wtime();

    // now we can release region lock
    dspaces_unlock_on_read(lock_name_regions, p_gcomm);
    sprintf(msg, "release the region read lock");
    my_message(msg, rank, LOG_WARNING);

    if(ret_get != 0){
        perror("get all regions error, now exit");
        exit(-1);
    }else{
        sprintf(msg, "read %d regions from dspaces, each has %ld bytes", num_region, region_memory_size);
        my_message(msg, rank, LOG_WARNING);
    }

    *p_buffer_all_regions = buffer_all_regions;
    *time_used = t2-t1;
}

void put_divs_buffer(int timestep,int pair_index_l, int pair_index_h ,int num_tasks, int rank, MPI_Comm *p_gcomm, float ** p_divs_this_rank, double* time_used){

    char msg[STRING_LENGTH];
    double t1, t2;
    // return values for div dspaces operations
    int ret_put = -1;
    uint64_t lb_div[3] = {0}, ub_div[3] = {0};
    // save divergence in a 1-d array! this will save space

    int ndim_div = 3;
    char var_name_div[128];
    sprintf(var_name_div, "div_data");

    //uint64_t gdim_div[3] = {10000,1,1};
    uint64_t gdim_div[3] = {num_tasks,1,1};
    dspaces_define_gdim(var_name_div, 3,gdim_div);


    char lock_name_divs[STRING_LENGTH];
    snprintf(lock_name_divs, STRING_LENGTH, "div_lock_t_%d", timestep);


    lb_div[0] = pair_index_l, ub_div[0] = pair_index_h;

    // div variable operation
    sprintf(msg, "try to acquired the div write lock %s",lock_name_divs);
    my_message(msg, rank, LOG_WARNING);
    dspaces_lock_on_write(lock_name_divs, p_gcomm);

    sprintf(msg, "get div write lock");
    my_message(msg, rank, LOG_WARNING);

    //pay attention to the order of dimensions!
    //
    sprintf(msg, "write divergence to index%d~%d",  pair_index_l, pair_index_h);
    my_message(msg, rank, LOG_WARNING);

    t1 = MPI_Wtime();
    ret_put = dspaces_put(var_name_div, timestep, sizeof(float), ndim_div, lb_div, ub_div, *p_divs_this_rank);
    t2 = MPI_Wtime();

    // write div information into div variable
    dspaces_unlock_on_write(lock_name_divs, p_gcomm);

    sprintf(msg, " div write lock is released");
    my_message(msg, rank, LOG_WARNING);

    // how about the symmetric part?
    if(ret_put == 0){
        sprintf(msg, "divergence of %d pairs have saved into dspaces",pair_index_h - pair_index_l +1);
    }else{
        sprintf(msg, "ERROR when storing divergence of region");
    }
    my_message(msg, rank, LOG_CRITICAL);
    *time_used = t2-t1;
}

// test segmentation fault
/*
   int validate_regions(float *buffer_a, int region_memory_size){
   int i;
   for(i = 0; i< region_memory_size;i++){
   buffer_a[i] = 0.3;
   }
   return 1;
   }
   */

int main(int argc, char **argv)
{
    int nprocs, rank;
    MPI_Comm gcomm;

    char result_path[STRING_LENGTH]="";

    // output results into a folded specified with a slurm jobid
    if(argc == 2){
        strcpy(result_path, argv[1]);
    }



    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
    gcomm = MPI_COMM_WORLD;

    // DataSpaces: Initalize and identify application
    // Usage: dspaces_init(num_peers, appid, Ptr to MPI comm, parameters)
    // Note: appid for get.c is 2 [for put.c, it was 1]
    dspaces_init(nprocs, 2, &gcomm, NULL);


    // Name our data.


    char msg[STRING_LENGTH];


    // region definition
    // those parameters are obtained after dividing 
    // !! there should be communication between two application

    /*
       int region_length = 10;
       int num_region = 400;
       */

    // wrap all info for regions
    Region_Def region_def;
    Region_Def *p_region_def = &region_def;

    fill_region_def(p_region_def);

    // some pre-definition
    int region_length, side_num_region,num_region,region_num_cell;
    size_t region_memory_size;

    // extract those 
    extract_region_def(p_region_def, &region_length, &side_num_region,&num_region,&region_num_cell, &region_memory_size);

    // !!!feng: we should know how many regions in total so that we can asign the pairs
    // MPI_receive here?
    int num_tasks = num_region*(num_region-1)/2;

    sprintf(msg, "now div variable has dimention %d",num_tasks);
    my_message(msg, rank, LOG_CRITICAL);


    // Each process will need to compute its DataSpace index
    int tasks_per_proc = num_tasks/nprocs;
    int tasks_left_over = num_tasks%nprocs;
    int pair_index_l, pair_index_h;


    // accumulated time for communication(read regions and write divs) and calculation(divs)
    double time_comm_regions = 0;
    double time_cal = 0;
    double time_comm_divs = 0;

    // assign tasks to each rank
    if(rank<tasks_left_over){
        pair_index_l = rank*(tasks_per_proc+1);
        pair_index_h = pair_index_l + tasks_per_proc;
    }else{
        pair_index_l = rank*tasks_per_proc+tasks_left_over;
        pair_index_h = pair_index_l + tasks_per_proc-1;
    }
    if(rank == 0)
        printf("--rank%d:num_tasks=%d, tasks_per_proc=%d,tasks_left_over=%d\n", rank,num_tasks, tasks_per_proc, tasks_left_over);

    // generate pair lookup table;
    int *table;
    generate_lookup_table(num_region, &table);
    sprintf(msg,"pair lookup table generated, I am responsible for P%d to P%d", pair_index_l, pair_index_h);
    my_message(msg, rank, LOG_CRITICAL);
    //    printf("index_h = %d, index_l = %d", pair_index_h, pair_index_l);

    float *divs_this_rank;

    // every time read all the regions
    float *buffer_all_regions;

    // k-nearest neighbours
    int k_npdiv = K_NPDIV;

    // use L2 divergence
    int div_func = 1;

    int timestep=0;

    double time_used;

    /*
       snprintf(lock_name_regions, STRING_LENGTH, "region_lock_same");
       snprintf(lock_name_divs, STRING_LENGTH, "div_lock_same");
       */
    while(timestep < MAX_VERSION){

        if(rank == 0){
            sprintf(msg, "\n********************timestep %d now start!\n",timestep);
            my_message(msg, rank, LOG_WARNING);
        }

        // read all regions in once
        // note that the regions a rank need may be not continous
        get_region_buffer(timestep, &region_def, rank, &gcomm, &buffer_all_regions, &time_used);
        time_comm_regions = time_used;

        sprintf(msg, "buffer starts with %.3f %.3f %.3f %.3f %.3f ", buffer_all_regions[0],buffer_all_regions[1],buffer_all_regions[2],buffer_all_regions[3],buffer_all_regions[4]);
        my_message(msg, rank, LOG_CRITICAL);

        // calculate divergence
        cal_local_divs(buffer_all_regions,&region_def,k_npdiv, div_func, table, pair_index_l, pair_index_h, rank, &divs_this_rank, &time_used);
        time_cal = time_used;
        if(buffer_all_regions){
            free(buffer_all_regions);
        }

        // put divs into dataspaces
        put_divs_buffer(timestep, pair_index_l, pair_index_h, num_tasks, rank, &gcomm, &divs_this_rank, &time_used);
        time_comm_divs = time_used;
        if(divs_this_rank != NULL){
            free(divs_this_rank);
        }

        // should wait until get all the divergence
        sprintf(msg,"--time eclapsed for read regions// calculation // put divs:%.4lf %.4lf, %.4f", time_comm_regions, time_cal, time_comm_divs);
        my_message(msg, rank, LOG_CRITICAL);

        // do we need this barrier?
        MPI_Barrier(gcomm);

        sprintf(msg,"--has reached barrier and yeild div read lock to producer");
        my_message(msg, rank, LOG_CRITICAL);
        timestep++;
    }

    if(table != NULL);
    free_lookup_table(table);

    sprintf(msg, "now finalize the dspaces and exit");
    my_message(msg, rank, LOG_CRITICAL);

    // DataSpaces: Finalize and clean up DS process
    dspaces_finalize();

    MPI_Barrier(gcomm);
    MPI_Finalize();

    return 0;
}
