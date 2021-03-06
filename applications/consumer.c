#include "consumer.h"
#define debug_1
//#define debug_2

// test segmentation fault
 int validate_regions(float *buffer_a, int region_memory_size){
   int i;
   for(i = 0; i< region_memory_size;i++){
   buffer_a[i] = 0.3;
   }
   return 1;
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

void cal_local_divs(float *buffer_all_regions, int region_length, int k_npdiv, int div_func, int *table, int  pair_index_l,int  pair_index_h,  int rank, float*divs_this_rank, double *p_time_used){
    int i;
    double t2, t3;
    // the index of the two pairs
    int a, b;
    float div = -1;
    float *buffer_a, *buffer_b;
    char msg[STRING_LENGTH];

    if(*p_time_used != 0){
        *p_time_used = 0;
    }


    t2 = MPI_Wtime();
    for(i = pair_index_l; i<= pair_index_h; i++){
        // get the region index of that pair

        //snprintf(msg, STRING_LENGTH,"i = %d, index_l = %d, index_h = %d\n", i, pair_index_l, pair_index_h);
        //my_message(msg, rank, LOG_WARNING);

        get_pair_index(table, i, &a, &b);

        //snprintf(msg, STRING_LENGTH,"try to access No.%d/%d pair, region %d and %d",i - pair_index_l, pair_index_h - pair_index_l +1, a, b);
        //my_message(msg, rank, LOG_WARNING);
#ifdef debug_2
        printf("\t\ttry to access No.%d/%d pair, region %d and %d\n",i - pair_index_l, pair_index_h - pair_index_l +1, a, b);
#endif

        buffer_a = buffer_all_regions + a*region_length*region_length*3;
        buffer_b = buffer_all_regions + b*region_length*region_length*3;


        div = get_divs( buffer_a , buffer_b, region_length, k_npdiv, div_func);


#ifdef debug_1
        //snprintf(msg, STRING_LENGTH,"got div of region %d and region %d: %.6f ",a, b, div);
        //my_message(msg, rank, LOG_WARNING);
#endif
        
        // save it into buffer first
        divs_this_rank[i - pair_index_l] = div;
        // record the time for communication and calculation
    }

    t3 = MPI_Wtime();

    *p_time_used = t3 -t2;
}

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

// start of dspaces configration

    // DataSpaces: Initalize and identify application
    // Usage: dspaces_init(num_peers, appid, Ptr to MPI comm, parameters)
    // Note: appid for get.c is 2 [for put.c, it was 1]
    //
    printf("trying init dspaces for %d process\n", nprocs);
    dspaces_init(nprocs, 2, &gcomm, NULL);

    char msg[STRING_LENGTH];

    /*
     * vel and pressure buffer
     */
    double time_used;

    char var_name_vel[STRING_LENGTH];
    char var_name_pres[STRING_LENGTH];
    sprintf(var_name_vel, "VEL");
    sprintf(var_name_pres, "PRES");

    uint64_t gdims_raw[2] = {POINTS_SIDE, POINTS_SIDE};
    dspaces_define_gdim(var_name_vel, 2,gdims_raw);
    dspaces_define_gdim(var_name_pres, 2,gdims_raw);


    int strip_size = POINTS_SIDE/nprocs;
    unsigned int dims[3] = {strip_size, POINTS_SIDE, 1};
    size_t elem_size_vel = sizeof(float)*3;
    size_t elem_size_pres = sizeof(float);
    unsigned int num_elems = dims[0]*dims[1]*dims[2];

    //x_min,y_min,z_min,x_max_y_max_z_max
    int bounds[6]={0};
    // xmin
    bounds[1] = strip_size*rank;
    // ymin = 0
    bounds[0] = 0;
    
    // xmax
    bounds[4] = strip_size*(rank+1)-1;
    // ymax
    bounds[3] = dims[1]-1;

    // prepare space
    float * vel_data = (float *)malloc(num_elems*elem_size_vel);
    if(vel_data == NULL){
          perror("vel data allocated error");
          exit(-1);
      }
    // prepare space for pres
    float * pres_data = (float *)malloc(num_elems*elem_size_pres);
    if(pres_data == NULL){
          perror("pres data allocated error");
          exit(-1);
      }

    // regions stripped to this rank
    int region_length = REGION_LENGTH;
    int num_region_row = dims[0]/region_length;
    int num_region_col =  dims[1]/region_length;
    int num_region = num_region_row*num_region_col; // this will be (201-1)/10 = 20

    int num_elems_region = num_region;
    size_t elem_size_region = (region_length)*(region_length)*3*sizeof(float);

    float *buffer_region = (float *)malloc(num_elems_region*elem_size_region);
    if(buffer_region== NULL){
        perror("    allocate space for striped  regions");
        exit(-1);
    }

    
    /*
     * aggregated sampled all data
     */
    // sampling related

    // it will operate on the same sample
    //char var_name_sample[STRING_LENGTH];
    //sprintf(var_name_sample_all, "sample");
    int nprocs_sim = (PROCS_PER_DIM*PROCS_PER_DIM);
    int sample_size = SAMPLE_SIZE;

    char var_name_sample[STRING_LENGTH];
    sprintf(var_name_sample, "sample");

    int bounds_sample_all[6]={0};
    // x_min
    bounds_sample_all[0] = 0;

    // x_max
    bounds_sample_all[3] = (nprocs_sim)*(sample_size) -1; 

    uint64_t gdims_sample[1] = {nprocs_sim*sample_size};
    dspaces_define_gdim(var_name_sample, 1,gdims_sample);


    int num_elems_sample_all = (bounds_sample_all[3]-bounds_sample_all[0] + 1);
    size_t elem_size_sample_all = (region_length)*(region_length)*3*sizeof(float);

    float *buffer_sample_all = (float *)malloc(num_elems_sample_all*elem_size_sample_all);
    if(buffer_sample_all== NULL){
        perror("    allocate space for sample all");
        exit(-1);
    }

    /*
     * divergence pairs
     * there will be num_region*num_region/2 divergences, each process will cal some of them
     */
    // how many clusters
    int cluster_k = NCLUSTERS;

	// how to compute divergence
	int k_npdiv = K_NPDIV;
    // use L2 divergence
    int div_func = 1;

    int num_tasks = num_elems_sample_all*(num_elems_sample_all-1)/2;

	// Each process will need to compute its DataSpace index
    int tasks_per_proc = num_tasks/nprocs;
    int tasks_left_over = num_tasks%nprocs;
    int pair_index_l, pair_index_h;
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

	int *table;
    //generate_lookup_table(num_region, &table);
    generate_lookup_table(num_elems_sample_all, &table);
    //
    sprintf(msg,"pair lookup table generated, I am responsible for P%d to P%d", pair_index_l, pair_index_h);
    my_message(msg, rank, LOG_CRITICAL);

    char var_name_divs[STRING_LENGTH];
    sprintf(var_name_divs, "divs");

    uint64_t gdims_divs[1] = {num_tasks};
    dspaces_define_gdim(var_name_divs, 1,gdims_divs);


    int bounds_divs[6]={0};

    // x_min
    bounds_divs[0] = pair_index_l;

    // x_max
    bounds_divs[3] = pair_index_h; 

    int num_elems_divs = (bounds_divs[3]-bounds_divs[0] + 1);
    size_t elem_size_divs = sizeof(float);

    float *buffer_divs = (float *)malloc(num_elems_divs*elem_size_divs);
    if(buffer_divs== NULL){
        perror("    allocate space for divs");
        exit(-1);
    }
    /*
     * medoids data
     * all process will receive the same medoids info
     */
    
    // generate 3 clusters in the end
    int medoids_size = NCLUSTERS;
    char var_name_medoids[STRING_LENGTH];
    sprintf(var_name_medoids, "medoids");

    uint64_t gdims_medoids[1] = {NCLUSTERS};
    dspaces_define_gdim(var_name_medoids, 1,gdims_medoids);


    int bounds_medoids[6]={0};
    // x_min
    bounds_medoids[0] = 0;

    // x_max
    bounds_medoids[3] = medoids_size -1; 

    int num_elems_medoids = bounds_medoids[3] - bounds_medoids[0] +1;
    // each element if just a id(from all sample) of the medoids
    size_t elem_size_medoids = sizeof(int);


    int *buffer_medoids = (int *)malloc(num_elems_medoids*elem_size_medoids);
    if(buffer_medoids== NULL){
        perror("    allocate space for medoids");
        exit(-1);
    }
    /*
     * clusterid data
     * all process will receive the same medoids info
     */
    char var_name_cluster[STRING_LENGTH];
    sprintf(var_name_cluster, "cluster");

    uint64_t gdims_cluster[1] = {nprocs*num_region};
    dspaces_define_gdim(var_name_cluster, 1, gdims_cluster);


    int bounds_cluster[6]={0};
    // x_min
    bounds_cluster[0] = rank*num_region;

    // x_max
    bounds_cluster[3] = (rank+1)*num_region -1;

    int num_elems_cluster = bounds_cluster[3] - bounds_cluster[0]+1;
    size_t elem_size_cluster = sizeof(float);

    float *buffer_cluster = (float *)malloc(num_elems_cluster*elem_size_cluster);
    if(buffer_cluster== NULL){
        perror("    allocate space cluster");
        exit(-1);
    }


// end of dspaces config

    // Name our data.

    // accumulated time for communication(read regions and write divs) and calculation(divs)
    //double time_comm_regions = 0;
    double t1,t2;
    // timer for end-to-end
    double t_start, t_end, time_latency;
    double time_comm_pres = 0;
    double time_comm_vel = 0;
    double time_comm_sample_all = 0;
    double time_comm_medoids =0;
    double time_comm_cluster= 0;
    double time_comm_divs= 0;
    double time_comp_divs =0;
    double time_comp_assign = 0;

    int timestep=0;
    while(timestep < MAX_VERSION){

        printf("********************timestep %d now start!\n",timestep);

        // do we need this barrier?
        MPI_Barrier(gcomm);

        // generated a list of random region ids
        // updated on March 2, three steps

        // 1. get stripped data, 
        get_common_buffer(timestep,2, bounds,rank, &gcomm, var_name_vel, (void **)&vel_data, elem_size_vel, &time_comm_vel);

        // this is actual start when enter dspaces read lock
        t_start = MPI_Wtime()- time_comm_vel;

        get_common_buffer(timestep,2, bounds,rank, &gcomm, var_name_pres, (void **)&pres_data, elem_size_pres, &time_comm_pres);

        // 3. divide into regions
        int num_region_2;
        divide(vel_data, dims,region_length,&num_region_2, buffer_region);
        printf("divide completed %d regions generated\n", num_region_2);
#ifdef debug_1
        int spacing = elem_size_region/sizeof(float);
        printf(" first data of all regions: %f %f %f \n", buffer_region[0], buffer_region[1], buffer_region[2]);
        printf(" last data of all regions: %f %f %f \n", buffer_region[spacing*num_region -3], buffer_region[spacing*num_region -2], buffer_region[spacing*num_region-1]);
#endif


        // 4. get aggregated sampled regions(dspaces get blocked if one applciation both writes and reds on the same variable)
        printf("start to gather global_sample\n");

        get_common_buffer(timestep,1,  bounds_sample_all, rank, &gcomm, var_name_sample, (void **)&buffer_sample_all, elem_size_region,  &time_comm_sample_all);
        printf("global_sample got\n");


#ifdef debug_1
        spacing = elem_size_region/sizeof(float);
        printf("\tfirst data of all samples: %f %f %f \n", buffer_sample_all[0], buffer_sample_all[1], buffer_sample_all[2]);
        printf("\tlast data of all samples: %f %f %f \n", buffer_sample_all[spacing*num_elems_sample_all -3], buffer_sample_all[spacing*num_elems_sample_all -2], buffer_sample_all[spacing*num_elems_sample_all-1]);
#endif

        // 5. calculate subset of divergence pairs based on sampled regions
        cal_local_divs(buffer_sample_all, region_length, k_npdiv,div_func,table, pair_index_l, pair_index_h, rank, buffer_divs, &time_comp_divs);
        printf("divergence calculated\n");

        // 6. send divergence to analyis
        put_common_buffer(timestep,1,  bounds_divs, rank, &gcomm, var_name_divs, (void **)&buffer_divs, elem_size_divs,  &time_comm_divs);
        printf("divergence sent\n");

        // wait here, for better performance those three steps should be put in next timestep
        // 6. get medoids info 
        get_common_buffer(timestep,1, bounds_medoids, rank, &gcomm, var_name_medoids, (void **)&buffer_medoids, elem_size_medoids, &time_comm_medoids);
        printf("medoids are %d %d %d\n",buffer_medoids[0], buffer_medoids[1], buffer_medoids[2]);

        // 7. assign clusterid
        t1 = MPI_Wtime();
        assign_clusterid(buffer_region, num_region, buffer_sample_all, num_elems_sample_all, buffer_medoids, cluster_k, region_length, k_npdiv, div_func,  buffer_cluster);
        t2 = MPI_Wtime();
        time_comp_assign = t2 - t1;
        printf("clusterid  assigned\n");

        // 7. put clusterid
        put_common_buffer(timestep,1,  bounds_cluster, rank, &gcomm, var_name_cluster, (void **)&buffer_cluster, elem_size_cluster, &time_comm_cluster);
        printf("clusterid  put to dspaces\n");


        // should wait until get all the divergence
        //sprintf(msg,"--time eclapsed for read regions// calculation // put divs:%.4lf %.4lf, %.4f", time_comm_vel, time_comp, time_comm_divs);

        sprintf(msg,"--has reached barrier and yeild div read lock to producer");
        my_message(msg, rank, LOG_CRITICAL);

        MPI_Barrier(gcomm);
        t_end = MPI_Wtime();
        time_latency = t_end-t_start;

        double time_comm_raw = time_comm_vel+ time_comm_pres;

        double global_time_comm_raw;
        double global_time_comm_sample;
        double global_time_comp_divs;
        double global_time_comm_divs;
        double global_time_comm_medoids;
        double global_time_comp_assign;
        double global_time_comm_cluster;

        MPI_Reduce(&time_comm_raw, &global_time_comm_raw, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comm_sample_all, &global_time_comm_sample, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comp_divs, &global_time_comp_divs, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comm_divs, &global_time_comm_divs, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comm_medoids, &global_time_comm_medoids, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comp_assign, &global_time_comp_assign, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comm_cluster, &global_time_comm_cluster, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);

        // Print the result
        if (rank == 0) {
          printf("%d comm raw Total %lf avg %lf\n",timestep,  global_time_comm_raw , global_time_comm_raw/ (nprocs));
          printf("%d comm sample Total %lf avg %lf\n",timestep,  global_time_comm_sample , global_time_comm_sample/ (nprocs));
          printf("%d comp divs Total %lf avg %lf\n",timestep,  global_time_comp_divs , global_time_comp_divs/ (nprocs));
          printf("%d comm divs Total %lf avg %lf\n",timestep,  global_time_comm_divs , global_time_comm_divs/ (nprocs));
          printf("%d comm medoids Total %lf avg %lf\n",timestep,  global_time_comm_medoids , global_time_comm_medoids/ (nprocs));
          printf("%d comp assign Total %lf avg %lf\n",timestep,  global_time_comp_assign , global_time_comp_assign/ (nprocs));
          printf("%d comm cluster Total %lf avg %lf\n",timestep,  global_time_comm_cluster , global_time_comm_cluster/ (nprocs));
          printf("%d all latency Total %lf\n",timestep,  time_latency );
        }
        timestep++;
    }

    // free buffer
    if(vel_data != NULL){
        free(vel_data);
        sprintf(msg,"--velocity buffer freed");
        my_message(msg, rank, LOG_CRITICAL);
    }
    if(pres_data != NULL){
        free(pres_data);
        sprintf(msg,"--pres buffer freed");
        my_message(msg, rank, LOG_CRITICAL);
    }

    if(table != NULL){
        free_lookup_table(table);
    }
    if(buffer_region != NULL){
        free(buffer_region);
        sprintf(msg,"-- buffer_region freed");
        my_message(msg, rank, LOG_CRITICAL);
    }

    if(buffer_sample_all != NULL){
        free(buffer_sample_all);
        sprintf(msg,"-- buffer_sample_all freed");
        my_message(msg, rank, LOG_CRITICAL);
    }
    if(buffer_divs != NULL){
        free(buffer_divs);
        sprintf(msg,"-- divs freed");
        my_message(msg, rank, LOG_CRITICAL);
    }
    if(buffer_medoids != NULL){
        free(buffer_medoids);
        sprintf(msg,"-- medoids freed");
        my_message(msg, rank, LOG_CRITICAL);
    }
    if(buffer_cluster != NULL){
        free(buffer_cluster);
        sprintf(msg,"-- buffer_cluster freed");
        my_message(msg, rank, LOG_CRITICAL);
    }

    sprintf(msg, "now finalize the dspaces and exit");
    my_message(msg, rank, LOG_CRITICAL);

    // DataSpaces: Finalize and clean up DS process
    dspaces_finalize();

    MPI_Barrier(gcomm);
    MPI_Finalize();
    return 0;
}

