#include "get_regions.h"
#define debug_1


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
    int err;
	int nprocs, rank;
	MPI_Comm gcomm;

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
	char var_name[STRING_LENGTH];
	sprintf(var_name, "region_data");


    char msg[STRING_LENGTH];
        
    
    // region definition
    // those parameters are obtained after dividing 
    // !! there should be communication between two application
    int i,j;
    
    /*
    int region_length = 10;
    int num_region = 400;
    */
    int region_length = REGION_LENGTH;
    int side_num_region = (POINTS_SIDE - 1)/REGION_LENGTH;
    int num_region =side_num_region*side_num_region ;
    


    // !!!feng: we should know how many regions in total so that we can asign the pairs
    // MPI_receive here?
    int num_tasks = num_region*(num_region-1)/2;

    sprintf(msg, "now div variable has dimention %d",num_tasks);
    my_message(msg, rank);


	// Each process will need to compute its DataSpace index
	int tasks_per_proc = num_tasks/nprocs;
	int tasks_left_over = num_tasks%nprocs;
    int pair_index_l, pair_index_h;
	int ds_lb_index, ds_ub_index;

    // timber
    double t1,t2,t3;

    // accumulated time for communication and calculation
    double time_comm = 0;
    double time_cal = 0;
	
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

    int region_num_cell = (region_length+1)*(region_length+1);
    // attention here: this is used only when malloc
    size_t  region_memory_size = region_num_cell*3*sizeof(float);

    // generate pair lookup table;
    int *table;
    generate_lookup_table(num_region, &table);
    sprintf(msg,"pair lookup table generated, I am responsible for P%d to P%d", pair_index_l, pair_index_h);
    my_message(msg, rank);
//    printf("index_h = %d, index_l = %d", pair_index_h, pair_index_l);

    // the index of the two pairs
    int a, b;
    float *buffer_a, *buffer_b;
    float *divs_this_rank;

    // Prepare LOWER and UPPER bound dimensions
    uint64_t lb[3] = {0}, ub[3] = {0};

    // Define the dimensionality of the data to be received 
    int ndim = 3;

    // parameter for second variable: div
    float div;
    uint64_t lb_div[3] = {0}, ub_div[3] = {0};
    // save divergence in a 1-d array! this will save space
    int ndim_div = 3;
	char var_name_div[128];
	sprintf(var_name_div, "div_data");

    
    //
    //uint64_t gdim_div[3] = {10000,1,1};
    uint64_t gdim_div[3] = {num_tasks ,1,1};
    dspaces_define_gdim(var_name_div, 3,gdim_div);

    // return values for regions dspaces operations
    int ret_get_0, ret_get_1;

    // return values for div dspaces operations
    int ret_put, ret_get;
    
    // every time read all the regions
    float *buffer_all_regions;

    // k-nearest neighbours
    int k_npdiv = K_NPDIV;

    // use L2 divergence
    int div_func = 1;

    int timestep=0;
    char lock_name_regions[STRING_LENGTH];
    char lock_name_divs[STRING_LENGTH];

    /*
    snprintf(lock_name_regions, STRING_LENGTH, "region_lock_same");
    snprintf(lock_name_divs, STRING_LENGTH, "div_lock_same");
        */
    while(timestep < MAX_VERSION){
        snprintf(lock_name_regions, STRING_LENGTH, "region_lock_t_%d", timestep);
        snprintf(lock_name_divs, STRING_LENGTH, "div_lock_t_%d", timestep);

        if(rank == 0){
            sprintf(msg, "\n********************timestep %d now start!\n",timestep);
            my_message(msg, rank);
        }
        // I should wait after all regions are placed into dspaces

        // return values for regions dspaces operations
        ret_get_0 = -1;
        ret_get_1 = -1;

        // return values for div dspaces operations
        ret_put = -1;
        ret_get = -1;


        // prepare buffer for divs
        int size_div = (pair_index_h - pair_index_l+1)*sizeof(float);
        sprintf(msg, "div buffer has size %d", size_div);
        my_message(msg, rank);

        divs_this_rank = (float *)malloc(size_div);
        if(divs_this_rank == NULL){
            perror("malloc error for div buffer, now exit");
            exit(-1);
        }


        // prepare to read regions from dataspaces
        lb[0] = 0;
        ub[0] = num_region - 1;
        

        // this is the to prefetch all regions
        // note that the regions a rank need may be not continous
        buffer_all_regions = (float*)malloc(region_num_cell*sizeof(float)*3*num_region);
        if(buffer_all_regions == NULL){
            perror("malloc error for regions buffer, now exit");
            exit(-1);
        }

        sprintf(msg, "try to acquired the region read lock %s", lock_name_regions );
        my_message(msg, rank);
        dspaces_lock_on_read(lock_name_regions, &gcomm);

        sprintf(msg, "get the  the region read lock");
        my_message(msg, rank);

        // read all regions in once
        t1 = MPI_Wtime();

        ret_get_0 = dspaces_get(var_name, timestep, region_memory_size, ndim, lb, ub, buffer_all_regions);

        t2 = MPI_Wtime();

        // now we can release region lock
        dspaces_unlock_on_read(lock_name_regions, &gcomm);
        sprintf(msg, "release the region read lock");
        my_message(msg, rank);

        if(ret_get_0 != 0){
            perror("get all regions error, now exit");
            exit(-1);
        }else{
            sprintf(msg, "read %d regions from dspaces, each has %ld bytes", num_region, region_memory_size);
            my_message(msg, rank);
        }

        // read all regions in once
        time_comm += t2 -t1;

        
        for(i = pair_index_l; i<= pair_index_h; i++){

            // get the region index of that pair
            get_pair_index(table, i, &a, &b);

#ifdef debug_1
            sprintf(msg, "try to access No.%d/%d pair, region %d and %d",i - pair_index_l, pair_index_h - pair_index_l +1, a, b);
            my_message(msg, rank);
#endif

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
                my_message(msg, rank);
            }
#endif
                
            // we can get the divergence now!
            div = get_divs( buffer_a , buffer_b, region_length, k_npdiv, div_func);

            sprintf(msg, "No.%d/%d pair, region %d and %d: %.3f",i - pair_index_l, pair_index_h - pair_index_l +1, a, b, div);
            my_message(msg, rank);

            t3 = MPI_Wtime();

            // save it into buffer first
            divs_this_rank[i - pair_index_l] = div;

            // record the time for communication and calculation
            time_cal += t3 -t2;
        }


        free(buffer_all_regions);


        lb_div[0] = pair_index_l, ub_div[0] = pair_index_h;

        // div variable operation
        sprintf(msg, "try to acquired the div write lock %s",lock_name_divs);
        my_message(msg, rank);
        dspaces_lock_on_write(lock_name_divs, &gcomm);

        sprintf(msg, "get div write lock");
        my_message(msg, rank);

        //pay attention to the order of dimensions!
        //
#ifdef debug_1
        sprintf(msg, "write divergence to index%d~%d",  pair_index_l, pair_index_h);
        my_message(msg, rank);
#endif

        t1 = MPI_Wtime();
        ret_put = dspaces_put(var_name_div, timestep, sizeof(float), ndim_div, lb_div, ub_div, divs_this_rank);
        t2 = MPI_Wtime();

        // write div information into div variable
        dspaces_unlock_on_write(lock_name_divs, &gcomm);

        sprintf(msg, " div write lock is released");
        my_message(msg, rank);

        // how about the symmetric part?
        if(ret_put == 0){
            sprintf(msg, "divergence of %d pairs have saved into dspaces",pair_index_h - pair_index_l +1);
        }else{
            sprintf(msg, "ERROR when storing divergence of region");
        }
        my_message(msg, rank);


        // now the divs buffer can be freed
        if(divs_this_rank != NULL){
            free(divs_this_rank);
        }

        // should wait until get all the divergence
        sprintf(msg,"--time eclapsed for read regions// calculation // put divs:%.4lf %.4lf, %.4f", time_comm, time_cal, t2 -t1);
        my_message(msg, rank);

        // do we need this barrier?
        MPI_Barrier(gcomm);

        sprintf(msg,"--has reached barrier and yeild div read lock to producer");
        my_message(msg, rank);

        timestep++;
    }
	
    if(table != NULL);
    free_lookup_table(table);

    sprintf(msg, "now finalize the dspaces and exit");
    my_message(msg, rank);

	// DataSpaces: Finalize and clean up DS process
	dspaces_finalize();

	MPI_Barrier(gcomm);
	MPI_Finalize();

	return 0;
}
