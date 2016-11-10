#include "get_regions.h"


void my_message(char *msg, int rank){
    printf("**rank %d: %s\n", rank, msg);
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
	char var_name[128];
	sprintf(var_name, "region_data");


    char msg[80];
        
    
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
    

    // how many neareast neighbours: 3
    int k = 3;

    // !!!feng: we should know how many regions in total so that we can asign the pairs
    // MPI_receive here?
    int num_tasks = num_region*(num_region-1)/2;

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
	
	if(rank<tasks_left_over){
		pair_index_l = rank*(tasks_per_proc+1);
		pair_index_h = pair_index_l + tasks_per_proc;
	}else{
		pair_index_l = rank*tasks_per_proc+tasks_left_over;
		pair_index_h = ds_lb_index + tasks_per_proc-1;
	}

    int region_memory_size = (region_length+1)*(region_length+1)*3*sizeof(float);

    // generate pair lookup table;
    int *table;
    generate_lookup_table(num_region, &table);
    my_message("pair lookup table generated", rank);


    int timestep=1;
    // I should wait after all regions are placed into dspaces

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
    uint64_t gdim_div[3] = {10000,1,1};
    dspaces_define_gdim(var_name_div, 3,gdim_div);

    // use L2 divergence
    int div_func = 1;

    // return values for regions dspaces operations
    int ret_get_0 = -1;
    int ret_get_1 = -1;

    // return values for div dspaces operations
    int ret_put = -1;
    int ret_get = -1;


    // prepare buffer for regions
    buffer_a = (float *)malloc(region_memory_size);
    buffer_b = (float *)malloc(region_memory_size);

    // prepare buffer for divs
    divs_this_rank = (float *)malloc((pair_index_h - pair_index_l+1)*sizeof(float));
    if(buffer_a == NULL || buffer_b == NULL || divs_this_rank == NULL){
        perror("malloc error");
        exit(-1);
    }


    // put the divergence into divergence matrix
    
    sprintf(msg, "try to acquired the region read lock");
    my_message(msg, rank);

	dspaces_lock_on_read("region_lock", &gcomm);

    sprintf(msg, "get the  the region read lock");
    my_message(msg, rank);



    for(i = pair_index_l; i<= pair_index_h; i++){

        // get the region index of that pair
        get_pair_index(table, i, &a, &b);

            
        // get buffer for left region
        lb[0]= a;
        ub[0] =a;


        t1 = MPI_Wtime();
        ret_get_0 = dspaces_get(var_name, 1, region_memory_size, ndim, lb, ub, buffer_a);

        // get buffer for right region
        lb[0]= b;
        ub[0] =b;
        ret_get_1 = dspaces_get(var_name, 1, region_memory_size, ndim, lb, ub, buffer_b);
        t2 = MPI_Wtime();

        if(ret_get_0 == 0 && ret_get_1 == 0){
            //sprintf(msg, "get regions %d and %d success",a, b);
        }else{
            sprintf(msg, "ERROR when getting regions %d and %d",a, b);
            my_message(msg, rank);
        }

        // we can get the divergence now!
        div = get_divs( buffer_a , buffer_b, region_length, k, div_func);

        t3 = MPI_Wtime();

        // save it into buffer first
        divs_this_rank[i - pair_index_l] = div;

        // record the time for communication and calculation
        time_comm += t2 -t1;
        time_cal += t3 -t2;
    }

    // now we can release velocity lock
	dspaces_unlock_on_read("region_lock", &gcomm);


    t2 = MPI_Wtime();

    free(buffer_a);
    free(buffer_b);


    dspaces_lock_on_write("div_lock", &gcomm);
    // write div information into div variable

    //!!!! pay attention to the order of dimensions!
    lb_div[0] = pair_index_l, ub_div[0] = pair_index_h;

    t1 = MPI_Wtime();
    ret_put = dspaces_put(var_name_div, timestep, sizeof(float), ndim_div, lb_div, ub_div, divs_this_rank);
    t2 = MPI_Wtime();

    // how about the symmetric part?

    if(ret_put == 0){
        sprintf(msg, "divergence of %d pairs have saved into dspaces",pair_index_h - pair_index_l +1);
    }else{
        sprintf(msg, "ERROR when storing divergence of region");
    }
    my_message(msg, rank);


    dspaces_unlock_on_write("div_lock", &gcomm);

    // now the divs buffer can be freed
    if(divs_this_rank != NULL){
        free(divs_this_rank);
    }

    // should wait until get all the divergence
    sprintf(msg,"--time eclapsed for read regions// calculation // put divs:%.4lf %.4lf, %.4f", time comm, time_cal, t2 -t1);
    my_message(msg, rank);
    MPI_Barrier(gcomm);

    sprintf(msg,"--has reached barrier and try to acquire div read lock");
    my_message(msg, rank);

    // every rank need at least acquire the lock
    dspaces_lock_on_read("div_lock", &gcomm);
    sprintf(msg, "acquired div read lock");
    my_message(msg, rank);
    
	// Report data to user
	if(rank==0){
        // get the clustering done here


        //reconstruct the divergence matrix, this is a ragged matrix, only le 
        //see the distancematrix function in cluster.c 
        double **matrix = malloc(num_region*sizeof(double *)); 
        if(matrix == NULL){
            perror("malloc for div matrix");
            exit(-1);
        }
        matrix[0] = NULL;
        for(i = 1; i< num_region; i++){
            matrix[i] = malloc(i*sizeof(double));
            if(matrix[i] == NULL) break;
        }
        if(i< num_region){
            for(j = 0 ;j< i; j++)
                free(matrix[j]);
            perror("malloc 2 for div matrix");
            exit(-2);
        }


        int count = 0;
        int error_flag = 0;

        sprintf(msg, "divergence matrix constructed, now try to fill all the values");
        my_message(msg, rank);

        t1 = MPI_Wtime();
        for(i = 1; i < num_region; i++){
            for(j = 0; j < i;j++){
                // its the same order as when its saved
                lb_div[0] = count;
                ub_div[0] = count;
		        ret_get = dspaces_get(var_name_div, timestep, sizeof(float), ndim_div, lb_div, ub_div, &(matrix[i][j]));
                if(ret_get != 0){
                    error_flag  = 1;
                    break;
                }
                count +=1;
            }
        }

        t2 = MPI_Wtime();


        sprintf(msg, "divergence matrix filled in %.3f s time", t2-t1);
        my_message(msg, rank);

        if(error_flag == 1){
            sprintf(msg, "ERROR when read divergence from Dspaces");
            my_message(msg, rank);
        }
        // if everything is fine we start clustering
        else{

            // clustering parameters
            int nclusters = 3;
            int npass = 100;
            int clusterid[num_region];
            double error;
            int ifound;


            sprintf(msg, "start clustering");
            my_message(msg, rank);

            t1 = MPI_Wtime();
            kmedoids(nclusters, num_region, matrix, npass, clusterid, &error, &ifound);

            t2 = MPI_Wtime();

            sprintf(msg, "finished clustering in %.3lf s  time", t2 -t1);
            my_message(msg, rank);

            // save cluster results into file
            char *output_path = "data/clusterid_201_1.txt";
            FILE * f_clusterid = fopen(output_path, "w");
            if(f_clusterid == NULL){
                perror("file open error");                                                                                                                                                            
                exit(-1);
            }

            for(i = 0; i < num_region; i++){
                fprintf(f_clusterid, "%d\n",clusterid[i]);
            }
            
            fclose(f_clusterid);
        }

        // free the divergence buffer
        for(i = 1; i< num_region; i++){
            if(matrix[i] != NULL) free(matrix[i]);
        }
        if(matrix != NULL) {
            free(matrix);
            printf("distance matrix freed\n");
        }
	}

	dspaces_unlock_on_read("div_lock", &gcomm);

    sprintf(msg, "divergence read lock released ");
    my_message(msg, rank);

    if(table != NULL);
    free_lookup_table(table);

	// DataSpaces: Finalize and clean up DS process
	dspaces_finalize();

	MPI_Barrier(gcomm);
	MPI_Finalize();

	return 0;
}
