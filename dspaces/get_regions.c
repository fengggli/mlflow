#include "get_regions.h"

void get_pair_index(i, range,a, b);

void generate_lookup_table(int num_region, int **p_table){
    int i, j, count;

    int num_pair = num_region*(num_region)/2;
    int *table = (int *)malloc(num_pair*sizeof(int)*2);
    if(tablle == NULL){
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
    // region definition
    // those parameters are obtained after dividing 
    // !! there should be communication between two application
    int region_length = 10;
    int num_region = 400;

	int err;
	int nprocs, rank, i;
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
	sprintf(var_name, "velocity_data");
	dspaces_lock_on_read("velocity_lock", &gcomm);

    // !!!feng: we should know how many regions in total so that we can asign the pairs
    // MPI_receive here?
    int num_region;
    int num_tasks = num_region*(num_region-1)/2;

	// Each process will need to compute its DataSpace index
	int tasks_per_proc = num_tasks/nprocs;
	int tasks_left_over = num_tasks%nprocs;
    int pair_index_l, pair_index_h;
	int ds_lb_index, ds_ub_index;
	
	if(rank<tasks_left_over){
		pair_index_l = rank*(tasks_per_proc+1);
		pair_index_h = pair_index_l + tasks_per_proc;
	}else{
		pair_index_l = rank*tasks_per_proc+tasks_left_over;
		pair_index_h = ds_lb_index + tasks_per_proc-1;
	}

    region_memory_size = (region_length+1)*(region_length+1)*3*sizeof(float);

    // generate pair lookup table;
    int *table;
    generate_lookup_table(num_region, &table);


    // I should wait after all regions are placed into dspaces
    MPI_Reiceive

    for(i = pair_index_l; i<= pair_index_h; i++){
        // the index of the two pairs
        int a, b;
        float div;

        // get the region index of that pair
        get_pair_index(table, i, &a, &b);

        float *buffer a = (float *)malloc(region_memory_size);
        float *buffer b = (float *)malloc(region_memory_size);
        // Define the dimensionality of the data to be received 
        int ndim = 1; 
            
        // Prepare LOWER and UPPER bound dimensions
        uint64_t lb[3] = {0}, ub[3] = {0};
        // get buffer for left region
        lb[0]= a;
        ub[0] =a;
        dspaces_get(var_name, 1, region_memory_size, ndim, lb, ub, buffer_a);

        // get buffer for right region
        lb[0]= b;
        ub[0] =b;
        dspaces_get(var_name, 1, region_memory_size, ndim, lb, ub, buffer_b);


        // we can get the divergence now!
        div = get_divs( buffer_a , buffer_b, region_length, k, div_func);

        free(buffer_a);
        free(buffer_b);

        // put the divergence into divergence matrix
        sprintf(div_name, "div_data");
		dspaces_lock_on_write("div_lock", &gcomm);

        // save divergence in a 1-d array! this will save space
        ndim = 1;
        lb[3] = {0}, ub[3] = {0};

        //!!!! pay attention to the order of dimensions!
        lb[0] = i, ub[0] = i;
		dspaces_put(var_name, timestep, sizeof(float), ndim, lb, ub, div);

        // how about the symmetric part?
		dspaces_unlock_on_write("div_lock", &gcomm);
    }
    // now we can release velocity lock
	dspaces_unlock_on_read("velocity_lock", &gcomm);

    // should wait until get all the divergence
    MPI_Barrior(gcomm);
	// Report data to user
	if(rank==0){
        // get the clustering done here
        sprintf(div_name, "div_data");
		dspaces_lock_on_write("div_lock", &gcomm);


        //reconstruct the divergence matrix, this is a ragged matrix, only le 
        //see the distancematrix function in cluster.c 
        double **matrix = num_region*sizeof(double *); 
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
            for(i = 0 ; i< j; i++)
                free(matrix[i]);
            perror("malloc 2 for div matrix");
            exit(-2);
        }


        ndim = 1;
        lb[3] = {0}, ub[3] = {0};

        int count = 0;
        for(i = 1, i < n, i++){
            for(j = 0; j < i, j++){
                // its the same order as when its saved
                lb[0] = count;
                ub[0] = count;
		        dspaces_get(var_name, timestep, sizeof(float), ndim, lb, ub, matrix[i][j]);
                count +=1;
            }
        }
		dspaces_unlock_on_write("div_lock", &gcomm);

        // clustering parameters
        int nclusters = 3;
        int npass = 100;
        int clusterid[num_region];
        double error;
        int ifound;

        kmdoids(nclusters, num_region, matrix, npass, clusterid, &error, &ifound);

        for(i = 1; i< num_region; i++){
            if(matrix[i] != NULL) free(matrix[i]);
        }
        if(matrix != NULL) {
            free(matrix);
            printf("distance matrix freed\n");
	}

    if(table != NULL);
    free_lookup_table(table);

	// DataSpaces: Finalize and clean up DS process
	dspaces_finalize();

	MPI_Barrier(gcomm);
	MPI_Finalize();

	return 0;
}
