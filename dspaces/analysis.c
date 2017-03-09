#include "analysis.h"
#include "ds_adaptor.h"

static void fill_div_matrix(double **matrix, float *buffer_divs, int num_region){
    int count =0;
    int i, j;
    for(i = 1; i < num_region; i++){
        for(j = 0; j < i;j++){
            // its the same order as when its saved
            matrix[i][j] = buffer_divs[count];
            // also write to divergence file
            count +=1;
        }
    }
}

// make sure this works
static void prepare_medoids(int *buffer_medoids, int *clusterids, int num_elems, int *cluster_k){
    int i,j,n;
    n = 0;
    for (i = 0; i < num_elems; ++i)
    {
        //printf("processing No %d element %d\n", i, clusterids[i]);
        
        int find = 0;
        for (j = 0; j < n; ++j)
        {
            //printf("\t %d component %d\n", j, clusterids[j]);
            // find new uniq value
            if (buffer_medoids[j] ==clusterids[i]){
                find =1;
               break;
            }
        }

        if (find == 0){
            buffer_medoids[n++] = clusterids[i];
            //printf("\t\t i=%d, j=%d add new uniq value %d\n", i, j, clusterids[i]);
        }
        // at most cluster_k values
    }
    *cluster_k = n;
}


int main(int argc, char **argv)
{
    // init dspaces
    int nprocs, rank;
    MPI_Comm gcomm;

    int nprocs_sim = (PROCS_PER_DIM*PROCS_PER_DIM);


    // MPI communicator
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
    gcomm = MPI_COMM_WORLD;

    // Initalize DataSpaces
    // # of Peers, Application ID, ptr MPI comm, additional parameters
    // # Peers: Number of connecting clients to the DS server
    // Application ID: Unique idenitifier (integer) for application
    // Pointer to the MPI Communicator, allows DS Layer to use MPI barrier func
    // Addt'l parameters: Placeholder for future arguments, currently NULL.
    char msg[STRING_LENGTH];

    printf("trying init dspaces for %d process\n", nprocs);
    dspaces_init(nprocs, 3, &gcomm, NULL);

    /*
     * divs data
     */
    // divs related
    int sample_size = SAMPLE_SIZE;

    int num_elems_sample_all = (nprocs_sim)*(sample_size); 
    int num_region = num_elems_sample_all;
    // divs tasks
    int num_tasks = num_elems_sample_all*(num_elems_sample_all-1)/2;

    char var_name_divs[STRING_LENGTH];
    sprintf(var_name_divs, "divs");

    uint64_t gdims_divs[1] = {num_tasks};
    dspaces_define_gdim(var_name_divs, 1,gdims_divs);

    int bounds_divs[6]={0};

    // x_min
    bounds_divs[0] = 0;

    // x_max
    bounds_divs[3] = num_tasks-1; 

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

    uint64_t gdims_medoids[1] = {10};
    dspaces_define_gdim(var_name_medoids, 1,gdims_medoids);

    int bounds_medoids[6]={0};
    // x_min
    bounds_medoids[0] = 0;

    // x_max
    bounds_medoids[3] = medoids_size -1; 

    int num_elems_medoids = medoids_size;
    size_t elem_size_medoids = sizeof(int);

    int *buffer_medoids = (int *)malloc(num_elems_medoids*elem_size_medoids);
    if(buffer_medoids== NULL){
        perror("    allocate space for sampled  regions");
        exit(-1);
    }

    /*
     * div matrix (ragged)
     */
    double **matrix = malloc(num_region*sizeof(double *));

    if(matrix == NULL){
        perror("malloc for div matrix");
        exit(-1);
    }
    matrix[0] = NULL;
    int i,j;
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

    /*
     *  further select a subset from all the samples
     */
    //int select_size = 48;
    //float *buffer_region_selected = (float*)malloc(select_size*elem_size_medoids);

    /*
     * npdiv parameter
     */

    sprintf(msg, "k = %d, max_timestep= %d", K_NPDIV,MAX_VERSION );
    my_message(msg, rank, LOG_CRITICAL);


    // Timestep notation left in to demonstrate how this can be adjusted
    int timestep=0;

    // timer
    double time_comm_medoids = 0;
    double time_comm_divs = 0;
    double time_comp =0;
    double t1, t2;

    // we will receive each timestamp
    while(timestep < MAX_VERSION){

        printf("********************timestep %d now start!\n",timestep);
        MPI_Barrier(gcomm);
        // updated on March 2
        // 1. read divs from all consumer procs
        get_common_buffer(timestep, 1, bounds_divs, rank, &gcomm, var_name_divs,(void **) &buffer_divs, elem_size_divs,  &time_comm_divs);


        // 1.1. using subset of them
        //select_sampled_buffer(buffer_region_sampled, buffer_region_selected,  select_size);

        // 2. save divs into ragmatrix
        fill_div_matrix(matrix, buffer_divs, num_region);

        // 3.  to decide k medoids
        int nclusters = NCLUSTERS;
        int npass = NPASS;
        int clusterids[num_region];
        double error;
        int ifound;

        printf("\tstart clustering\n");
        my_message(msg, rank, LOG_WARNING);

        t1 = MPI_Wtime();
        kmedoids(nclusters, num_region, matrix, npass, clusterids, &error, &ifound);
        t2 = MPI_Wtime();
        time_comp = t2-t1;

        printf("\tfinished clustering in %.3lf s  time\n", t2 -t1);
        printf("\terror is %.3lf, %d times/ %d passes give the best results\n", error, ifound, npass);
        
        // 4. put the k medoids(k region ids) to dspaces
        int ncluster_2;
        prepare_medoids(buffer_medoids, clusterids, num_region, &ncluster_2);
        printf("number of cluster %d, medoids: %d %d %d \n", ncluster_2, buffer_medoids[0], buffer_medoids[1], buffer_medoids[2]);

        MPI_Barrier(gcomm);
        put_common_buffer(timestep,1,  bounds_medoids, rank, &gcomm, var_name_medoids, (void **)&buffer_medoids, elem_size_medoids, &time_comm_medoids);
        
        double global_time_comm_divs;
        double global_time_comp_medoids;
        double global_time_comm_medoids;

        MPI_Reduce(&time_comm_divs, &global_time_comm_divs, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comp, &global_time_comp_medoids, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comm_medoids, &global_time_comm_medoids, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);

        // Print the result
        if (rank == 0) {
          printf("%d comm divs Total %lf avg %lf\n",timestep,  global_time_comm_divs , global_time_comm_divs/ (nprocs));
          printf("%d comp medoids Total %lf avg %lf\n",timestep,  global_time_comp_medoids, global_time_comp_medoids/ (nprocs));
          printf("%d comm medoids %lf avg %lf\n",timestep,  global_time_comm_medoids , global_time_comm_medoids/ (nprocs));
        }
        timestep++;
    }
    
    sprintf(msg, "now finalize the dspaces and exit");
    my_message(msg, rank, LOG_CRITICAL);

    // DataSpaces: Finalize and clean up DS process
    dspaces_finalize();

    MPI_Barrier(gcomm);
    MPI_Finalize();

    return 0;
}
