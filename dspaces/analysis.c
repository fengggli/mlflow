#include "analysis.h"
#include "ds_adaptor.h"

int main(int argc, char **argv)
{
    // init dspaces
    int nprocs, rank;
    MPI_Comm gcomm;

    char result_path[STRING_LENGTH]="";

    // output results into a folded specified with a slurm jobid
    if(argc == 2){
        strcpy(result_path, argv[1]);
    }


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
     * sampled vel data
     */
    // sampling related
    int sample_size = 40;
    int nprocs_consumer = 32;

    char var_name_sample[STRING_LENGTH];
    sprintf(var_name_sample, "sample");

    int bounds_sampled[6]={0};
    // x_min
    bounds_sampled[0] = 0;

    // x_max
    bounds_sampled[3] = (nprocs_consumer)*(sample_size) -1; 

    int num_elems_sample = (bounds_sampled[3] - bounds_sampled[0])*(bounds_sampled[4] - bounds_sampled[1])*(bounds_sampled[5] - bounds_sampled[2]);
    size_t elem_size_sample = (region_length)*(region_length)*3*sizeof(float);

    float *buffer_regions_sampled = (float *)malloc(num_elems_sample*elem_size_sample);
    if(buffer_regions_sampled== NULL){
        perror("    allocate space for sampled  regions");
        exit(-1);
    }

    /*
     * medoids data
     * all process will receive the same medoids info
     */
    // generate 3 clusters in the end
    int medoids_size = 3;
    char var_name_medoids[STRING_LENGTH];
    sprintf(var_name_medoids, "medoids");

    int bounds_medoids[6]={0};
    // x_min
    bounds_medoids[0] = 0;

    // x_max
    bounds_medoids[3] = medoids_size -1; 

    int num_elems_medoids = medoids_size;
    size_t elem_size_medoids = (region_length)*(region_length)*3*sizeof(float);

    float *buffer_medoids = (float *)malloc(num_elems_medoids*elem_size_medoids);
    if(buffer_medoids== NULL){
        perror("    allocate space for sampled  regions");
        exit(-1);
    }

    /*
     *  further select a subset from all the samples
     */
    int select_size = 48;
    float *buffer_region_selected = (float*)malloc(select_size*elem_size_medoids);

    int num_region = NUM_REGION;
    /*
     * npdiv parameter
     */
    int k_npdiv = K_NPDIV;


    sprintf(msg, "k = %d, max_timestep= %d", K_NPDIV,MAX_VERSION );
    my_message(msg, rank, LOG_CRITICAL);

    sprintf(msg, "ds init ok, result path:%s, address: %p,len: %zu\n",result_path, result_path, strlen(result_path) );
    my_message(msg, rank, LOG_CRITICAL);


    // Timestep notation left in to demonstrate how this can be adjusted
    int timestep=0;

    // timer
    double time_comm_sampled = 0;
    double time_comm_cluster = 0;
    double time_comp =0;

    // we will receive each timestamp
    while(timestep < MAX_VERSION){

        // updated on March 2
        // 1. read samples from all consumer procs
        put_common_buffer(timestep, bounds_sampled, rank, &gcomm, var_name_sampled, &buffer_region_sampled, elem_size_region,  &time_comm_sampled);

        // 2. using subset of them
        select_sampled_buffer(buffer_region_sampled, buffer_region_selected,  select_size);

        // 3.  to decide k medoids
        int nclusters = NCLUSTERS;
        int npass = NPASS;
        int clusterid[num_region];
        double error;
        int ifound;


        sprintf(msg, "start clustering");
        my_message(msg, rank, LOG_WARNING);

        t1 = MPI_Wtime();
        kmedoids(nclusters, num_region, matrix, npass, clusterid, &error, &ifound);

        t2 = MPI_Wtime();
        time_comp = t2-t1;

        sprintf(msg, "finished clustering in %.3lf s  time", t2 -t1);
        my_message(msg, rank, LOG_WARNING);

        sprintf(msg, "error is %.3lf, %d times/ %d passes give the best results\n", error, ifound, npass);
        my_message(msg, rank, LOG_WARNING);
        
        // 4. put the k medoids(k regions) to dspaces
        prepare_medoids(buffer_region_selected, buffer_medoids, clusterid);

        put_common_buffer(timestep, bounds_medoids, rank, &gcomm, var_name_medoids, &buffer_medoids, elem_size_medoids, &time_comm_medoids);
        
        double global_time_comm_cluster;
        double global_time_comm_divs;
        double global_time_comp;

        MPI_Reduce(&time_comm_cluster, &global_time_comm_cluster, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comm_divs, &global_time_comm_divs, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);
        MPI_Reduce(&time_comp, &global_time_comp, 1, MPI_DOUBLE, MPI_SUM, 0, gcomm);

        // Print the result
        if (rank == 0) {
          printf("%d Computation Total %lf avg %lf\n",timestep,  global_time_comp , global_time_comp/ (nprocs));
          printf("%d cluster Total %lf avg %lf\n",timestep,  global_time_comm_cluster , global_time_comm_cluster/ (nprocs));
          printf("%d divs Total %lf avg %lf\n",timestep,  global_time_comm_divs , global_time_comm_divs/ (nprocs));
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
