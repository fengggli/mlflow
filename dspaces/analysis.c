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

    int num_region = NUM_REGION;

#ifdef INCLUDE_ML
    // start of ds_adaptor preparation

    // data layout
    //int dims[3] = {1, POINTS_SIDE, POINTS_SIDE};
    //int num_points = dims[0]*dims[1]*dims[2];



    char var_name_cluster[STRING_LENGTH];
    sprintf(var_name_cluster, "CLUSTER");

    // prepare space, 
    float * cluster_data = (float *)malloc(num_region* sizeof(float));
    // end of dspaces preparation
#endif

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

    char lock_name_regions[STRING_LENGTH];
    char lock_name_divs[STRING_LENGTH];


    char output_path[STRING_LENGTH];
    /*
       snprintf(lock_name_regions, STRING_LENGTH, "region_lock_same");
       snprintf(lock_name_divs, STRING_LENGTH, "div_lock_same");
       */


    // save all the divergence
    char divs_path[STRING_LENGTH];

    int side_num_region = (POINTS_SIDE -1)/REGION_LENGTH; // this will be (201-1)/10 = 20

    // timer
    double time_comm_divs = 0;
    double time_comm_cluster = 0;
    double time_comp =0;

    // we will receive each timestamp
    while(timestep < MAX_VERSION){
        snprintf(lock_name_regions, STRING_LENGTH, "region_lock_t_%d", timestep);
        snprintf(lock_name_divs, STRING_LENGTH, "div_lock_t_%d", timestep);


        if(argc == 2){
            snprintf(divs_path, STRING_LENGTH,"%s/all_divs/%d_k_%d_t_%d.txt",result_path,  POINTS_SIDE,k_npdiv, timestep);
        }
        else{
            snprintf(divs_path, STRING_LENGTH,"data/parallel/divs_results/all_dist_%d_k_%d_t_%d.txt", POINTS_SIDE,k_npdiv, timestep);

        }
        FILE * f_divs = fopen(divs_path, "w");

        if(f_divs == NULL){
            printf("file %s not found, exit", divs_path);
            exit(-1);
        }

        if(rank == 0){
            sprintf(msg, "\n********************timestep %d now start!\n",timestep);
            my_message(msg, rank, LOG_WARNING);
        }


        // parameter for second variable: div
        uint64_t lb_div[3] = {0}, ub_div[3] = {0};
        // save divergence in a 1-d array! this will save space
        int ndim_div = 3;
        char var_name_div[STRING_LENGTH];
        sprintf(var_name_div, "div_data");
        int num_tasks = num_region*(num_region-1)/2;

        uint64_t gdim_div[3] = {num_tasks,1,1};
        dspaces_define_gdim(var_name_div, 3,gdim_div);

        sprintf(msg, "now div variable has dimention %d",num_tasks);
        my_message(msg, rank, LOG_WARNING);

        double t1, t2, t3;
        int i, j, ret_get;


        // now anlysis
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
        my_message(msg, rank, LOG_WARNING);

        // get data from dspaces
        // how many pairs

        float* all_divs = (float*)malloc(num_tasks*sizeof(float));

        lb_div[0] = 0; 
        ub_div[0] = num_tasks-1;

        // every rank need at least acquire the lock
        sprintf(msg, "try to acquired div read lock %s", lock_name_divs);
        my_message(msg, rank, LOG_WARNING);

        dspaces_lock_on_read(lock_name_divs, &gcomm);
        sprintf(msg, "acquired div read lock");
        my_message(msg, rank, LOG_WARNING);



        t1 = MPI_Wtime();
        ret_get = dspaces_get(var_name_div, timestep, sizeof(float), ndim_div, lb_div, ub_div, all_divs);
        t2 = MPI_Wtime();

        time_comm_divs=t2-t1;

        dspaces_unlock_on_read(lock_name_divs, &gcomm);

        sprintf(msg, "divergence read lock released ");
        my_message(msg, rank, LOG_WARNING);

        // reconstruct the matrix
        if(ret_get != 0){
            error_flag  = 1;
            exit(-1);
        }
        sprintf(msg, "all the divergence is read from dataspaces");
        my_message(msg, rank, LOG_WARNING);


        for(i = 1; i < num_region; i++){
            for(j = 0; j < i;j++){
                // its the same order as when its saved
                matrix[i][j] = all_divs[count];
                // also write to divergence file
                fprintf(f_divs, "%d\t %d\t %f\n", i, j, all_divs[count]);
                count +=1;
            }
        }

        if(f_divs != NULL)
            fclose(f_divs);

        t3 = MPI_Wtime();


        sprintf(msg, "divergence matrix read %.3f s,filled in %.3f s time,also saved in %s", t2-t1, t3- t2, divs_path);
        my_message(msg, rank, LOG_WARNING);

        if(error_flag == 1){
            sprintf(msg, "ERROR when read divergence from Dspaces");
            my_message(msg, rank, LOG_CRITICAL);
        }
        // if everything is fine we start clustering
        else{

            // clustering parameters
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

            // save cluster results into file
            if(argc == 2){
                strcpy(result_path, argv[1]);
            }
            sprintf(msg, "tring to write clusterid to parent path:%s,address: %p,len: %zu\n",result_path, result_path, strlen(result_path) );
            my_message(msg, rank, LOG_CRITICAL);

            if(argc == 2){
                snprintf(output_path, STRING_LENGTH,"%s/clusterids/%d_k_%d_t_%d.txt",result_path,  POINTS_SIDE,k_npdiv, timestep);
            }
            else{
                snprintf(output_path, STRING_LENGTH,"data/parallel/clustering_results/clusterid_%d_k_%d_t_%d.txt", POINTS_SIDE,k_npdiv, timestep);
            }
            // write cluster id to dspaces
            // cast to float


#ifdef INCLUDE_ML
            for(i = 0; i < num_region;i++ ){
                cluster_data[i] = clusterid[i];
            }
            put_cluster_buffer(timestep, &num_region ,rank, &gcomm, var_name_cluster,  &cluster_data, &time_comm_cluster);
#endif
            // dspaces put

            /* no need to write here
            // write into filesystem
            FILE * f_clusterid = fopen(output_path, "w");
            if(f_clusterid == NULL){
                printf("file %s not found, exit\n", output_path);
                exit(-1);
            }
            else{
                sprintf(msg, "clustering results saved  in %s", output_path);
                my_message(msg, rank, LOG_CRITICAL);
            }

            //print the header
            fprintf(f_clusterid, "x y z clusterid\n");

            int count1 = 0;
            for(i = 0; i < side_num_region;i++ ){
                for(j = 0; j< side_num_region; j++){
                    fprintf(f_clusterid, "%d %d %d %d\n",i, j, 0, clusterid[count1++]);
                }
            }

            fclose(f_clusterid);
            */
        }

        // free the divergence buffer
        for(i = 1; i< num_region; i++){
            if(matrix[i] != NULL) free(matrix[i]);
        }
        if(matrix != NULL) {
            free(matrix);
            printf("distance matrix freed\n");
        }


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
