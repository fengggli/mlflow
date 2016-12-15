#include "analysis.h"

void my_message(char *msg, int rank){
    printf("**rank %d: %s\n", rank, msg);
}


int main(int argc, char **argv)
{
	int err;
	int nprocs, rank;
	MPI_Comm gcomm;

    
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
	dspaces_init(1, 3, &gcomm, NULL);

    /*
     * npdiv parameter
     */
    int k_npdiv = K_NPDIV;

    char msg[STRING_LENGTH];

    sprintf(msg, "k = %d, max_timestep= %d", K_NPDIV,MAX_VERSION );
    my_message(msg, rank);

    sprintf(msg, "dataspaces init successfully");
    my_message(msg, rank);

        /*
    sprintf(msg, "dataspaces init complete");
    my_message(msg, rank);
    */

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
    FILE * f_divs;

    // we will receive each timestamp
    while(timestep < MAX_VERSION){
        snprintf(lock_name_regions, STRING_LENGTH, "region_lock_t_%d", timestep);
        snprintf(lock_name_divs, STRING_LENGTH, "div_lock_t_%d", timestep);

        snprintf(divs_path, STRING_LENGTH,"data/parallel/divs_results/all_dist_%d_k_%d_t_%d.txt", POINTS_SIDE,k_npdiv, timestep);

        FILE * f_divs = fopen(divs_path, "w");

        if(f_divs == NULL){
            perror("cannot access the dist_path file");
            exit(-1);
        }

            if(rank == 0){
                sprintf(msg, "\n********************timestep %d now start!\n",timestep);
                my_message(msg, rank);
            }

            int num_region = NUM_REGION;

                        // parameter for second variable: div
            float div;
            uint64_t lb_div[3] = {0}, ub_div[3] = {0};
            // save divergence in a 1-d array! this will save space
            int ndim_div = 3;
            char var_name_div[STRING_LENGTH];
            sprintf(var_name_div, "div_data");
            int num_tasks = num_region*(num_region-1)/2;

            uint64_t gdim_div[3] = {num_tasks,1,1};
            dspaces_define_gdim(var_name_div, 3,gdim_div);

            sprintf(msg, "now div variable has dimention %d",num_tasks);
            my_message(msg, rank);

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
            my_message(msg, rank);

            // get data from dspaces
            // how many pairs

            float* all_divs = (float*)malloc(num_tasks*sizeof(float));

            lb_div[0] = 0; 
            ub_div[0] = num_tasks-1;

            // every rank need at least acquire the lock
            sprintf(msg, "try to acquired div read lock %s", lock_name_divs);
            my_message(msg, rank);

            dspaces_lock_on_read(lock_name_divs, &gcomm);
            sprintf(msg, "acquired div read lock");
            my_message(msg, rank);



            t1 = MPI_Wtime();
            ret_get = dspaces_get(var_name_div, timestep, sizeof(float), ndim_div, lb_div, ub_div, all_divs);
            t2 = MPI_Wtime();

            dspaces_unlock_on_read(lock_name_divs, &gcomm);

            sprintf(msg, "divergence read lock released ");
            my_message(msg, rank);
            
            // reconstruct the matrix
            if(ret_get != 0){
                        error_flag  = 1;
                        exit(-1);
                    }
            sprintf(msg, "all the divergence is read from dataspaces");
            my_message(msg, rank);


            for(i = 1; i < num_region; i++){
                for(j = 0; j < i;j++){
                    // its the same order as when its saved
                    matrix[i][j] = all_divs[count];
                    // also write to divergence file
                    fprintf(f_divs, "%d\t %d\t %f\n", i, j, all_divs[count]);
                    count +=1;
                }
            }

            t3 = MPI_Wtime();


            sprintf(msg, "divergence matrix read %.3f s,filled in %.3f s time", t2-t1, t3- t2);
            my_message(msg, rank);

            if(error_flag == 1){
                sprintf(msg, "ERROR when read divergence from Dspaces");
                my_message(msg, rank);
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
                my_message(msg, rank);

                t1 = MPI_Wtime();
                kmedoids(nclusters, num_region, matrix, npass, clusterid, &error, &ifound);

                t2 = MPI_Wtime();

                sprintf(msg, "finished clustering in %.3lf s  time", t2 -t1);
                my_message(msg, rank);

                // save cluster results into file
                snprintf(output_path, STRING_LENGTH,"data/parallel/clustering_results/clusterid_%d_k_%d_t_%d.txt", POINTS_SIDE,k_npdiv, timestep);
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

        if(f_divs != NULL)
            fclose(f_divs);
            

        timestep++;
	}

    sprintf(msg, "now finalize the dspaces and exit");
    my_message(msg, rank);

	// DataSpaces: Finalize and clean up DS process
	dspaces_finalize();

	MPI_Barrier(gcomm);
	MPI_Finalize();

	return 0;
}
