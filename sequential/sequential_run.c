#include "sequential_run.h"
// read a block cutout from hdf file
// divide into regions
// for each region get the density list
// for each pair of regions, get the divergence using their density lists
// output all the divergence into a txt file

/*

extern int read_data(const char* file_name, float **pressure, float **velocity, int *dim1, int *dim2, int *dim3);
extern int free_data(float *pressure, float *velocity);
extern float get_divs(float *A, float *B, int region_length, int k, int div_func);
*/

int get_all_divs(float *regions,int num_region, int d2, int d3, int region_length,double **matrix,char *divs_path, int k_npdiv, int div_func){
    int i, j;

    float  div;

    FILE * f_divs = fopen(divs_path, "w");
    
    for(i = 1; i< num_region; i++){

        //fprintf(stderr, "\t processing divs between region %d and others\n", i);
        for(j = 0; j < i ; j++){

            printf("***processing divs between region %d and %d\n", i, j);

            // starting address of each region as input
            div = get_divs( regions + i*d2*d3 , regions + j*d2*d3, region_length, k_npdiv, div_func);

#ifdef USE_SYNTHETIC
            div = div;
#endif

            printf("***divergence between region %d and %d is %.3f\n", i, j, div);
           // printf("\t %.3lf s time is used for div\n", time_current - time_old);
           
            //fprintf(f_dist,"%f\n",div);


            matrix[i][j] = div;
            fprintf(f_divs, "%d\t %d\t %f\n", i, j, div);
        }
    }

        //fclose(f_dist);
    fclose(f_divs);
    return 0;
}


int run_clustering (int nclusters, int num_region, double **matrix,
  int npass, char *output_path, double* error, int* ifound){
    int i, j;
    double div;
    char str[20];
    int clusterid[num_region];

    printf("***starting to run kmedoids algorithm\n");
    printf("nclusters = %d\nnum of elements %d\n npass %d\n. error = %.3lf\n", nclusters,num_region,npass,*error);
    
    kmedoids (nclusters, num_region, matrix, npass,clusterid,error,ifound);

    // save cluster results into file
    FILE * f_clusterid = fopen(output_path, "w");
    if(f_clusterid == NULL){
        perror("file open error");
        exit(-1);
    }

    for(i = 0; i < num_region; i++){
        fprintf(f_clusterid, "%d\n",clusterid[i]);
    }
    
    fclose(f_clusterid);
    return 0;
}


int main(){
    
    printf(" start.\n");

    /**********************
     * user definition starts here
     */

    int timestep = 0;

    // find k=5 nearest neighbours to determine the density
    int k_npdiv = K_NPDIV;

    // input file
    char file_name[STRING_LENGTH];

    for(timestep = 0; timestep< MAX_VERSION; timestep++  ){

        fprintf(stderr, "timestep %d started \n", timestep);


#ifdef USE_SYNTHETIC
       // divs output file
        char divs_path[STRING_LENGTH];
        sprintf(divs_path, "data/sequential/all_divs_%d_k_%d_t_%d_syn.txt", POINTS_SIDE, k_npdiv, timestep);

        // clustering output file
        char output_path[STRING_LENGTH]; 
        sprintf(output_path, "data/sequential/clusterid_%d_k_%d_t_%d_syn.txt", POINTS_SIDE, k_npdiv, timestep);

#else
        sprintf(file_name,"data/isotropic_%d_%d_1_t_%d.h5",POINTS_SIDE ,POINTS_SIDE, timestep);

        // divs output file
        char divs_path[STRING_LENGTH];
        sprintf(divs_path, "data/sequential/all_divs_%d_k_%d_t_%d_ragged.txt", POINTS_SIDE, k_npdiv, timestep);

        // clustering output file
        char output_path[STRING_LENGTH]; 
        sprintf(output_path, "data/sequential/clusterid_%d_k_%d_t_%d_ragged.txt", POINTS_SIDE, k_npdiv, timestep);
#endif

        /*
         * user definition ends here
         ***************************/

        // timer
        double time_start,time_old,time_current; /* * Part I obtain all the divs */ // set 1: 201*201*1 data cutout, region_length 10 // set 2: 601*601*1 data cutout, region_length 30 //const char *file_name = "data/test_1_2_3_4.h5"; // 
        // divide configrations, length of 10 region will have 11*11 121 points
        int region_length;

        
        float *pressure, *velocity;
        int dim1, dim2, dim3;

        int i,j, p, q, ret;

        region_length = REGION_LENGTH;


        // read data to buffer

#ifdef USE_SYNTHETIC
        dim1 = 41;
        dim2 = 41;
        dim3 = 1;

        printf("synthetic data dimension is %d * %d * %d\n", dim1, dim2, dim3);
#else
        read_data(file_name, &pressure, &velocity, &dim1, &dim2, &dim3);
        //print_p_data(pressure, dim1, dim2, dim3);
        //save visualized result
        
        printf("data is read, dimension is %d * %d * %d\n", dim1, dim2, dim3);
        fprintf(stderr, "data is read, dimension is %d * %d * %d\n", dim1, dim2, dim3);

#endif
        
        // split the slice into regions
        // since we know each region's exact size, we only need record the starting address
        int num_region;

        // num_regions* points_in_region*3
        float  *regions;

        // dim1 = dim2 = 201, in a region, 11 points in each side, thus length 10, the whole block has length 201-1 = 200. So (200/10)^2 400 regions
        // d1 d2 d3 is the 3-dimension matrix 
        int d1 = dim1;
        // fixed on Dec 13
        //int d2 = region_length*region_length;;
        int d2 = (region_length + 1)*(region_length+1);
        int d3 = 3;



#ifdef USE_SYNTHETIC
        divide_synthetic(d1, region_length, &num_region, &regions);
#else
        divide(velocity, d1, region_length, &num_region, &regions);
#endif

        if(num_region != NUM_REGION){
            printf("you should define the NUM_REGION in conf file");
            exit(-1);
        }

        printf("generate %d regions\n", num_region);

        
        // caculate the divergence between all regions
        // distance matrix can be very large
        
        // use L-2 divergence
        int div_func = 1;

        // allocate space the distance matrix
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


        // timer for divergence calculation
        time_start = get_cur_time();
        //time_current = time_start;


        ret =  get_all_divs(regions, num_region,  d2, d3, region_length, matrix, divs_path, k_npdiv, div_func);
        if(ret != 0){
            perror("divergence calculation error");
            exit(-1);
        }else{
            // give time
            
            time_current = get_cur_time();
            printf("distance matrix is saved in matrix");
            printf("\t %.3lf s time is used for all divsdiv\n", time_current - time_start);
            fprintf(stderr, "\t %.3lf s time is used for all divsdiv\n", time_current - time_start);
        }

        //FILE * f_dist = fopen(dist_path, "w");
        // if regions are just pointers, no need to free 
        if(regions != NULL)
            free(regions);
        // free buffer
#ifndef USE_SYNTHETIC
        if(pressure != NULL|| velocity != NULL)
            free_data(pressure, velocity);
#endif
        // also free pdata block


        /*
         * Part II kmedoids clustering
         */


        int nclusters ,npass, ifound;
        double error;

        // timer

        nclusters = NCLUSTERS;
        num_region =NUM_REGION ; 
        npass =NPASS;

        //int clusterid[num_region];


        time_start = get_cur_time();

        //call the divergence function
        printf("starting clustering , num_regions is %d\n",  num_region );
        run_clustering (nclusters,num_region, matrix,npass,output_path, &error,&ifound);

        time_current = get_cur_time();
        printf("clustering completed, result  is stored in %s\n", output_path);
        fprintf(stderr, "clustering completed, result  is stored in %s\n", output_path);
        printf("%.4f time is used for clustering\n", time_current - time_start );
        printf("error is %.3lf, %d times/ %d passes give the best results\n", error, ifound, npass);


        // free the divergence buffer
        for(i = 1; i< num_region; i++){
            if(matrix[i] != NULL) free(matrix[i]);
        }
        if(matrix != NULL) {
            free(matrix);
            printf("distance matrix freed\n");
        }

        // end of the timestep
    }

   	    
    return 1;
}
