/*
 * * Get datatype and dataspace identifiers and then query
 * * read only one time slice
 * * dataset class, order, size, rank and dimensions.
 * */

#include "hdf5.h"
#include <stdlib.h>
//#include "read_file.h"
// #define FILE        "../data/test_1_2_3_4.h5"
#define FILE        "../data/isotropic_255_255_5.h5"
#define DATASET_NAME
#define DATASET_P "p00000" // dataset name for pressure
#define DATASET_U "u00000" // dataset name for velocity 

#define RANK         4 // four dims: x y z and time 
#define p_malloc_error(I) {printf("malloc failuren at %d",I);exit(-1);} 


/*
int get_data(const char *file_name, const char *dataset_name, float ****pressure, float *****velocity){}
int free_data(float ****pressure, float *****velocity){}
*/


int main(){
    hid_t       file;                        /* handles */
    hid_t       dataset_p;  
    hid_t       dataset_u;  
    hid_t       filespace_p;                   
    hid_t       filespace_u;                   
    hid_t       memspace_p;                  
    hid_t       memspace_u;                  
    hid_t       cparms_p;                   
    hid_t       cparms_u;                   
    hsize_t     dims_p[4];                     /* dataset and chunk dimensions*/ 
    hsize_t     dims_u[4];                     /* dataset and chunk dimensions*/ 
    //hsize_t     chunk_dims[2];
    //hsize_t     col_dims[1];

    //hsize_t     nx=4, ny=3,nz=2;

    herr_t      status, status_n;                             

    //int         chunk_out[2][5];   /* buffer for chunk to be read */
    //int         column[10];        /* buffer for column to be read */
    int         rank_p;
    int         rank_u;
    hsize_t i, j, k, m;
    

 
    /*
     * Open the file and the dataset.
     */
    file = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset_p = H5Dopen(file, DATASET_P, H5P_DEFAULT);
    dataset_u = H5Dopen(file, DATASET_U, H5P_DEFAULT);
 
    /*
     * Get dataset rank and dimension for pressure data.
     */
    filespace_p = H5Dget_space(dataset_p);    /* Get filespace handle first. */
    rank_p      = H5Sget_simple_extent_ndims(filespace_p);
    //printf("dataset rank %d \n", rank_p);
    status_n  = H5Sget_simple_extent_dims(filespace_p, dims_p, NULL);
    printf("pressure dataset rank %d, dimensions %lu x %lu x %lu x %lu \n",
       rank_p, (unsigned long)(dims_p[0]),(unsigned long)(dims_p[1]),(unsigned long)(dims_p[2]), (unsigned long)(dims_p[3]));

    float  data_p[dims_p[0]][dims_p[1]][dims_p[2]];  /* buffer for pressure to be read */

    
    /*
     * Get dataset rank and dimension for velocity data.
     */

    filespace_u = H5Dget_space(dataset_u);    /* Get filespace handle first. */
    rank_u      = H5Sget_simple_extent_ndims(filespace_u);
    //printf("dataset rank %d \n", rank_u);
    status_n  = H5Sget_simple_extent_dims(filespace_u, dims_u, NULL);
    printf(" velocity dataset rank %d, dimensions %lu x %lu x %lu x %lu \n",
       rank_u, (unsigned long)(dims_u[0]),(unsigned long)(dims_u[1]),(unsigned long)(dims_u[2]), (unsigned long)(dims_u[3]));

    float  data_u[dims_u[0]][dims_u[1]][dims_u[2]][3];  /* buffer for velocity to be read */

    /*
     * Get creation properties list.
     */
    cparms_u = H5Dget_create_plist(dataset_u); /* Get properties handle first. */
    cparms_p = H5Dget_create_plist(dataset_p); /* Get properties handle first. */

    /* 
     * Check if dataset is chunked.
     */

    if (H5D_CHUNKED == H5Pget_layout(cparms_p))  {
        printf("pressure data is  thunked\n");

    /*
     * Get chunking information: rank and dimensions
     */
    /*
    rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
    printf("chunk rank %d, dimensions %lu x %lu\n", rank_chunk,
           (unsigned long)(chunk_dims[0]), (unsigned long)(chunk_dims[1]));
    */
    }
    if (H5D_CHUNKED == H5Pget_layout(cparms_u))  {
        printf("velocity data is  thunked\n");
    }

        /*
     * allocate space for dataout
     * do we need dynamic memory?
     */
    /*
    data_out = (float ***)malloc(dims[0]*sizeof(float **));
    if(data_out!=NULL){
        for(i = 0; i< dims[0]; i++){
            data_out[i] = (float **)malloc(dims[1]*sizeof(float *));
                if(data_out[i] !=NULL){
                    for(j = 0; j < dims[1];j++){
                        data_out[i][j] = (float*)malloc(dims[2]*sizeof(float));
                        printf("allocate %d float for data[%d][%d]\n", dims[2],i,j);
                        if(data_out[i][j]==NULL)
                            p_malloc_error(1);
                    }
                }
                else p_malloc_error(2);
        }
    }
    else
        p_malloc_error(3);
        
    */
    
 
     /*
     * Define the memory space to read dataset.
     */
    memspace_p = H5Screate_simple(RANK,dims_p,NULL);
    memspace_u = H5Screate_simple(RANK,dims_u,NULL);


     /*
     * Read dataset back and display.
     */

    
    status = H5Dread(dataset_p, H5T_NATIVE_FLOAT, memspace_p, filespace_p,
             H5P_DEFAULT, data_p);
    status = H5Dread(dataset_u, H5T_NATIVE_FLOAT, memspace_u, filespace_u,
             H5P_DEFAULT, data_u);
    printf("\n");
    printf(" pressure Dataset: \n");

    /*
     *
    split_region(data_p, data_u, block_size, overlap =0);
    get_divs(xbags, num_bags, div_func, results);
    kmedoids (int nclusters, int nelements, double** distance,
        int npass, int clusterid[], double* error, int* ifound);
    */

    
    H5Dclose(dataset_p);
    H5Sclose(filespace_p);
    H5Sclose(memspace_p);

    H5Dclose(dataset_u);
    H5Sclose(filespace_u);
    H5Sclose(memspace_u);

    H5Fclose(file);
    /*
    
    for(i =0; i<  dims[0]; i++){
        for(j = 0; j < dims[1]; j++){
            free(data_out[i][j]);
        }
        free(data_out[i]);
    }
    free(data_out);
    */
    printf("all buffer is freed\n");

    return 0;
}


	

