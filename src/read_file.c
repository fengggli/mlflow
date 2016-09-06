/*
 * * Get datatype and dataspace identifiers and then query
 * * read only one time slice
 * * dataset class, order, size, rank and dimensions.
 * */

#include "hdf5.h"
#include <stdlib.h>
//#include "read_file.h"
#define FILE        "../data/test_1_2_3_4.h5"
#define DATASET_NAME
#define DATASET_P "p00000" // dataset name for pressure
#define DATASET_U "u00000" // dataset name for velocity 

/*
 * dimensions of dataset
 */
#define NX 4
#define NY 3
#define NZ 2

#define RANK         4 // four dims: x y z and time 
#define p_malloc_error(I) {printf("malloc failuren at %d",I);exit(-1);} 



int main(){
    hid_t       file;                        /* handles */
    hid_t       dataset_p;  
    hid_t       dataset_u;  
    hid_t       filespace_u;                   
    hid_t       memspace_u;                  
    hid_t       cparms;                   
    hsize_t     dims[4];                     /* dataset and chunk dimensions*/ 
    //hsize_t     chunk_dims[2];
    //hsize_t     col_dims[1];

    herr_t      status, status_n;                             

    float  data_p[4][3][2];  /* buffer for pressure to be read */
    float  data_u[4][3][2][3];  /* buffer for velocity to be read */
    //int         chunk_out[2][5];   /* buffer for chunk to be read */
    //int         column[10];        /* buffer for column to be read */
    int         rank;
    hsize_t i, j, k, m;
    

 
    /*
     * Open the file and the dataset.
     */
    file = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset_p = H5Dopen(file, DATASET_P, H5P_DEFAULT);
    dataset_u = H5Dopen(file, DATASET_U, H5P_DEFAULT);
 
    /*
     * Get dataset rank and dimension.
     */

    filespace_u = H5Dget_space(dataset_u);    /* Get filespace handle first. */
    rank      = H5Sget_simple_extent_ndims(filespace_u);
    printf("dataset rank %d \n", rank);
    status_n  = H5Sget_simple_extent_dims(filespace_u, dims, NULL);
    printf("dataset rank %d, dimensions %lu x %lu x %lu x %lu \n",
       rank, (unsigned long)(dims[0]),(unsigned long)(dims[1]),(unsigned long)(dims[2]), (unsigned long)(dims[3]));

    /*
     * Get creation properties list.
     */
    cparms = H5Dget_create_plist(dataset_u); /* Get properties handle first. */

    /* 
     * Check if dataset is chunked.
     */

    if (H5D_CHUNKED == H5Pget_layout(cparms))  {
        printf("its thunked\n");

    /*
     * Get chunking information: rank and dimensions
     */
    /*
    rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
    printf("chunk rank %d, dimensions %lu x %lu\n", rank_chunk,
           (unsigned long)(chunk_dims[0]), (unsigned long)(chunk_dims[1]));
    */
    }

    /*
     * Define the memory space to read dataset.
     */
    memspace_u = H5Screate_simple(RANK,dims,NULL);

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
     * Read dataset back and display.
     */
    
    status = H5Dread(dataset_u, H5T_NATIVE_FLOAT, memspace_u, filespace_u,
             H5P_DEFAULT, data_u);
    printf("\n");
    printf("Dataset: \n");
    

    for (i = 0; i < dims[0]; i++) {
        for (j = 0; j < dims[1]; j++){ 
            for (k = 0; k < dims[2]; k++){
                printf("(");
                for(m=0;m<dims[3];m++)
                    printf("%.6f ", data_u[i][j][k][m]);
                printf(")\t");

            }
        printf("\n");
        }
        printf("\n");
    }
    
    

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


	

