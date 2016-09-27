#include "read_file.h"

/*
 * read data from hdf5 file
 * input:
 *      file_name: path of the hdf5 file
 * output:
 *      pressure: pressure of this block 
 *      velocity:
 *      dim1, dim2, dim3: size of this block
 *
 */

int read_data(const char* file_name, float **pressure, float **velocity, int *dim1, int *dim2, int *dim3){
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
    file = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
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

    //float  data_p[dims_p[0]][dims_p[1]][dims_p[2]];  /* buffer for pressure to be read */
    float  *data_p = (float *)malloc(sizeof(float)*dims_p[0]*dims_p[1]*dims_p[2]);
    *pressure = data_p;

    
    /*
     * Get dataset rank and dimension for velocity data.
     */

    filespace_u = H5Dget_space(dataset_u);    /* Get filespace handle first. */
    rank_u      = H5Sget_simple_extent_ndims(filespace_u);
    //printf("dataset rank %d \n", rank_u);
    status_n  = H5Sget_simple_extent_dims(filespace_u, dims_u, NULL);
    printf(" velocity dataset rank %d, dimensions %lu x %lu x %lu x %lu \n",
       rank_u, (unsigned long)(dims_u[0]),(unsigned long)(dims_u[1]),(unsigned long)(dims_u[2]), (unsigned long)(dims_u[3]));

    *dim1 = dims_u[0];
    *dim2 = dims_u[1];
    *dim3 = dims_u[2];

    // float  data_u[dims_u[0]][dims_u[1]][dims_u[2]][3];  /* buffer for velocity to be read */
    float * data_u = (float *)malloc(sizeof(float)*dims_u[0]*dims_u[1]*dims_u[2]*3);  /* buffer for velocity to be read */
    *velocity =data_u;

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

    return 0;
}
