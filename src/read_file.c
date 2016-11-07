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
    hsize_t     dims_p[4];                     /* dataset and chunk dimensions*/ 
    hsize_t     dims_u[4];                     /* dataset and chunk dimensions*/ 

    /*
     * Open the file and the dataset.
     */
    //float  data_p[dims_p[0]][dims_p[1]][dims_p[2]];  /* buffer for pressure to be read */
    float  *data_p = (float *)malloc(sizeof(float)*dims_p[0]*dims_p[1]*dims_p[2]);

    dims_u[0]= 201;
    dims_u[1]= 201;
    dims_u[2]= 1;
    dims_u[3]= 3;

    //printf("dataset rank %d \n", rank_u);
    printf(" velocity dataset, dimensions %lu x %lu x %lu x %lu \n",
        (unsigned long)(dims_u[0]),(unsigned long)(dims_u[1]),(unsigned long)(dims_u[2]), (unsigned long)(dims_u[3]));

    *dim1 = dims_u[0];
    *dim2 = dims_u[1];
    *dim3 = dims_u[2];

    // float  data_u[dims_u[0]][dims_u[1]][dims_u[2]][3];  /* buffer for velocity to be read */
    float * data_u = (float *)malloc(sizeof(float)*dims_u[0]*dims_u[1]*dims_u[2]*3);  /* buffer for velocity to be read */
    *velocity =data_u;

        return 0;
}

int free_data(float *pressure, float *velocity){
    printf("    try to free pressure and velocity\n");
    if(pressure != NULL){
        free(pressure);
    }else{
        printf("    pressure free error\n");
    }
    if(velocity != NULL){
        free(velocity);
    }else
        printf("    velocity free error\n");
    printf("buffer freed\n");
    return 1;
}
