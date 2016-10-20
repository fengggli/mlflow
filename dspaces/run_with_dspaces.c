#include "run.h"
#include "run_with_dspaces.h"
#include "divide.h"

int generate_regions(char *hdfpath, int region_length, int *p_num_region, float **p_regions){
    
    printf(" read data to buffer\n");

    int dim1, dim2, dim3;
    
    float *pressure, *velocity;

    // read data into buffer
    read_data(file_name, &pressure, &velocity, &dim1, &dim2, &dim3);

    int num_regions;
    float *regions;

    // dim1 is the dimension of x, we assume that the datacut is a square: dim1 = dim2
    divide(velocity, dim1, region_length, &num_region, &regions);

    int region_length;

    data
    float *velocity;
    int dim1, dim2, dim3;

    int num_region;

    // feed the the address of regions here
    float  *regions = ;

    // dim1 = dim2 = 201, in a region, 11 points in each side, thus length 10, the whole block has length 201-1 = 200. So (200/10)^2 400 regions
    // d1 d2 d3 is the 3-dimension matrix 
    int d1 = dim1;
    int d2 = region_length*region_length;;
    int d3 = 3;
    divide(velocity, d1, region_length, &num_region, &regions);
    *p_regions = regions;

    
    // free buffer
    free_data(pressure, velocity);
    return 1;
}
