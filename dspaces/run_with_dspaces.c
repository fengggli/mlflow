#include "run_with_dspaces.h"
#include "divide.h"
#include "read_file.h"

int generate_regions(char *hdfpath, int region_length, int *p_num_region, float **p_regions){
    
    printf("    read data from %s to buffer\n", hdfpath);

    int dim1, dim2, dim3;
    
    float *pressure, *velocity;

    // read data into buffer
    // sythetic data to verify
#ifdef USE_SYTHETIC
    dim1 = POINTS_SIDE; // x and y all have dimension size of 201

    // generate fake data
    divide_synthetic(&velocity, dim1, regions_length, p_num_regions, p_regions);

#else
    read_data(hdfpath, &pressure, &velocity, &dim1, &dim2, &dim3);

    // dim1 is the dimension of x, we assume that the datacut is a square: dim1 = dim2
    divide(velocity, dim1, region_length, p_num_region, p_regions);

    printf("    all the regions are divided\n");
    
    // free buffer
    free_data(pressure, velocity);
    printf("    tmp buffer freeed, now prepare to flush into dataspaces\n ");
#endif

    
    return 1;
}


