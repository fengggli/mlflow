#include "run_with_dspaces.h"
#include "divide.h"
#include "read_file.h"

int generate_regions(char *hdfpath, int region_length, int *p_num_region, float **p_regions){
    
    printf("    read data to buffer\n");

    int dim1, dim2, dim3;
    
    float *pressure, *velocity;

    // read data into buffer
    read_data(hdfpath, &pressure, &velocity, &dim1, &dim2, &dim3);

    // dim1 is the dimension of x, we assume that the datacut is a square: dim1 = dim2
    divide(velocity, dim1, region_length, p_num_region, p_regions);

    printf("    all the regions are divided\n");
    
    // free buffer
    free_data(pressure, velocity);
    printf("    tmp buffer freeed, now prepare to flush into dataspaces\n ");
    return 1;
}

double get_cur_time() {
  struct timeval   tv;
  struct timezone  tz;
  double cur_time;

  gettimeofday(&tv, &tz);
  cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;

  return cur_time;
}
