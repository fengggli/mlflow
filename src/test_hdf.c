#include "read_file.h"

int main(){
    printf(" read data to buffer\n");
    //const char *file_name = "data/test_1_2_3_4.h5";
    const char *file_name = "data/isotropic_255_255_5.h5";
    float *pressure, *velocity;
    int dim1, dim2, dim3;

    // read data to buffer
    read_data(file_name, &pressure, &velocity, &dim1, &dim2, &dim3);
    //print_p_data(pressure, dim1, dim2, dim3);
    //save visualized result
    
    printf("data is read, dimension is %d * %d * %d\n", dim1, dim2, dim3);
    
    // free buffer
    free_data(pressure, velocity);



    return 1;
}
