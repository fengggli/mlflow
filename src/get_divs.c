#include "get_divs.h"
//#define debug_1

// parallel version and sequential version have different divergence
#define DEBUG_REGION_MAPPING



float  get_k_order(float *tmp_array, int length, int k){
    int i, j;
    float min;
    int min_index = -1;
    for(i=0;i<k-1;i++){
        min = LARGE_NUMBER;
        for(j = 0; j< length; j++){
            if(tmp_array[j] < min){
                min = tmp_array[j];
                min_index = j;
            }
        }

        tmp_array[min_index] = LARGE_NUMBER;
    }
    min = LARGE_NUMBER;
    for(j = 0; j< length; j++){
            if(tmp_array[j] < min){                                                                                                                                                                                                           
                min = tmp_array[j];
            }
    }
    return min;
    // last round get smallest value
}

void print_matrix(float *Dismatrix, int l, int x0, int y0, int x1, int y1){
    printf("print matrix\n");
    int i,j;
    for(i = x0; i < x1; i++){
        for(j = y0; j < y1; j++){
            printf("%.6f\t", Dismatrix[i*l + j]);
        }
        printf("\n");
    }
}

/* 
 * inside the region, get the distance to a point's k-neares neighbours
 * input:
 *  i, the index of current point
 *  DisMatrix, the distance lookup table
 *  length: length of the distance matrix
 *  k, how many nearest neighbours
 * return:
 *  the distance to k-nearest neighbours
 */

// there should be a bug here, if the point comes from current region, what should i do? consistency with the (n-1) in the density estimation
float get_bound_dist(int i, float *A, int length, int k){
    // actually this is a task to sort all the elements in the i-th row in DisMatrix

    // the index of k-th least distance
    int j;
    // values of this point
    float tmp_x = *(A + 3*i);
    float tmp_y = *(A + 3*i +1);
    float tmp_z = *(A + 3*i +2);

    // get distances with all other points
    float *tmp_array = (float *)malloc(sizeof(float)*length);
    if(tmp_array == NULL){
        perror("malloc errr");
        exit(-1);
    }

    float dist_tmp;

    for(j = 0; j< length; j++){
        if(j == i){
            dist_tmp = LARGE_NUMBER;
        }
        else{
           dist_tmp = pow(tmp_x - *(A + 3*j),2) + pow(tmp_y - *(A + 3*j +1), 2)+ pow(tmp_z - *(A+3*j+2),2);
        }
        if(dist_tmp == 0){
            printf("A[%d]:(%f, %f, %f), B[%d]:(%f, %f, %f)\n", i, tmp_x, tmp_y, tmp_z, j, *(A + 3*j),*(A + 3*j+1),*(A + 3*j+2));
        }
        tmp_array[j] = dist_tmp;
        //printf("\tset %d th distance\n", j);
    }
    //float *array_to_sort = DisMatrix + i*length;
    float ret =  get_k_order(tmp_array, length, k);
    free(tmp_array);
    return ret;
}


/*
 * calculate the divergence of two regions
 * input:
 *      A and B: starting address of region A and region B; note that each region should be continuous in memory 
 *      region_length: size length of each region: num of points in each side
 *      k the nearest neibours need to check in divergence
 *      div_func: divergence function to use:
 * return: divergence  
 */

float get_divs(float *A, float *B, int region_length, int k, int div_func){
    // there are (region_length+1)^2 points in each region
    // we can treat data as if there are in 1-demension
    int i;
    int num_cell = (region_length)*(region_length);

    // distance to the k-nearest neighbours
    float k_dist_a;
    float k_dist_b;

    // u is for gaussian
    float div= 0, tmp = 0;

    // the estimated density 
    float den_a, den_b;

    // start to accumulate the divgence(linear kernel here)
    for(i = 0; i< num_cell; i++){
        // for each point in both regions, get the distance to its k-nearest neighbours
#ifdef debug_1
        printf("\t\ti = %i, A address %p, B address %p, number of cell %d, k= %d\n", i, A, B, num_cell, k);
#endif
        k_dist_a = get_bound_dist(i, A, num_cell, k);
        k_dist_b = get_bound_dist(i, B, num_cell, k);

        // estimate the density
        den_a = k/((4/3)*(num_cell -1)*PI*pow(k_dist_a, 3));
        den_b = k/((4/3)*(num_cell -1)*PI*pow(k_dist_b,3));

#ifdef debug_1
        printf("\tpoint %d dist in lh %lf, dist in rh %lf;density in lh: %lf, in rh: %lf\n", i, k_dist_a,k_dist_b,den_a, den_b);
#endif
        // Now we can free the distance lookup table(matrix)
        // allert if the value is too small, use logorithm here
        if(div_func == 0){
            // linear kernel
            tmp = den_a*den_b; 
            div += tmp;
        }else if(div_func == 1){
            // L_2 divergence
            tmp += pow(den_a - den_b, 2);
        }else if(div_func == 2){
            // KL divergence
            div += den_a*log(den_a/den_b);
        }
    }

    // for L_2 divergence, square root is required
    if(div_func == 1){
        div = sqrt(tmp);
#ifdef debug_1
        printf("get div = sqrt(%lf) = %lf\n", tmp, div);
#endif
    }
    if(!isfinite(div)){
            //fprintf(stderr, "ERR:got infinite div %f, density_a= %f, density_b=%f\n", div, den_a, den_b);
            float noise;
            float noise_range = 0.05;
            noise = (float)rand()/(float)(RAND_MAX/noise_range);
            return (1+1000*noise);
        }
        if(isnan(div)){
            // this can happen when velocity data are just 0
            //fprintf(stderr, "ERR:got nan  div %f", div);
            float noise;
            float noise_range = 0.05;
            noise = (float)rand()/(float)(RAND_MAX/noise_range);
            return (1+1000*noise);
        }
    return div;
}
