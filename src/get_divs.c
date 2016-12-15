#include "get_divs.h"
//#define debug_1

// parallel version and sequential version have different divergence
#define DEBUG_REGION_MAPPING

/*
 * this is basically a sort function
 * insertion sort here
 * input/ output: array_to_sort
 */
void get_kth_dist(float *array_to_sort, int length_to_sort, int k){
    int i,j;
    float key;
    for(j = 1; j< length_to_sort; j++){
        key = array_to_sort[j];
        for(i=j-1;array_to_sort[i] > key&i>=0;i--){
            array_to_sort[i+1] = array_to_sort[i];
        }
        array_to_sort[i+1] = key;
    }
#ifdef debug
    print_matrix(array_to_sort, 1, length_to_sort);
#endif
    //return array_to_sort[k-1];
}

void print_matrix(float *Dismatrix, int num_row, int num_col){
    printf("print matrix\n");
    int i,j;
    for(i = 0; i < num_row; i++){
        for(j = 0; j < num_col; j++){
            printf("%.3f\t", Dismatrix[i*num_col + j]);
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
float get_bound_dist(int i, float *DisMatrix, int length, int k){
    // actually this is a task to sort all the elements in the i-th row in DisMatrix

    // the index of k-th least distance
    int m,k_dist;

    // start point and length of array to sort
    // bug fixed here: used to modify origin matrix!
    float tmp_array[length];
    for(m = 0; m< length; m++)
        tmp_array[m] = *(DisMatrix+ i*length +m);

    //float *array_to_sort = DisMatrix + i*length;
    int length_to_sort = length;
    get_kth_dist(tmp_array, length_to_sort, k);

    // if all the tipple values are the same, nearest neibour distance can be 0
    return tmp_array[k-1];
}


/*
 * calculate the divergence of two regions
 * input:
 *      A and B: starting address of region A and region B; note that each region should be continuous in memory 
 *      region_length: size length of each region: num of points in each side -1
 *      k the nearest neibours need to check in divergence
 *      div_func: divergence function to use:
 * return: divergence  
 */

float get_divs(float *A, float *B, int region_length, int k, int div_func){
    // there are (region_length+1)^2 points in each region
    // we can treat data as if there are in 1-demension
    int i, j;
    int num_cell = (region_length+1)*(region_length +1);
    //printf("there are %d points in each region\n", num_cell);
    //fprintf(stderr, "there are %d points in each region\n", num_cell);

    float* DisMatrixA = (float *)malloc(sizeof(float)*num_cell*num_cell);
    float* DisMatrixB = (float *)malloc(sizeof(float)*num_cell*num_cell);

    // distance to the k-nearest neighbours
    float dist_tmp = 0;
    float k_dist_a;
    float k_dist_b;

    // u is for gaussian
    float div= 0, tmp = 0, tmp2, u;

    // the estimated density 
    float den_a, den_b;


#ifdef debug_1
    printf("matrix A:\n");
#endif

    // get all the pair-wise distances between all the points 
    for(i = 0; i < num_cell; i++){
        // bug here, but why!!!!!
#ifdef debug_1
        printf("\t point %d: (%.3f %.3f %.3f)\n", i, *(A + 3*i+ 0), *(A +3*i + 1),*(A + 3*i +2)); 
#endif

        for(j = i; j < num_cell; j++){
            if(j == i) 
                DisMatrixA[i*num_cell+j] = LARGE_NUMBER ;
            else{
                dist_tmp = (pow(*(A + 3*i) - *(A + 3*j),2) + pow(*(A + 3*i +1) - *(A + 3*j +1), 2)+ pow(*(A+3*i+2) - *(A+3*j+2),2) );
                DisMatrixA[i* num_cell+j] = dist_tmp;
                DisMatrixA[j* num_cell+i] = dist_tmp; // fixed the bug k_dist_a == 0
            }
            //printf("A: distance(%d, %d)is %.4f\n",i,j, DisMatrixA[i*num_cell +j]);
            //print_matrix(DisMatrixA, num_cell);
        }
    }

    // also get distances in region B
#ifdef debug
    printf("distance region matrix at %p is calculated\n", (void *)A);
#endif
    for(i = 0; i < num_cell; i++)
        for(j = i; j < num_cell; j++){
            if(j == i) 
                DisMatrixB[i*num_cell+j] = LARGE_NUMBER ;
            else{
                dist_tmp = (pow(*(B + 3*i) - *(B + 3*j),2) + pow(*(B + 3*i +1) - *(B + 3*j +1), 2)+ pow(*(B+3*i+2) - *(B+3*j+2),2) );
                DisMatrixB[i* num_cell+j] = dist_tmp;
                DisMatrixB[j* num_cell+i] = dist_tmp;
            }
            //printf("B: distance(%d, %d)is %.4f\n",i,j, DisMatrixB[i*num_cell +j]);
        }

#ifdef debug
    printf("distance matrix region at %p is calculated\n", (void *)(B));
    print_matrix(DisMatrixA, num_cell, num_cell);
    print_matrix(DisMatrixB, num_cell, num_cell);
#endif

    // start to accumulate the divgence(linear kernel here)
    for(i = 0; i< num_cell; i++){
        // for each point in both regions, get the distance to its k-nearest neighbours
        k_dist_a = get_bound_dist(i, DisMatrixA, num_cell, k);
        k_dist_b = get_bound_dist(i, DisMatrixB, num_cell, k);

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

#ifdef debug
    printf("div between regions at %p and %p is calculated\n", (void*)(A), (void *)(B));
#endif

    // for L_2 divergence, square root is required
    if(div_func == 1){
        div = sqrt(tmp);
    }
    free(DisMatrixA);
    free(DisMatrixB);
    return div;
}
