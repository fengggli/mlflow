/*
 * Sep 7 by Feng Li, IUPUI
 * calculate the divergence of two square 2-d regions
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI (3.1415926)
#define LARGE_NUMBER (999999)
//#define debug

/*
 * this is basically a sort function
 * insertion sort here
 */
void get_kth_dist(double *array_to_sort, int length_to_sort, int k){
    int i,j;
    double key;
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


/* print the contant of matrix
 */
print_matrix(double *Dismatrix, int num_row, int num_col){
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
double get_bound_dist(int i, double *DisMatrix, int length, int k){
    // actually this is a task to sort all the elements in the i-th row in DisMatrix

    // the index of k-th least distance
    int m,k_dist;

    // start point and length of array to sort
    // bug here: try to modify origin matrix!
    double tmp_array[length];
    for(m = 0; m< length; m++)
        tmp_array[m] = *(DisMatrix+ i*length +m);

    //double *array_to_sort = DisMatrix + i*length;
    int length_to_sort = length;
    get_kth_dist(tmp_array, length_to_sort, k);

    // if all the tipple values are the same, nearest neibour distance can be 0
    return tmp_array[k-1];
}

/*
 * calculate the divergence of two regions
 * input:
 *      region A and region B; note that each region should be continuous in memory 
 *      length of region
 *      k the nearest neibours need to check in divergence
 * return: divergence  
 */
double get_divs(double (*A)[3], double (*B)[3], int region_length, int k, int div_func){
    // there are region_length^2 points in each region
    // we can treat data as if there are in 1-demension
    int i, j;
    int num_elem = pow(region_length,2);
    double* DisMatrixA = (double *)malloc(sizeof(double)*num_elem*num_elem);
    double* DisMatrixB = (double *)malloc(sizeof(double)*num_elem*num_elem);

    // distance to the k-nearest neighbours
    double dist_tmp = 0;
    double k_dist_a;
    double k_dist_b;

    // u is for gaussian
    double div= 0, tmp = 0, tmp2, u;

    // the estimated density 
    double den_a, den_b;

    // get all the pair-wise distances between all the points 
    for(i = 0; i < num_elem; i++)
        for(j = i; j < num_elem; j++){
            if(j == i) 
                DisMatrixA[i*num_elem+j] = LARGE_NUMBER ;
            else{
                dist_tmp = (pow(A[i][0] - A[j][0],2) + pow(A[i][1] - A[j][1], 2)+ pow(A[i][2] - A[j][2],2) );
                DisMatrixA[i* num_elem+j] = dist_tmp;
                DisMatrixA[j* num_elem+i] = dist_tmp; // fixed the bug k_dist_a == 0
            }
            //printf("A: distance(%d, %d)is %.4f\n",i,j, DisMatrixA[i*num_elem +j]);
            //print_matrix(DisMatrixA, num_elem);
        }

    // also get distances in region B
    for(i = 0; i < num_elem; i++)
        for(j = i; j < num_elem; j++){
            if(j == i) 
                DisMatrixB[i*num_elem+j] = LARGE_NUMBER ;
            else{
                dist_tmp = (pow(B[i][0] - B[j][0],2) + pow(B[i][1] - B[j][1],2)+ pow(B[i][2] - B[j][2], 2) );
                DisMatrixB[i* num_elem+j] = dist_tmp;
                DisMatrixB[j* num_elem+i] = dist_tmp;
            }
            //printf("B: distance(%d, %d)is %.4f\n",i,j, DisMatrixB[i*num_elem +j]);
        }
#ifdef debug
    print_matrix(DisMatrixA, num_elem, num_elem);
    print_matrix(DisMatrixB, num_elem, num_elem);
#endif

    // start to accumulate the divgence(linear kernel here)
    for(i = 0; i< num_elem; i++){
        // for each point in both regions, get the distance to its k-nearest neighbours
        k_dist_a = get_bound_dist(i, DisMatrixA, num_elem, k);
        k_dist_b = get_bound_dist(i, DisMatrixB, num_elem, k);

        // estimate the density
        den_a = k/((4/3)*(num_elem -1)*PI*pow(k_dist_a, 3));
        den_b = k/((4/3)*(num_elem)*PI*pow(k_dist_b,3));

//#ifdef debug
        printf("\tpoint %d dist in lh %lf, dist in rh %lf;density in lh: %lf, in rh: %lf\n", i, k_dist_a,k_dist_b,den_a, den_b);
//#endif
        
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
    }
    free(DisMatrixA);
    free(DisMatrixB);
    return div;
}

