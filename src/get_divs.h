/*
 * Sep 7 by Feng Li, IUPUI
 * calculate the divergence of two square 2-d regions
 */
#include <stdio.h>
#include <math.h>

/*
 * this is basically a sort function
 * use bubble sort since we only need the kth-least distance
 */
int k_index = get_kth_index(double *array_to_sort, int length_to_sort, int k){
    // bubble sort here
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
double k_dist = get_bound_dist(int i, double * DisMatrix, int length, int k){
    // actually this is a task to sort all the elements in the i-th row in DisMatrix

    // the index of k-th least distance
    int k_index;

    // start point and length of array to sort
    double *array_to_sort = DisMatrix + i*length;
    int length_to_sort = length;
    k_index = get_kth_index(array_to_sort, length_to_sort, k);

    // then return that distance
    return *(array_to_sort + k_index);
}

/*
 * calculate the divergence of two regions
 * input:
 *      region A and region B; note that each region should be continuous in memory 
 *      length of region
 * return: divergence  
 */
double get_divs(double (*A)[3], double (*B)[3], int region_length){
    // there are region_length^2 points in each region
    // we can treat data as if there are in 1-demension
    int i, j;
    int num_elem = region_length^2;
    double DisMatrixA = (double *) malloc(sizeof(double)*elem_num*elem_num);
    double DisMatrixB = (double *) malloc(sizeof(double)*elem_num*elem_num);

    // distance to the k-nearest neighbours
    double k_dist_a;
    double k_dist_b;

    double div= 0, tmp;

    // the estimated density 
    double den_a, den_b;

    // get all the pair-wise distances between all the points 
    for(i = 0; i < num_elem; i++)
        for(j = i+1; j < region_length^2; j++)
            DisMatrixA[i* num_elem+j] = sqrt((A[i][0] - A[j][0]^2) +(A[i][1] - A[j][1]^2)+ (A[i][2] - A[j][2]^2) );

    // also get distances in region B
    for(i = 0; i < num_elem; i++)
        for(j = i+1; j < region_length^2; j++)
            DisMatrixB[i* num_elem+j] = sqrt((B[i][0] - B[j][0]^2) +(B[i][1] - B[j][1]^2)+ (B[i][2] - B[j][2]^2) );

    // start to accumulate the divgence(linear kernel here)
    for(i = 0; i< num_elem; i++i){
        // for each point in both regions, get the distance to its k-nearest neighbours
        k_dist_a = get_bound_dist(i, DisMatrixA, elem_num, k);
        k_dist_b = get_bound_dist(i, DisMatrixB, elem_num, k);

        // estimate the density
        den_a = k/((4/3)*(num_elem -1)*PI*k_dist_a^3);
        den_b = k/((4/3)*(num_elem)*PI*k_dist_b^3);
        
        // Now we can free the distance lookup table(matrix)
        // allert if the value is too small, use logorithm here
        tmp = den_a*den_b; 
        /*
         * if tmp is too small give some allerts
        if(tmp <<)
        */
        div +=tmp;
        }
    return div;
}

               
