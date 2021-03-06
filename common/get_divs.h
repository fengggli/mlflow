#ifndef GET_DIVS_H
#define GET_DIVS_H
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

void print_matrix(float *Dismatrix, int l, int x0, int y0, int x1, int y1);

// get the k-th smallest value from size length vector
float  get_k_order(float *tmp_array, int length, int k);
        

float get_bound_dist(int i, float *DisMatrix, int length, int k);

float get_divs(float *A, float *B, int region_length, int k, int div_func);

#endif

