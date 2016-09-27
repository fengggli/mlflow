#ifndef DIVIDE_H 
#define DIVIDE_H 

#include <stdio.h>
#include <stdlib.h>
/*
 * Sep 7, by Feng Li, IUPUI
 * Devide a slice of turbulence data into regions
 * before splitting, we have a block of a*b*1, default value of c is 1
 * each region is a block with a length of l, ovelapping length default value is 0
 */


/*
 * input
 *  pdata, pointer to data buffer
 *  dim: num of points in x and y dimension 
 *  l: region_length: for 11 points in each size, the length will be 10
 * output:
 *  p_num_region; number of regions
 *  p_regions: the tripple features(vx, vy, dc). there will be num_regions* points_in_each_region*3 float numbers,note that features are stored in flat format.
 */
void divide(float *pdata, int dim, int l, int *p_num_region, float **p_regions){
    int p,q,ii,jj, side_num_region;

    // index inside the logic block 
    int index_x, index_y;

    // index if translated into the linear bubuffer
    int linear_index;
    
    // velocity of dimension x and y
    float ux, uy, *tripple_address, dist2;

    if((dim - 1)%l != 0){
        printf("not pefect division, try different region size of datacut size\n");
        exit(-1);
    }

    // actual num_regions will be square of this
    side_num_region = (dim -1)/l;
    *p_num_region = side_num_region* side_num_region;

    int d1 = *p_num_region;
    int d2 = (l+1)*(l+1);
    int d3 = 3;

    float *regions = (float *)malloc(d1*d2*d3*sizeof(float));
    if(regions == NULL){
        perror("allocate space for regions");
        exit(-1);
    }else{
        printf("%d x %d x %d space is allocated to region\n", d1, d2, d3);
    }

    for(p = 0; p < side_num_region; p++){
        for( q = 0; q < side_num_region; q++){
            // for each region
            // this region will have corner:
            //  (pl, ql) , (pl, (q+1)l)
            //  ((p+1)l,ql), ((p+1)l,(q+1)l)
            //  also the region all have the center (pl+l/2, ql+l/2)
            for(ii = 0; ii < l+1; ii++){
                for(jj = 0; jj < l+1; jj++){
                    // for each point inside the region
                    index_x = p*l+ii;
                    index_y = q*l+jj;

                    // mapped the logic address into linear address
                    linear_index = index_x*dim + index_y;
                    tripple_address = pdata+ 3*linear_index;
                    ux = *tripple_address;
                    uy = *(tripple_address + 1);

                    // square of dist to center
                    dist2 = (ii - l/2.0)*(ii - l/2.0) + (jj - l/2.0)*(jj - l/2.0); 

                    // add this tripple into region
                    *(regions + (p*side_num_region+q)*d2*d3 + (ii*(l+1) +jj)*d3 + 0) = ux;
                    *(regions + (p*side_num_region+q)*d2*d3 + (ii*(l+1) +jj)*d3 + 1) = uy;
                    *(regions + (p*side_num_region+q)*d2*d3 + (ii*(l+1) +jj)*d3 + 2) = dist2;
                }
            }
        }
    }
    *p_regions = regions;
}

#endif

