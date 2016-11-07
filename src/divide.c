# include "divide.h"

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
        perror("    allocate space for regions");
        exit(-1);
    }else{
        printf("    %d x %d x %d x sizeof(float) space is allocated to region\n", d1, d2, d3);
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

                    int shift_tuple = (p*side_num_region+q)*d2*d3 + (ii*(l+1) +jj)*d3;
                    float *index_tuple = regions + shift_tuple;


                    *(index_tuple + 0) = ux;
                    *(index_tuple + 1) = uy;
                    *(index_tuple + 2) = dist2;
                }
            }
            //printf("    region %d-%d is tupled\n", p, q);
        }
    }
    *p_regions = regions;
}
