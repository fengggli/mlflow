#include "consumer.h"
#include <float.h>



void assign_clusterid(float *buffer_region,int num_region, float * buffer_sample_all, int num_sample_all,int *buffer_medoids,int k,int region_length,int k_npdiv, int div_func,  float *clusterids){
    int i, j;
    float clusterid;
    float divs, divs_tmp;
    float *buffer_a, *buffer_b;

    int spacing = region_length*region_length*3;

    for(i = 0; i< num_region; i++){
        // get divergence to each medoids
        clusterid = -1;
        divs = FLT_MAX;
        buffer_a = buffer_region + i*spacing;

        for( j = 0; j <k; j++){
            // compare this region with one medoids
            buffer_b = buffer_sample_all + buffer_medoids[j]*spacing;
            divs_tmp = get_divs( buffer_a , buffer_b, region_length, k_npdiv, div_func);
            //printf("divergence to medoids %d is %f\n", buffer_medoids[j], divs_tmp);
            if(divs_tmp < divs){
                divs =divs_tmp;
                clusterid = j;
            }
        }
        //printf("clusterid for region %d is %f\n", i, clusterid);
        clusterids[i] = clusterid;
    }
}

