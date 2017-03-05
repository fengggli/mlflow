#include "consumer.h"

// fake one, only select first num_elem_sample
void prepare_sampled_buffer(float *buffer_region, float *buffer_sample, int num_elems_region, int num_elems_sample, int region_length){
    int i;
    int spacing = region_length*region_length*3;
    for(i = 0; i< num_elems_sample; i++){
        memcpy(buffer_region+ spacing, buffer_sample + spacing, spacing*sizeof(float));
    }
}

void assign_clusterid(float *buffer_region,int num_region, float * buffer_sample_all, int num_sample_all,int *buffer_medoids,int k,int region_length,int k_npdiv, int div_func,  float *clusterids){
    int i, j;
    float clusterid;
    float divs, divs_tmp;
    float *buffer_a, *buffer_b;

    int spacing = region_length*region_length*3;

    for(i = 0; i< num_region; i++){
        // get divergence to each medoids
        clusterid = -1;
        divs = -1;
        buffer_a = buffer_region + i*spacing;

        for( j = 0; j <k; j++){
            // compare this region with one medoids
            buffer_b = buffer_sample_all + buffer_medoids[j]*spacing;
            divs_tmp = get_divs( buffer_a , buffer_b, region_length, k_npdiv, div_func);
            if(divs_tmp < divs){
                divs =divs_tmp;
                clusterid = j;
            }
            if(divs_tmp <0){
                printf("negative divergence !\n");
                exit(-1);
            }
        }
        clusterids[i] = clusterid;
    }

}

