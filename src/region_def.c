#include "region_def.h"

void fill_region_def(Region_Def *p_region_def){
    // some pre-definition
    p_region_def->region_length = REGION_LENGTH;
    p_region_def->side_num_region = (POINTS_SIDE - 1)/REGION_LENGTH;
    p_region_def-> num_region =p_region_def->side_num_region*p_region_def->side_num_region ;

    p_region_def->region_num_cell = (p_region_def->region_length+1)*(p_region_def->region_length+1);
    // attention here: this is used only when malloc
    p_region_def->region_memory_size = p_region_def->region_num_cell*3*sizeof(float);
}

void extract_region_def(Region_Def *p_region_def, int *p_region_length,int * p_side_num_region, int *p_num_region,int * p_region_num_cell, size_t *p_region_memory_size){
    *p_region_length = p_region_def->region_length;
    *p_side_num_region = p_region_def->side_num_region;
    *p_num_region = p_region_def->num_region ;
    *p_region_num_cell = p_region_def->region_num_cell;
    // attention here: this is used only when malloc
    *p_region_memory_size = p_region_def->region_memory_size;
}
