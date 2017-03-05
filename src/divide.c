# include "divide.h"


void fill_region(float* to, float *from, size_t region_memory_size){
    // add some noise instead of simply memcoy
    int i;
    float noise;
    float noise_range = 0.05;

    // better performance if thunck by thunck
    for(i = 0; i*sizeof(float) < region_memory_size; i++){

        //noise
        noise = (float)rand()/(float)(RAND_MAX/noise_range);

        *to = *from+ noise;
    }
}

void gen_region_samples(float **sample_regions, int num_types, int region_length){
    // allocate
    int i;

    float A[27] = 
    {  
        0.4,0,2,
        0.4,0,1,
        0.4,0,2,
        0.4,0,1,
        0.4,0,0, 
        0.4,0,1,
        0.4,0,2, 
        0.4,0,1,
        0.4,0,2};

    // a sample vortex
    float B[27] = {
        0.2,0.2, 2,
        0.4, 0 ,1,
        0.2, -0.2, 2, 
        0, 0.4,1,
        0 ,0 ,0,
        0, -0.4, 1,
        -0.2 , 0.2, 2,
        -0.4, 0, 1,
        -0.2, -0.2, 2};

    // also a turbulence
    float C[27];
    for(i = 0; i < 27 ; i++){
        if(i%3!=2){
            C[i] = -B[i];
        }
        else{ 
            C[i] = B[i];
        }
    }

    // another pattern
    float D[27] = {
        0.2,0.2, 2,
        0.4, 0 ,1,
        0.2, -0.2, 2,
        0.4, 0,1,
        0.4 ,0 ,0,
        0.4, 0, 1,
        0.2 , -0.2, 2,
        0.4, 0, 1,
        0.2, 0.2, 2};


    float * cur_region;
    size_t region_memory_size = (region_length+1)*(region_length+1)*3*sizeof(float);

    
    for(i =0; i < num_types; i++){
        sample_regions[i] = (float *)malloc(region_memory_size);
        if(sample_regions[i] == NULL){
            printf("sample region No. %d  allocate error\n", i);
            exit(-1);
        }
    }

    /*
     * type 0: one-direction flow
     */
    cur_region = sample_regions[0];
    memcpy(cur_region, A, region_memory_size);

    /* follow this for more complexed region
    for(ii = 0; ii < region_length+1; ii++){
        for(jj = 0; jj < region_length+1; jj++){
            // for each point inside the region

            ux = MAX_VELOCITY;
            uy = MAX_VELOCITY;

            // square of dist to center
            dist2 = (ii - region_length/2.0)*(ii - region_length/2.0) + (jj - region_length/2.0)*(jj - region_length/2.0); 

            // access this triple
            float* index_tuple = cur_region + (ii*(length+1) +jj)*d3;

            *(index_tuple + 0) = ux;
            *(index_tuple + 1) = uy;
            *(index_tuple + 2) = dist2;
        }
    }
    */



    /*
     * type 1: circles
     */

    
    cur_region = sample_regions[1];

    memcpy(cur_region, B, region_memory_size);

    /*
     * type 2: division in center
     */

    cur_region = sample_regions[2];

    memcpy(cur_region, D, region_memory_size);
}

void free_region_samples(float **sample_regions, int num_types){
    int i;
    for(i = 0; i< num_types; i++){
        if(sample_regions[i]){
            free(sample_regions[i]);
            printf("sample region %d is freed\n", i);
        }
    }
    if(sample_regions){
        free(sample_regions);
        printf("all sample regions are freed");
    }
}



void divide_synthetic(int dim, int region_length, int *p_num_region, float **p_regions){
    int p,q, side_num_region;


    if((dim - 1)%region_length != 0){
        printf("not pefect division, try different region size of datacut size\n");
        exit(-1);
    }

    // actual num_regions will be square of this
    side_num_region = (dim -1)/region_length; // this will be (201-1)/10 = 20
    *p_num_region = side_num_region* side_num_region;

    int d1 = *p_num_region;

    size_t region_memory_size = (region_length+1)*(region_length+1)*3*sizeof(float);
    
    int type_region;
    

    float *this_region;

    // alocate space for regions
    float *regions = (float *)malloc( d1 * region_memory_size);
    if(regions == NULL){
        perror("    allocate space for regions");
        exit(-1);
    }else{
        printf("    %d x %ld  bytes space is allocated to region\n", d1, region_memory_size );
    }

    // this is how we divide all the region areas
    int index_1_3 = side_num_region/3.0;//6  
    int index_2_3 = 2*index_1_3;
    fprintf(stderr, "two index: %d and %d\n", index_1_3, index_2_3);


    // generate all the sample regions
    float **sample_regions;
    int num_types = 3;

    sample_regions = (float **)malloc(sizeof(float *)* num_types);
    if(sample_regions == NULL){
        printf("sample region allocate error\n");
        exit(-1);
    }
    else{
        printf("sample regions address saved in %p\n", sample_regions);
    }

    gen_region_samples(sample_regions, num_types, region_length);

    this_region = regions;

    int count =0;
    for(p = 0; p < side_num_region; p++){

        //assign regions for sythetic data
        //types of regions:
        //   regions have three types;
        //  0. one direction
        //          -> -> ->
        //          -> -> ->
        //          -> -> ->
        //  1. circle
        //         
        //  2. division in center but still almost in one direction



        if(p < index_1_3){
            fprintf(stderr, "row %d is type 0, one direction\n", p);
            type_region = 0;}
        else if(p < index_2_3){
            fprintf(stderr, "row %d is type 1, circle\n", p);
            type_region = 1;}
        else{
            fprintf(stderr, "row %d is type 2, division\n", p);
            type_region = 2;
        }
        for( q = 0; q < side_num_region; q++){
            // for each region
            // this region will have corner:
            //  (pl, ql) , (pl, (q+1)l)
            //  ((p+1)l,ql), ((p+1)l,(q+1)l)
            //  also the region all have the center (pl+l/2, ql+l/2)

            // fill different patterns in each region
            memcpy(this_region, sample_regions[type_region], region_memory_size);
            // I need to add some noise here:
            fill_region(this_region, sample_regions[type_region], region_memory_size);
            //
            fprintf(stderr, "No.%d region has type %d\n", count++, type_region);
            
           //get the start address of this region
            this_region += (region_length+1)*(region_length+1)*3 ;
        }
    }
    free_region_samples(sample_regions, num_types);

    *p_regions = regions;

}

/*
 *  a quick demonstration:
 *  before mapping:
 *      [200][400][3], row and colomn dimension are 200 and 400, each position has vx, vy, vz
 *
 *      use region size 20x20
 *  after mapping:
 *      [10][20][400][3], there will be 10*20 regions, each region has 400 points, each points will have vx, vy, and distance to center
 *
 *
 *  for simulation generated data
 *  before mapping
 *      slowest in the x and fastest in the z directions
 *      vel data will be aligned in format as: vx0, vy0, vz0, vx1, vy1, vz1
 *  after mapping
 *      
 */
void divide(float *pdata, int dims[3], int l, int *p_num_region, float *buffer_region){
    int p,q,ii,jj;
    //side_num_region;

    // index inside the logic block 
    int index_x, index_y;

    // index if translated into the linear bubuffer
    int linear_index;
    
    // velocity of dimension x and y
    float ux, uy, *tripple_address, dist2;

    if((dims[0])%l != 0|| dims[1]%l != 0){
        printf("not pefect division, try different region size of datacut size\n");
        exit(-1);
    }

    // actual num_regions will be square of this
    int num_region_row = dims[0]/l;
    int num_region_col =  dims[1]/l;
    *p_num_region = dims[0]*dims[1]/(l*l); // this will be (201-1)/10 = 20

    //int d1 = *p_num_region;
    int d2 = (l)*(l);
    int d3 = 3;

    float *regions = buffer_region;
    

    for(p = 0; p < num_region_row; p++){
        for( q = 0; q < num_region_col; q++){
            // for each region
            // this region will have corner:
            //  (pl, ql) , (pl, (q+1)l)
            //  ((p+1)l,ql), ((p+1)l,(q+1)l)
            //  also the region all have the center (pl+l/2, ql+l/2)
            for(ii = 0; ii < l; ii++){
                for(jj = 0; jj < l; jj++){
                    // for each point inside the region
                    index_x = p*l+ii;
                    index_y = q*l+jj;

                    // mapped the logic address into linear address
                    linear_index = index_x*dims[1] + index_y;
                    tripple_address = pdata+ 3*linear_index;
                    ux = *tripple_address;
                    uy = *(tripple_address + 1);

                    // square of dist to center
                    dist2 = (ii - (l-1)/2.0)*(ii - (l-1)/2.0) + (jj - (l-1)/2.0)*(jj - (l-1)/2.0); 

                    // add this tripple into region

                    // inefficient
                    int shift_tuple = (p*num_region_col+q)*d2*d3 + (ii*(l) +jj)*d3;
                    float *index_tuple = regions + shift_tuple;

                    *(index_tuple + 0) = ux;
                    *(index_tuple + 1) = uy;
                    *(index_tuple + 2) = dist2;
                }
            }
            //printf("    region %d-%d is tupled\n", p, q);
        }
    }
}
