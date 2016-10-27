#include "do_clustering.h"
#include "run.h"
#include <stdlib.h>

int main(){
    //char *dist_path = "data/all_dist.txt";
    //char *output_path = "data/clusterid.txt";
    char *dist_path = "data/all_dist_201_lg.txt";
    char *output_path = "data/clusterid_201_lg_1.txt";
    int nclusters, num_region,npass, ifound;
    double error;

    // timer
    double time_start, time_current;

    nclusters = 3;
    num_region = 400; 
    npass =100;

    //int clusterid[num_region];


    time_start = get_cur_time();

    //call the divergence function
    printf("read divergence from %s, num_regions is %d\n", dist_path, num_region );
    run_clustering (nclusters,num_region, dist_path,npass,output_path, &error,&ifound);

    time_current = get_cur_time();
    printf("clustering completed, result  is stored in %s\n", output_path);
    printf("%.4f time is used for clustering\n", time_current - time_start );
    printf("error is %.3lf, num of runs is %d\n", error, ifound);

    return 0;
}
