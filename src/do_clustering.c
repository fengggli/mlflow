#include "do_clustering.h"

int main(){
    char *dist_path = "data/all_dist.txt";
    char *output_path = "data/clusterid.txt";
    int nclusters, num_region,npass, ifound;
    double error;

    nclusters = 3;
    num_region = 400; 
    npass =100;

    //int clusterid[num_region];

    //call the divergence function
    run_clustering (nclusters,num_region, dist_path,npass,output_path, &error,&ifound);
    printf("clustering completed, result  is stored in %s\n", output_path);
    printf("error is %.3lf, num of runs is %d\n", error, ifound);

    return 0;
}
