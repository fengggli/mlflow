#ifndef DO_CLUSTERING_H
#define DO_CLUSTERING_H
/*
 * Sep 10, by Feng Li, IUPUI
 * based on the divergences beteen the regions, do the clustering job
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cluster.h"

/*
 *  kmedoid methods is implemented in cluster.h
 *
 */
void run_clustering (int nclusters, int num_region, char *dist_path,
  int npass, char *output_path, double* error, int* ifound);

#endif
