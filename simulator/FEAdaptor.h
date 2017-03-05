#ifndef FEADAPTOR_HEADER
#define FEADAPTOR_HEADER

class Attributes;
class Grid;

namespace FEAdaptor
{
  void Initialize(int numScripts, char* scripts[]);

  void Finalize();

  void CoProcess(Grid& grid, Attributes& attributes, double time,
                 unsigned int timeStep, bool lastTimeStep);
}

// this will map cluster id buffer into the same size of pressure buffer
// this function is similar to divide function
// dim is the how many points in each side// 201

static void map_regions(float *cluster_raw, float *cluster_data, const unsigned int dim[3], const unsigned int region_length);

#endif
