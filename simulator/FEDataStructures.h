// define grid info and field info
//  toy simulation definition is also from here

#ifndef FEDATASTRUCTURES_HEADER
#define FEDATASTRUCTURES_HEADER

#include <cstddef>
#include <vector>

// this will store all the local information in each partition
class Grid
{
public:
  Grid();
  void Initialize(const unsigned int numPoints[3], const double spacing[3]);
  // still need to provide partion info
  unsigned int GetNumberOfLocalPoints();
  unsigned int GetNumberOfLocalCells();
  void GetLocalPoint(unsigned int pointId, double* point);
  unsigned int* GetNumPoints();
  unsigned int* GetExtent();
  double* GetSpacing();
private:
  unsigned int RegionLength;
  unsigned int NumPoints[3];

  // they should use same grid but for regions there I need zoom to fit large regions
  unsigned int Extent[6];
  double Spacing[3];
};

// this is for each partition
class Attributes
{
// A class for generating and storing point and cell fields.
// Velocity is stored at the points and pressure is stored
// for the cells. The current velocity profile is for a
// shearing flow with U(y,t) = y*t, V = 0 and W = 0.
// Pressure is constant through the domain.
public:
  Attributes();
  void Initialize(Grid* grid);
  void UpdateFields(double time);
  void UpdateFields(float *vel, float*pres);
  void UpdateFields(float *vel, float*pres, float*cluster_data);
  // commented by feng
  // double* GetVelocityArray();
 
  float* GetVelocityArray();
  float* GetPressureArray();
  float* GetClusterIdArray();

private:
  std::vector<float> Velocity;
  std::vector<float> Pressure;

  // this will have smaller size
  std::vector<float> ClusterId;
  Grid* GridPtr;
};
#endif
