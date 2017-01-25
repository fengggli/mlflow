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
  unsigned int GetNumberOfLocalPoints();
  unsigned int GetNumberOfLocalCells();
  void GetLocalPoint(unsigned int pointId, double* point);
  unsigned int* GetNumPoints();
  unsigned int* GetExtent();
  double* GetSpacing();
private:
  unsigned int NumPoints[3];
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
  // commented by feng
  // double* GetVelocityArray();
 
  float* GetVelocityArray();
  float* GetPressureArray();

private:
  std::vector<float> Velocity;
  std::vector<float> Pressure;
  Grid* GridPtr;
};
#endif
