#include "FEDataStructures.h"

#include <iostream>
#include <mpi.h>
#include <assert.h>
#include "region_def.h"

Grid::Grid()
{
  this->NumPoints[0] = this->NumPoints[1] = this->NumPoints[2] = 0;
  this->Spacing[0] = this->Spacing[1] = this->Spacing[2] = 0;
}


// generate the grid
// partition happens here
void Grid::Initialize(const unsigned int numPoints[3], const double spacing[3] )
{
  if(numPoints[0] == 0 || numPoints[1] == 0 || numPoints[2] == 0)
    {
    std::cerr << "Must have a non-zero amount of points in each direction.\n";
    }
  for(int i=0;i<3;i++)
    {
    this->NumPoints[i] = numPoints[i];
    this->Spacing[i] = spacing[i];
    }
  /*
  int mpiRank = 0, mpiSize = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

  // feng this is a partion in x direction
  this->Extent[0] = mpiRank*numPoints[0]/mpiSize;
  this->Extent[1] = (mpiRank+1)*numPoints[0]/mpiSize;
  if(mpiSize != mpiRank+1)
    {
    // if it's not the last 
    this->Extent[1]++;
    }
    */
  this->Extent[0] = 0;
  this->Extent[1] = numPoints[0]-1;

  this->Extent[2] = this->Extent[4] = 0;
  this->Extent[3] = numPoints[1] -1;
  this->Extent[5] = numPoints[2] -1;
}


// how many number of points in this partition
unsigned int Grid::GetNumberOfLocalPoints()
{
  return (this->Extent[1]-this->Extent[0]+1)*(this->Extent[3]-this->Extent[2]+1)*
    (this->Extent[5]-this->Extent[4]+1);
}


unsigned int Grid::GetNumberOfLocalCells()
{
  return (this->Extent[1]-this->Extent[0] )*(this->Extent[3]-this->Extent[2])*
    (this->Extent[5]-this->Extent[4] );
}


// feng
// point id: local id
// relative point id in this partition
void Grid::GetLocalPoint(unsigned int pointId, double* point)
{
  unsigned int logicalX = pointId%(this->Extent[1]-this->Extent[0]+1);
  // feng this should be 
  assert(logicalX <= this->Extent[1]);
  point[0] = this->Spacing[0]*logicalX;
  unsigned int logicalY = pointId%((this->Extent[1]-this->Extent[0]+1)*(this->Extent[3]-this->Extent[2]+1));
  logicalY /= this->Extent[1]-this->Extent[0]+1;
  assert(logicalY <= this->Extent[3]);
  point[1] = this->Spacing[1]*logicalY;
  unsigned int logicalZ = pointId/((this->Extent[1]-this->Extent[0]+1)*
                                   (this->Extent[3]-this->Extent[2]+1));
  assert(logicalZ <= this->Extent[5]);
  point[2] = this->Spacing[2]*logicalZ;
}

unsigned int* Grid::GetNumPoints()
{
  return this->NumPoints;
}

unsigned int* Grid::GetExtent()
{
  return this->Extent;
}

double* Grid::GetSpacing()
{
  return this->Spacing;
}

Attributes::Attributes()
{
  this->GridPtr = NULL;
}

void Attributes::Initialize(Grid* grid)
{
  this->GridPtr = grid;
}


// this is only a example
// update points of this partition
void Attributes::UpdateFields(double time)
{
    // how many points in this partition
  unsigned int numPoints = this->GridPtr->GetNumberOfLocalPoints();
  this->Velocity.resize(numPoints*3);
  for(unsigned int pt=0;pt<numPoints;pt++)
    {
    double coord[3];
    
    // coord is the relative position in each partition
    this->GridPtr->GetLocalPoint(pt, coord);

    // velocity is y* time
    this->Velocity[pt] = coord[1]*time;
    }

  // Feng fill the second and third dim  with 0
  std::fill(this->Velocity.begin()+numPoints, this->Velocity.end(), 0.);
  //unsigned int numCells = this->GridPtr->GetNumberOfLocalCells();

  //saved as point data, pressure is always 1
  this->Pressure.resize(numPoints);
  std::fill(this->Pressure.begin(), this->Pressure.end(), 1.);
}

// update points of this partition
void Attributes::UpdateFields(float * vel, float * pres)
{
    // how many points in this partition
  unsigned int numPoints = this->GridPtr->GetNumberOfLocalPoints();

  this->Velocity.resize(numPoints*3);
  this->Pressure.resize(numPoints);
  this->Velocity.assign(vel, vel + 3*numPoints);
  this->Pressure.assign(pres, pres + numPoints);
}

void Attributes::UpdateFields(float * vel, float * pres, float *cluster_data){
    // how many points in this partition
  unsigned int numPoints = this->GridPtr->GetNumberOfLocalPoints();

  this->Velocity.resize(numPoints*3);
  this->Pressure.resize(numPoints);
  this->ClusterId.resize(numPoints);
  this->Velocity.assign(vel, vel + 3*numPoints);
  this->Pressure.assign(pres, pres + numPoints);
  this->ClusterId.assign(cluster_data, cluster_data + numPoints);
}

float* Attributes::GetVelocityArray()
{
  if(this->Velocity.empty())
    {
    return NULL;
    }
  return &this->Velocity[0];
}

float* Attributes::GetPressureArray()
{
  if(this->Pressure.empty())
    {
    return NULL;
    }
  return &this->Pressure[0];
}
float* Attributes::GetClusterIdArray()
{
  if(this->ClusterId.empty())
    {
    return NULL;
    }
  return &this->ClusterId[0];
}
