#include <iostream>
#include "FEAdaptor.h"
#include "FEDataStructures.h"

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>

#include "region_def.h"

namespace
{
  vtkCPProcessor* Processor = NULL;
  // I need two grid here
  vtkImageData* VTKGrid = NULL;

  void BuildVTKGrid(Grid& grid)
  {
    // The grid structure isn't changing so we only build it
    // the first time it's needed. If we needed the memory
    // we could delete it and rebuild as necessary.
    if(VTKGrid == NULL)
      {
      VTKGrid = vtkImageData::New();
      int extent[6];
      for(int i=0;i<6;i++)
        {
        extent[i] = grid.GetExtent()[i];
        }
      VTKGrid->SetExtent(extent);
      VTKGrid->SetSpacing(grid.GetSpacing());

      }
  }

  void UpdateVTKAttributes(Grid& grid, Attributes& attributes)
  {
    if(VTKGrid->GetPointData()->GetNumberOfArrays() == 0)
      {
      // velocity array
      vtkNew<vtkFloatArray> velocity;
      velocity->SetName("velocity");
      velocity->SetNumberOfComponents(3);
      velocity->SetNumberOfTuples(static_cast<vtkIdType>(grid.GetNumberOfLocalPoints()));
      VTKGrid->GetPointData()->AddArray(velocity.GetPointer());


      // pressure array
      vtkNew<vtkFloatArray> pressure;
      pressure->SetName("pressure");
      pressure->SetNumberOfComponents(1);
      VTKGrid->GetPointData()->AddArray(pressure.GetPointer());

#ifdef INCLUDE_ML
      // cluster array
      vtkNew<vtkFloatArray> cluster;
      cluster->SetName("cluster");
      cluster->SetNumberOfComponents(1);
      VTKGrid->GetPointData()->AddArray(cluster.GetPointer());

#endif

      }

    /* commented by feng li, store all info in cell array
    if(VTKGrid->GetCellData()->GetNumberOfArrays() == 0)
      {
      // pressure array
      vtkNew<vtkFloatArray> pressure;
      pressure->SetName("pressure");
      pressure->SetNumberOfComponents(1);
      //VTKGrid->GetCellData()->AddArray(pressure.GetPointer());
      }
      */
    vtkFloatArray* velocity = vtkFloatArray::SafeDownCast(
      VTKGrid->GetPointData()->GetArray("velocity"));

    // feng: I will store velocity attributes as vtk format instead

    // The velocity array is ordered as vx0,vx1,vx2,..,vy0,vy1,vy2,..,vz0,vz1,vz2,..
    // so we need to create a full copy of it with VTK's ordering of
    // vx0,vy0,vz0,vx1,vy1,vz1,..

    // feng: this is the data from
    float* velocityData = attributes.GetVelocityArray();

    vtkIdType numTuples = velocity->GetNumberOfTuples();
    for(vtkIdType i=0;i<numTuples;i++)
      {
      float values[3] = {velocityData[3*i], velocityData[3*i+1],
                          velocityData[3*i+2]};
      velocity->SetTupleValue(i, values);
      //velocity->SetTypedTuple(i, values);
      }

    vtkFloatArray* pressure = vtkFloatArray::SafeDownCast(
      VTKGrid->GetPointData()->GetArray("pressure"));
    // The pressure array is a scalar array so we can reuse
    // memory as long as we ordered the points properly.
    float* pressureData = attributes.GetPressureArray();
    pressure->SetArray(pressureData, static_cast<vtkIdType>(grid.GetNumberOfLocalPoints()), 1);

#ifdef INCLUDE_ML
    vtkFloatArray* cluster = vtkFloatArray::SafeDownCast(
      VTKGrid->GetPointData()->GetArray("cluster"));
    // The cluster array is a scalar array so we can reuse
    // memory as long as we ordered the points properly.
    // need zoom here
    float* clusterData = attributes.GetClusterIdArray();
    cluster->SetArray(clusterData, static_cast<vtkIdType>(grid.GetNumberOfLocalPoints()), 1);
#endif
  }

  void BuildVTKDataStructures(Grid& grid, Attributes& attributes)
  {
    BuildVTKGrid(grid);
    UpdateVTKAttributes(grid, attributes);
  }
}

namespace FEAdaptor
{

  void Initialize(int numScripts, char* scripts[])
  {
    if(Processor == NULL)
      {
      Processor = vtkCPProcessor::New();
      Processor->Initialize();
      }
    else
      {
      Processor->RemoveAllPipelines();
      }
    for(int i=1;i<numScripts;i++)
      {
      vtkNew<vtkCPPythonScriptPipeline> pipeline;
      pipeline->Initialize(scripts[i]);
      Processor->AddPipeline(pipeline.GetPointer());
      }
  }

  void Finalize()
  {
    if(Processor)
      {
      Processor->Delete();
      Processor = NULL;
      }
    if(VTKGrid)
      {
      VTKGrid->Delete();
      VTKGrid = NULL;
      }
  }

  void CoProcess(Grid& grid, Attributes& attributes, double time,
                 unsigned int timeStep, bool lastTimeStep)
  {

     // grid is defined outside
    vtkNew<vtkCPDataDescription> dataDescription;
    
    // add input
    //dataDescription->AddInput("input");
    dataDescription->AddInput("input");
    //dataDescription->AddInput("input");
    dataDescription->SetTimeData(time, timeStep);
    if(lastTimeStep == true)
      {
      // assume that we want to all the pipelines to execute if it
      // is the last time step.
      dataDescription->ForceOutputOn();
      }
    if(Processor->RequestDataDescription(dataDescription.GetPointer()) != 0)
      {
      BuildVTKDataStructures(grid, attributes);
      dataDescription->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
      int wholeExtent[6];
      for(int i=0;i<6;i++)
        {
        wholeExtent[i] = grid.GetNumPoints()[i];
        }

      dataDescription->GetInputDescriptionByName("input")->SetWholeExtent(wholeExtent);

      Processor->CoProcess(dataDescription.GetPointer());
      }
  }
} // end of Catalyst namespace
