cmake_minimum_required(VERSION 2.8)
set(USE_CATALYST ON CACHE BOOL "Link the simulator with Catalyst")
#set(CMAKE_LIBRARY_PATH "${DS_DIR}/lib" CACHE STRING "library path")
#set(CMAKE_FIND_LIBRARY_SUFFIXES .a)


# a toy simulator
#add_executabl(sim_gen sim_gen.c)
#target_link_libraries(sim_gen ds_adaptor)


# set catalyst environment
set(USE_CATALYST ON CACHE BOOL "Link the simulator with Catalyst")
#set(ParaView_DIR "/N/soft/cle4/paraview/ParaView-v4.1.0/Build-mpi" CACHE PATH "The path to paraview.")
if(USE_CATALYST)
  find_package(ParaView 4.1 REQUIRED COMPONENTS vtkPVPythonCatalyst)
  include("${PARAVIEW_USE_FILE}")
  set(Adaptor_SRCS
    FEAdaptor.cxx
    )
  add_library(CxxImageDataExampleAdaptor ${Adaptor_SRCS})
  target_link_libraries(CxxImageDataExampleAdaptor vtkPVPythonCatalyst)
  #target_link_libraries(CxxImageDataExampleAdaptor vtkCoProcessorImplementation)
  add_definitions("-DUSE_CATALYST")
  if(NOT PARAVIEW_USE_MPI)
    message(SEND_ERROR "ParaView must be built with MPI enabled")
  endif()
else()
  find_package(MPI REQUIRED)
  #include_directories(${MPI_CXX_INCLUDE_PATH})
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif()


# catalyst is to read from dataspaces and connect to paraview server
add_executable(catalyst FEDriver.cxx FEDataStructures.cxx)
if(USE_CATALYST)

  #set_target_properties(catalyst PROPERTIES LINK_FLAGS "-dynamic")
  target_link_libraries(catalyst LINK_PRIVATE ds_adaptor CxxImageDataExampleAdaptor)

  include(vtkModuleMacros)
  include(vtkMPI)
  vtk_mpi_link(catalyst)
else()
  target_link_libraries(catalyst LINK_PRIVATE ${MPI_LIBRARIES})
endif()

