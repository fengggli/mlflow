cmake_minimum_required(VERSION 2.8 FATAL_ERROR)                                    

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -dynamic -std=c++11" CACHE STRING "cppflags")
SET( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -dynamic " CACHE STRING "cflags")

project("myIcofoam")                                                              
add_definitions ("-Wall -Wextra")                                                  
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")                                                                                                                                                                                         
#set(CMAKE_PREFIX_PATH "/N/u/lifen/BigRed2/OpenFOAM/OpenFOAM-v1612+/platforms/linux64GccDPInt32Opt/lib")
                                                                                   
#set(ParaView_DIR "/N/soft/cle4/paraview/ParaView-v4.1.0/Build-mpi" CACHE PATH "The path to paraview.")                                
set(ParaView_DIR "~/envs/Superbuild/build_mpi_static_syspython_catalyst/install/lib/cmake/paraview-5.2/" CACHE PATH "The path to paraview.")
set(DS_DIR "/N/u/lifen/BigRed2/envs/Dataspacesroot" CACHE PATH "The path to dataspaces.")


add_executable(myIcoFoam myIcoFoam.C)
set(OpenFoamRoot "/N/u/lifen/BigRed2/OpenFOAM/OpenFOAM-v1612+" CACHE PATH "The path to Openfoam.")                                
set(OpenFoamLib "${OpenFoamRoot}/platforms/linux64GccDPInt32Opt/lib" CACHE PATH "The path to Openfoam lib")                                

include_directories("${OpenFoamRoot}/src/finiteVolume/lnInclude")
include_directories("${OpenFoamRoot}/src/meshTools/lnInclude")
include_directories("${OpenFoamRoot}/src/OpenFOAM/lnInclude")
include_directories("${OpenFoamRoot}/src/OSspecific/POSIX/lnInclude")

include_directories()

include_directories(${DS_DIR}/include ${CMAKE_SOURCE_DIR}/cluster)                                
find_library(DSPACES dspaces ${DS_DIR}/lib)
find_library(DSCOMMON dscommon ${DS_DIR}/lib)
find_library(DART dart ${DS_DIR}/lib)

find_library(FV finiteVolume  ${OpenFoamLib} REQUIRED)
find_library(MESH meshTools  ${OpenFoamLib} REQUIRED)

target_link_libraries(myIcoFoam ${DSPACES} ${DSCOMMON} ${DART} ${FV} ${MESH})
