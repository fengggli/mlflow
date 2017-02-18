cmake_minimum_required(VERSION 2.8 FATAL_ERROR)                                    


# this is openfoam cavaty generator
add_executable(myIcoFoam myIcoFoam.C)
set(OpenFoamRoot "/N/u/lifen/Karst/OpenFOAM/OpenFOAM-v1612+" CACHE PATH "The path to Openfoam.")                                
set(OpenFoamLib "${OpenFoamRoot}/platforms/linux64GccDPInt32Opt/lib" CACHE PATH "The path to Openfoam lib")                                

include_directories("${OpenFoamRoot}/src/finiteVolume/lnInclude")
include_directories("${OpenFoamRoot}/src/meshTools/lnInclude")
include_directories("${OpenFoamRoot}/src/OpenFOAM/lnInclude")
include_directories("${OpenFoamRoot}/src/OSspecific/POSIX/lnInclude")

find_library(FV finiteVolume  ${OpenFoamLib} REQUIRED)
find_library(MESH meshTools  ${OpenFoamLib} REQUIRED)

target_link_libraries(myIcoFoam ds_adaptor ${FV} ${MESH})
