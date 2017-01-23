cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_EXE_LINKER_FLAGS "-dynamic ${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "link flags")
project("DataBroker")
add_definitions ("-Wall -Wextra")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

set(ParaView_DIR "/N/soft/cle4/paraview/ParaView-v4.1.0/Build-mpi" CACHE PATH "The path to paraview.")
set(DS_DIR "/N/u/lifen/BigRed2/Dataspacesroot" CACHE PATH "The path to dataspaces.")                                                                  

include_directories("${CMAKE_SOURCE_DIR}/src")     

add_subdirectory(src)
add_subdirectory(simulator)
add_subdirectory(dspaces)
add_subdirectory(cluster)


