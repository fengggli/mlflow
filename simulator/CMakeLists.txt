cmake_minimum_required(VERSION 2.8)
set(CMAKE_EXE_LINKER_FLAGS "-dynamic ${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "link flags")
set(USE_CATALYST ON CACHE BOOL "Link the simulator with Catalyst")
#set(CMAKE_LIBRARY_PATH "${DS_DIR}/lib" CACHE STRING "library path")



# project(Simulator)

add_executable(sim_gen sim_gen.c)                                                                                                                     
include_directories("${CMAKE_SOURCE_DIR}/src" ${DS_DIR}/include)                                                                                                       
find_library(DSPACES dspaces ${DS_DIR}/lib)
find_library(DSCOMMON dscommon ${DS_DIR}/lib)
find_library(DART dart ${DS_DIR}/lib)

target_link_libraries(sim_gen ${DSPACES} ${DSCOMMON} ${DART} UTILITY)



