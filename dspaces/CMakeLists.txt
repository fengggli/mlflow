cmake_minimum_required(VERSION 2.8)
#set(CMAKE_LIBRARY_PATH "${DS_DIR}/lib" CACHE STRING "library path")


add_executable(analysis analysis.c)                                                                                                                     
add_executable(get_regions get_regions.c )                                                                                                                     

add_executable(consumer consumer.c )                                                                                                                     
add_executable(put_regions put_regions.c run_with_dspaces.c)
#add_executable(catalyst_2 catalyst.c)


set(DS_DIR "/N/u/lifen/BigRed2/Dataspacesroot" CACHE PATH "The path to dataspaces.")                                                                  
#include_directories(${DS_DIR}/include ${CMAKE_SOURCE_DIR}/cluster)                                                                                                       
find_library(DSPACES dspaces ${DS_DIR}/lib)
find_library(DSCOMMON dscommon ${DS_DIR}/lib)
find_library(DART dart ${DS_DIR}/lib)

#Dynamic/Shared Libs
#Static start
#set_target_properties(put_regions PROPERTIES LINK_SEARCH_START_STATIC 1)
#set_target_properties(put_regions PROPERTIES LINK_SEARCH_END_STATIC 1)
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#Static Libs

target_link_libraries(put_regions ${DSPACES} ${DSCOMMON} ${DART} UTILITY )
target_link_libraries(analysis ${DSPACES} ${DSCOMMON} ${DART} UTILITY CLUSTER)
target_link_libraries(get_regions ${DSPACES} ${DSCOMMON} ${DART} UTILITY )

#target_link_libraries(catalyst_2 ${DSPACES} ${DSCOMMON} ${DART} UTILITY )

target_link_libraries(consumer ${DSPACES} ${DSCOMMON} ${DART} UTILITY )


