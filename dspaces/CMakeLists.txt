cmake_minimum_required(VERSION 2.8)


add_executable(analysis analysis.c)                                                                                                                     
target_link_libraries(analysis ds_adaptor CLUSTER)

add_executable(get_regions get_regions.c )                                                                                                                     
target_link_libraries(get_regions ds_adaptor )

add_executable(consumer consumer.c sample.c)                                                                                                                     
target_link_libraries(consumer ds_adaptor )

add_executable(put_regions put_regions.c run_with_dspaces.c)
target_link_libraries(put_regions ds_adaptor )
#add_executable(catalyst_2 catalyst.c)

#Dynamic/Shared Libs
#Static start
#set_target_properties(put_regions PROPERTIES LINK_SEARCH_START_STATIC 1)
#set_target_properties(put_regions PROPERTIES LINK_SEARCH_END_STATIC 1)
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#Static Libs


#target_link_libraries(catalyst_2 ${DSPACES} ${DSCOMMON} ${DART} UTILITY )



