cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_EXE_LINKER_FLAGS "-dynamic ${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "link flags")
project("DataBroker")

add_subdirectory(src)
add_subdirectory(simulator)


