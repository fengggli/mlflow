cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# c and C++ flags
set(CMAKE_C_COMPILER "/N/soft/rhel6/openmpi/gnu/1.8.4/bin/mpicc" CACHE STRING "C compiler")
set(CMAKE_C_FLAGS "-std=gnu99 ${ADD_FLAGS}" CACHE STRING "cflags")
set(CMAKE_CXX_COMPILER "/N/soft/rhel6/gcc/4.9.4/bin/g++" CACHE STRING "C compiler")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${ADD_FLAGS} -std=c++11 -DOPENFOAM_PLUS=1612 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -DOMPI_SKIP_MPICXX -fPIC" CACHE STRING "cppflags")

project("DataBroker")
add_definitions ("-Wall -Wextra")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

#set(ParaView_DIR "/N/soft/cle4/paraview/ParaView-v4.1.0/Build-mpi" CACHE PATH "The path to paraview.")
set(ParaView_DIR "/N/soft/rhel6/paraview/ParaView-v4.3.1-source/Build" CACHE PATH "The path to paraview.")
#set(ParaView_DIR "~/envs/Superbuild/build_mpi_static_syspython_catalyst/install/lib/cmake/paraview-5.2/" CACHE PATH "The path to paraview.")
set(DS_DIR "/N/u/lifen/Karst/Dataspacesroot" CACHE PATH "The path to dataspaces.")                                                                  

include_directories(${DS_DIR}/include ${CMAKE_SOURCE_DIR}/cluster ${CMAKE_SOURCE_DIR}/simulator ${CMAKE_SOURCE_DIR}/common)                                                                                                       

add_subdirectory(common)
add_subdirectory(simulator)
add_subdirectory(applications)
add_subdirectory(cluster)
add_subdirectory(catalyst)
