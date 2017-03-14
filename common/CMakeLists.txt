set(H5_DIR "/opt/cray/hdf5/1.8.16/CRAY/8.3" CACHE PATH "The path to dataspaces.")                                                                  

include_directories("${H5_DIR}/include")

find_library(DSPACES dspaces ${DS_DIR}/lib)
find_library(DSCOMMON dscommon ${DS_DIR}/lib)
find_library(DART dart ${DS_DIR}/lib)

# this is the libary to interact with dspaces
add_library(ds_adaptor ds_adaptor.c)
target_link_libraries(ds_adaptor ${DSPACES} ${DSCOMMON} ${DART} UTILITY)

add_library(UTILITY common_utility.c divide.c get_divs.c read_file.c region_def.c)

#find_library(H5 hdf5 ${H5_DIR}/lib/ static)

target_link_libraries(UTILITY hdf5)



