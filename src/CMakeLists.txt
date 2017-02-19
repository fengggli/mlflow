set(H5_DIR "/opt/cray/hdf5/1.8.16/CRAY/8.3" CACHE PATH "The path to dataspaces.")                                                                  

include_directories("${H5_DIR}/include")



add_library(UTILITY common_utility.c divide.c get_divs.c read_file.c region_def.c)

#find_library(H5 hdf5 ${H5_DIR}/lib/ static)

target_link_libraries(UTILITY hdf5)



