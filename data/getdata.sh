#!/bin/bash

# how many points in each side
POINTS_SIDE=201
for timestep in {0..10}
do
    file_name=isotropic_${POINTS_SIDE}_${POINTS_SIDE}_1_t_${timestep}.h5
    vtk_name=isotropic_${POINTS_SIDE}_${POINTS_SIDE}_1_t_${timestep}.vti
    if [ ! -f $file_name ] ; then
        wget http://dsp033.pha.jhu.edu/jhtdb/getcutout/com.gmail.lf921227-069b89fb/isotropic1024coarse/p,u/${timestep},1/0,1/0,${POINTS_SIDE}/0,${POINTS_SIDE}/hdf5 -O $file_name
        echo "download ${POINTS_SIDE}*${POINTS_SIDE}*1 h5 data cut in timestep ${timestep}"
    fi

    if [ ! -f $vtk_name ] ; then
        wget http://dsp033.pha.jhu.edu/jhtdb/getcutout/com.gmail.lf921227-069b89fb/isotropic1024coarse/p,u/${timestep},1/0,1/0,${POINTS_SIDE}/0,${POINTS_SIDE}/vtk -O $vtk_name
        echo "download ${POINTS_SIDE}*${POINTS_SIDE}*1 vtk data cut in timestep ${timestep}"
    fi
done

echo 'all data files are present in data/ folder'
