#!/bin/bash
# download file from JHTDB, both large and small dataset are used
LARGE_FILE=isotropic_255_255_128.h5
MID_FILE=isotropic_255_255_5.h5
SMALL_FILE=isotropic_2_2_2.h5
# first visualization at Sep 27(side length = 10 for each region)
FILE_201=isotropic_201_201_1.h5
# Oct 1, second visualization at S(side length = 30 for each region)
FILE_601=isotropic_601_601_1.h5

# for debugging purpose
FILE_51=isotropic_51_51_1.h5

if [ ! -f  $LARGE_FILE ];then
        echo 'download 25*255*128(large) with velocity and pressure'
        wget http://dsp033.pha.jhu.edu/jhtdb/getcutout/com.gmail.lf921227-069b89fb/isotropic1024coarse/p,u/0,1/0,256/0,256/0,128/hdf5 -O $LARGE_FILE
fi


if [ ! -f $MID_FILE ] ; then
        echo 'download 25*255*5(mid) with velocity and pressure'
        wget http://dsp033.pha.jhu.edu/jhtdb/getcutout/com.gmail.lf921227-069b89fb/isotropic1024coarse/p,u/0,1/0,256/0,256/0,5/hdf5/ -O $MID_FILE
fi

if [ ! -f $SMALL_FILE ] ; then
        echo 'download 2*2*2(small) with velocity and pressure'
        wget http://dsp033.pha.jhu.edu/jhtdb/getcutout/com.gmail.lf921227-069b89fb/isotropic1024coarse/p,u/0,1/0,2/0,2/0,2/hdf5/ -O $SMALL_FILE
fi

if [ ! -f $FILE_201 ] ; then
        echo 'download 201*201*1 data cut'
        wget http://dsp033.pha.jhu.edu/jhtdb/getcutout/com.gmail.lf921227-069b89fb/isotropic1024coarse/p,u/0,1/0,1/0,201/0,201/hdf5 -O $FILE_201 
fi

if [ ! -f $FILE_51 ] ; then
        echo 'download 51*51*1 data cut'
        wget http://dsp033.pha.jhu.edu/jhtdb/getcutout/com.gmail.lf921227-069b89fb/isotropic1024coarse/p,u/0,1/0,1/0,51/0,51/hdf5 -O $FILE_51
fi

if [ ! -f $FILE_601 ] ; then
        echo 'download 601*601*1 data cut'
        wget http://dsp033.pha.jhu.edu/jhtdb/getcutout/com.gmail.lf921227-069b89fb/isotropic1024coarse/p,u/0,1/0,1/0,601/0,601/hdf5 -O $FILE_601 
fi

# how many points in each side
POINTS_SIDE=201
for timestep in {1..10}
do
    file_name=isotropic_${POINTS_SIDE}_${POINTS_SIDE}_1_t_${timestep}.h5
    if [ ! -f $file_name ] ; then
        wget http://dsp033.pha.jhu.edu/jhtdb/getcutout/com.gmail.lf921227-069b89fb/isotropic1024coarse/p,u/${timestep},1/0,1/0,${POINTS_SIDE}/0,${POINTS_SIDE}/hdf5 -O $file_name
        echo "download ${POINTS_SIDE}*${POINTS_SIDE}*1 data cut in timestep ${timestep}"
    fi
done

echo 'all data files are present in data/ folder'
