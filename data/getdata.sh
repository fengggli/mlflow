#!/bin/bash
# download file from JHTDB, both large and small dataset are used
LARGE_FILE=isotropic_255_255_128.h5
MID_FILE=isotropic_255_255_5.h5
SMALL_FILE=isotropic_2_2_2.h5
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

echo 'all data files are present in data/ folder'
