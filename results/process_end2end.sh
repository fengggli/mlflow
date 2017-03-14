#!/bin/bash
if [ "$#" -ne 1 ]; then
   echo "specify jobdir"
   exit
fi

JOBDIR=$1
tmp_dir=./tmp
mkdir -pv ${tmp_dir}

less $JOBDIR/myicofoam.log.1.00 |grep 'comp sim'|awk '{print $NF}' > ${tmp_dir}/sim_comp_sim.tmp
less $JOBDIR/myicofoam.log.1.00 |grep 'all latency'|awk '{print $NF}' > ${tmp_dir}/sim_all_latency.tmp

less $JOBDIR/consumer.log.1.000 |grep 'comp divs'|awk '{print $NF}' > ${tmp_dir}/consumer_comp_divs.tmp
less $JOBDIR/consumer.log.1.000 |grep 'comp assign'|awk '{print $NF}' > ${tmp_dir}/consumer_comp_assign.tmp

less $JOBDIR/catalyst.log.1.0 |grep 'comp cat'|awk '{print $NF}' > ${tmp_dir}/catalyst_comp_cat.tmp
less $JOBDIR/catalyst.log.1.0 |grep 'all latency'|awk '{print $NF}' > ${tmp_dir}/catalyst_all_latency.tmp
