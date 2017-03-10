#!/bin/bash
if [ "$#" -ne 1 ]; then
   echo "specify jobdir"
   exit
fi

JOBDIR=$1
tmp_dir=./tmp
mkdir -pv ${tmp_dir}

less $JOBDIR/myicofoam.log.1.00 |grep 'comm raw'|awk '{print $NF}' > ${tmp_dir}/sim_comm_raw.tmp
less $JOBDIR/myicofoam.log.1.00 |grep 'comm sample'|awk '{print $NF}' > ${tmp_dir}/sim_comm_sample.tmp
less $JOBDIR/myicofoam.log.1.00 |grep 'comp sim'|awk '{print $NF}' > ${tmp_dir}/sim_comp_sim.tmp

less $JOBDIR/consumer.log.1.00 |grep 'comm raw'|awk '{print $NF}' > ${tmp_dir}/consumer_comm_raw.tmp
less $JOBDIR/consumer.log.1.00 |grep 'comm sample'|awk '{print $NF}' > ${tmp_dir}/consumer_comm_sample.tmp
less $JOBDIR/consumer.log.1.00 |grep 'comm divs'|awk '{print $NF}' > ${tmp_dir}/consumer_comm_divs.tmp
less $JOBDIR/consumer.log.1.00 |grep 'comm medoids'|awk '{print $NF}' > ${tmp_dir}/consumer_comm_medoids.tmp
less $JOBDIR/consumer.log.1.00 |grep 'comm cluster'|awk '{print $NF}' > ${tmp_dir}/consumer_comm_cluster.tmp
less $JOBDIR/consumer.log.1.00 |grep 'comp divs'|awk '{print $NF}' > ${tmp_dir}/consumer_comp_divs.tmp
less $JOBDIR/consumer.log.1.00 |grep 'comp assign'|awk '{print $NF}' > ${tmp_dir}/consumer_comp_assign.tmp

less $JOBDIR/analysis.log.1.0 |grep 'comm divs'|awk '{print $NF}' > ${tmp_dir}/analysis_comm_divs.tmp
less $JOBDIR/analysis.log.1.0 |grep 'comm medoids'|awk '{print $NF}' > ${tmp_dir}/analysis_comm_medoids.tmp
less $JOBDIR/analysis.log.1.0 |grep 'comp medoids'|awk '{print $NF}' > ${tmp_dir}/analysis_comp_medoids.tmp

less $JOBDIR/catalyst.log.1.0 |grep 'comm raw'|awk '{print $NF}' > ${tmp_dir}/catalyst_comm_raw.tmp
less $JOBDIR/catalyst.log.1.0 |grep 'comm cluster'|awk '{print $NF}' > ${tmp_dir}/catalyst_comm_cluster.tmp
less $JOBDIR/catalyst.log.1.0 |grep 'comp cat'|awk '{print $NF}' > ${tmp_dir}/catalyst_comp_cat.tmp
