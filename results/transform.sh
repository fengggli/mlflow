#/bin/bash

for timestep in {0..99}
do
    INPUT=clusterids/201_k_5_t_${timestep}.txt
    OUTPUT=cate_id/201_k_5_t_${timestep}.txt
    python transform.py $INPUT > $OUTPUT
done




