#!/bin/bash --login

surf=$1

cpu_id=${SLURM_PROCID}
energy=1e20
theta=10
ndivs=32
angle=`echo ${SLURM_ARRAY_TASK_ID} \* 10 + 5 | bc`
results_dir=diff_surfaces/Incoming_${angle}/Surface_${surf}/Results
input_dir=diff_surfaces/Incoming_${angle}/Surface_${surf}/InputFiles/rough_${cpu_id}.in
output_dir=diff_surfaces/Incoming_${angle}/Surface_${surf}/OutputFiles/rough_${cpu_id}.out

cp diff_surfaces/Incoming_${angle}/Surface_${surf}/rough_template.in ${input_dir}

sed -i "s|TTTT|..\/Cascades\/E${energy}_T${theta}_N${ndivs}_${cpu_id}.dat|g" ${input_dir}

sed -i "s|RRRP|1 "${results_dir}"\/E${energy}_T${theta}_N${ndivs}_${cpu_id}.phi|g" ${input_dir}

sed -i "s|RRRT|1 "${results_dir}"\/E${energy}_T${theta}_N${ndivs}_${cpu_id}.theta|g" ${input_dir}

#modify input files

./rough.exe < ${input_dir} > ${output_dir}
