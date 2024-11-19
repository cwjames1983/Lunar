#!/bin/bash --login

#SBATCH --job-name=thetas_surfaces
#SBATCH --account=mwavcs
#SBATCH --partition=workq
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=task_%x.out
#SBATCH --array=0-8

module load fftw/3.3.8

surf=$1

angle=`echo ${SLURM_ARRAY_TASK_ID} \* 10 + 5 | bc`

if [ ! -d diff_surfaces/Incoming_${angle} ]
then

	mkdir diff_surfaces/Incoming_${angle}
	cp rough_template.in diff_surfaces/Incoming_${angle}/

	sed -i "s/QQQQ/0 ${angle} 0.0 0.0 10.0/g" diff_surfaces/Incoming_${angle}/rough_template.in

fi

#for surf in $(seq 0 8)
#do

	if [ ! -d diff_surfaces/Incoming_${angle}/Surface_${surf} ]
	then

		mkdir diff_surfaces/Incoming_${angle}/Surface_${surf}
		mkdir diff_surfaces/Incoming_${angle}/Surface_${surf}/InputFiles
		mkdir diff_surfaces/Incoming_${angle}/Surface_${surf}/OutputFiles
		mkdir diff_surfaces/Incoming_${angle}/Surface_${surf}/Results

	fi

	cp diff_surfaces/Incoming_${angle}/rough_template.in diff_surfaces/Incoming_${angle}/Surface_${surf}/
	cd ../Surfaces/

	cp surf.in surf_${surf}.in

	sed -i "s/SSSS/1 rough_50m_50m_0.05m_${surf}.dat/g" surf_${surf}.in
	./surf.exe < surf_${surf}.in

	cd ../FacetPropagation/

	sed -i "s|SSSS|2 ../Surfaces/rough_50m_50m_0.05m_${surf}.dat|g" diff_surfaces/Incoming_${angle}/Surface_${surf}/rough_template.in

	srun run_actually.sh ${surf}

#done

