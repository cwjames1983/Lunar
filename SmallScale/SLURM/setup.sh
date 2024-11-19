#!/bin/bash

for surf in $(seq 0 8)
do

	sbatch run.sh ${surf}
	
done
