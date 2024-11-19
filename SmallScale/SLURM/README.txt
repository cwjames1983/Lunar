The code here was developed by Lisa Holland (2024).

This directory contains code to generate many cascade templates (e.g. for a mass job submission).

REQUIREMENTS

Looks for surfaces in
../Surfaces/ (will automatically create surfaces in here)
../Cascades/ (Cascade will have form of E${energy}_T${theta}_N${ndivs}_${cpu_id}.dat)


##### setup.sh ####
High level program to submit many SLURM jobs, with each corresponding to a single surfave

####### run.sh ########

This is a SLURM submission script for Garrawarla. Run with:
sbatch run.sh [SURFACE_ID]

It takes as input a surface file id (e.g. 1--9)

It generates the directory structure:
    - diff_surfaces/Incoming_${angle}/Surface_${surf}
    - diff_surfaces/Incoming_${angle}/Surface_${surf}/InputFiles
    - diff_surfaces/Incoming_${angle}/Surface_${surf}/OutputFiles
    - diff_surfaces/Incoming_${angle}/Surface_${surf}/Results



It then generates a surface file according to "surf.in" (copied to surf_ID.in), and modifies this to specify the output files.

It also changes a line in rough_template (which is copied locally) to specify the new input surface

Finally, it calls "run_actually" which calls rough.exe


###### run_actually.sh ##########
Script to actually run once instance of rough.exe on Garrawarla (or equivalent SLURM managed system). Does this for a single track for a single surface on a single input angle

Needs to know CPU ID to know which track file to input

Identifies itself with SLURM_PROCID via angle=`echo ${SLURM_ARRAY_TASK_ID} \* 10 + 5 | bc`

Modifies input file "input_dir" with
    - output file
    - .phi output file
    - .theta output file

Calls rough.exe, directing to output file.

