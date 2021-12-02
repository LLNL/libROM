#!/bin/bash
#SBATCH -N 1
#SBATCH -J pore_collapse_dmd
#SBATCH -t 1:00:00
#SBATCH -p pdebug
#SBATCH -o pore_collapse_dmd.log
#SBATCH --open-mode truncate

DATA_DIR="/usr/workspace/nlrom/poreCollapse"

rm -rf pore_collapse_list
mkdir pore_collapse_list

for i in $(seq 11 1 20)
do
	ls $DATA_DIR/${i}gpa/run_036.* -I run_036.00001 > pore_collapse_list/${i}gpa
        echo ${i}gpa >> pore_collapse_list/training_gpa
done

cp pore_collapse_list/training_gpa pore_collapse_list/testing_gpa
srun -n 4 ../../build/examples/dmd/pore_collapse
