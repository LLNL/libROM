#!/bin/bash

#SBATCH -N 3
#SBATCH -J rebaseline-libROM
#SBATCH -t 0:02:00
#SBATCH -p pdebug
#SBATCH -A modlred
#SBATCH -o rebaseline-libROM.out

# NOTE(oxberry1@llnl.gov): The execution time value is set
# deliberately low so that a job clears the queue quickly,
# meaning it can be executed repeatedly without charging
# banks too much time. In testing, the optimized version
# takes ~6s to execute, so a debug version shouldn't take
# more than 1-2 minutes.

set -x # Echo every command used in the script for safety
set -e # Halt script if any command returns nonzero exit code

# Delete any baselines kicking around in the output directory.
REPO_ROOT="$(git rev-parse --show-toplevel)"
UNEVEN_TEST=${REPO_ROOT}/uneven_dist
SMOKE_TEST=${REPO_ROOT}/smoke_test
RANDOM_TEST=${REPO_ROOT}/random_test
BASELINES_DIR=${REPO_ROOT}/BASELINES
OUTPUT_DIR=${REPO_ROOT}/NEW-BASELINES
mkdir -p ${OUTPUT_DIR}
rm -f ${OUTPUT_DIR}/*.out ${OUTPUT_DIR}/*.diff

# NOTE(oxberry1@llnl.gov): This script intended to be run using an
# sxterm. Quartz and RZTopaz are each 36-core machines so this loop
# runs process counts that are valid for a single node of a 36-core
# machine.
#
# As noted above, the script halts if any command returns a nonzero
# exit code. However, diff will return an exit code of 1 if the
# files are different, and will return an exit code greater than
# 1 if there is an error. To suppress spurious halting errors,
# we add || [[ $? == 1 ]] to the end of the diff command, which
# checks to see if the exit code of diff is 1. If the exit code of
# diff is 1, that logical clause will be true, so the compound
# command will return an exit code of zero, and the script will
# continue. 
for num_procs in 1 2 3 4 5 6;
do
    srun -N1 -n${num_procs} ${UNEVEN_TEST} \
	 > ${OUTPUT_DIR}/uneven_${num_procs}proc.out
    diff ${BASELINES_DIR}/uneven_${num_procs}proc.out \
	 ${OUTPUT_DIR}/uneven_${num_procs}proc.out \
	 > ${OUTPUT_DIR}/uneven_${num_procs}proc.diff \
	 || [[ $? == 1 ]]
done

for num_procs in 1 2 3 6;
do
    srun -N1 -n${num_procs} ${SMOKE_TEST} \
	 > ${OUTPUT_DIR}/smoke_${num_procs}proc.out
    diff ${BASELINES_DIR}/smoke_${num_procs}proc.out \
	 ${OUTPUT_DIR}/smoke_${num_procs}proc.out \
	 > ${OUTPUT_DIR}/smoke_${num_procs}proc.diff \
	 || [[ $? == 1 ]]
done

for num_procs in 1 2 4 5 10 20 25;
do
    srun -N1 -n${num_procs} ${RANDOM_TEST} \
	 > ${OUTPUT_DIR}/random_${num_procs}proc.out
    diff ${BASELINES_DIR}/random_${num_procs}proc.out \
	 ${OUTPUT_DIR}/random_${num_procs}proc.out \
	 > ${OUTPUT_DIR}/random_${num_procs}proc.diff \
	 || [[ $? == 1 ]]
done

srun -N2 -n50 ${RANDOM_TEST} > ${OUTPUT_DIR}/random_50proc.out
diff ${BASELINES_DIR}/random_50proc.out ${OUTPUT_DIR}/random_50proc.out \
     > ${OUTPUT_DIR}/random_50proc.diff || [[ $? == 1 ]]
srun -N3 -n100 ${RANDOM_TEST} > ${OUTPUT_DIR}/random_100proc.out
diff ${BASELINES_DIR}/random_100proc.out ${OUTPUT_DIR}/random_100proc.out \
     > ${OUTPUT_DIR}/random_100proc.diff || [[ $? == 1 ]]

