#!/bin/bash
#SBATCH --job-name=PartMC_case
#SBATCH -n 101
#SBATCH -p sesempi 
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4000
# Email if failed run
#SBATCH --mail-type=FAIL
# Email when finished
#SBATCH --mail-type=END
# My email address
#SBATCH --mail-user=zzheng25@illinois.edu

# the parent_path contains cases/ folder
export SLURM_SUBMIT_DIR=/data/keeling/a/zzheng25/scenario_generator
# case folder name under cases/case_*
export case=case_6
# number of (scenarios+1)
export scenario_num_plus_1=101
# PartMC path
export PMC_PATH=/data/keeling/a/zzheng25/partmc
# Path for results
export WORK_DIR=/data/keeling/a/zzheng25/d/partmc-mam4-high-lat

## INSTRUCTION
# module load gnu/openmpi-3.0.0-gnu-6.1.0
# cd Scheduler
# mpif90 -o scheduler.x scheduler.F90
# mv scheduler.x ..
# chmod 744 partmc_submit.sh
# sed -i 's/\r//g' partmc_submit.sh 
# Within partmc_submit.sh:
# - define #SBATCH -n XXX (XXX should be the same as scenario_num_plus_1)
# - define the variables in export XXXX
# - define the user email
# sbatch partmc_submit.sh

######## Do not Change ########
# The job script can create its own job-ID-unique directory 
# to run within.  In that case you'll need to create and populate that 
# directory with executables and inputs
mkdir -p $WORK_DIR/$SLURM_JOB_ID
cd $WORK_DIR/$SLURM_JOB_ID
cp -r $PMC_PATH/build build
cp -r $PMC_PATH/src src

# Copy the scenario directory that holds all the inputs files
cp -r $SLURM_SUBMIT_DIR/cases/$case/scenarios .

# Copy things to run this job
# Need the scheduler and the joblist
cp $SLURM_SUBMIT_DIR/scheduler.x .
cp $SLURM_SUBMIT_DIR/cases/$case/joblist .

# Run the library. One core per job plus one for the master.
mpirun -np $scenario_num_plus_1 ./scheduler.x joblist /bin/bash -noexit -nostdout > log
