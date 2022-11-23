#!/bin/bash
# Job name:
#SBATCH --job-name=2_force
#
# Partition:
#SBATCH --partition=128s
#
# Processors:
#SBATCH --ntasks=32
#SBATCH --nodes=2
#
# Wall clock limit:
#SBATCH --time=96:00:00
#
# Email at end:
#SBATCH --mail-user=sxd190113@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o 2_force.o%j
#------------------------------------------------------

echo "Starting job 2_force"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=2_force.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job 2_force complete"
