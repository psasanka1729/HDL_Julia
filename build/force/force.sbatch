#!/bin/bash
# Job name:
#SBATCH --job-name=force
#
# Partition:
#SBATCH --partition=CMT
#
# Processors:
#SBATCH --ntasks=32
#SBATCH --nodes=2
#
# Wall clock limit:
#SBATCH --time=7-00:00:00
#
# Email at end:
#SBATCH --mail-user=sxd190113@utdallas.edu
#SBATCH --mail-type=END
#
# Specify output files:
#SBATCH -o force.o%j
#------------------------------------------------------

echo "Starting job force"

module load launcher

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=force.task
export LAUNCHER_SCHED=interleaved

$LAUNCHER_DIR/paramrun

echo "Job force complete"
