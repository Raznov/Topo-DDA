#!/bin/bash
#SBATCH --job-name=GpuEvo
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10240MB
#SBATCH --time=00:01:00
#SBATCH --partition=commons
#SBATCH --mail-user=yz128@rice.com
#SBATCH --output="SlurmOut/job-%j.out"
#SBATCH --gres=gpu:1
echo "I Ran on:"

cat $SLURM_NODEFILE
#
cd $SLURM_SUBMIT_DIR

echo "I started on:"
date

./te

echo "I finished on:"
date

