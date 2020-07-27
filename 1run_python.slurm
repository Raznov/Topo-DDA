#!/bin/bash
#SBATCH --job-name=python
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=10240MB
#SBATCH --time=23:59:00
#SBATCH --partition=nordland
#SBATCH --account=nordland
#SBATCH --mail-type=ALL
#SBATCH --mail-user=slurmjob@outlook.com
#SBATCH --output="SlurmOut/job-%j.out"
echo "I Ran on:"

cat $SLURM_NODEFILE
#
cd $SLURM_SUBMIT_DIR

echo "I started on:"
date

~/.conda/envs/major/bin/python plot.py 828212-41416-25-508-01-PointE-Adam

echo "I finished on:"
date

