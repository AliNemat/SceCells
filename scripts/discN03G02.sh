#!/bin/bash
#SBATCH -o DiscSimu.%j.%N.out 
##SBATCH -D 
#SBATCH -J DiscSimu_N03_2 
#SBATCH --ntasks=1 
#SBATCH --mail-type=end 
#SBATCH --mail-user=wsun2@nd.edu 
#SBATCH --time=99:00:00 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --share
#SBATCH --gres=gpu:1
#SBATCH --nodelist=gpu03
srun --gres=gpu:1 ./bin/runDiskSimulation -slurm N03_2
