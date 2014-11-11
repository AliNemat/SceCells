#!/bin/bash
#SBATCH -o myjob.%j.%N.out 
##SBATCH -D 
#SBATCH -J tutorial_1 
#SBATCH --ntasks=1 
#SBATCH --mail-type=end 
#SBATCH --mail-user=@nd.edu 
#SBATCH --time=00:01:00 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --share
#SBATCH --gres=gpu:1

sleep 10
srun -l ./bin/runDiskSimulation
