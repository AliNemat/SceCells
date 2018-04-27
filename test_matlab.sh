#!/bin/csh

#$ -M anematba@nd.edu	 # Email address for job notification
#$ -m  abe		 # Send mail when job begins, ends and aborts
#$ -q  *@@acms_gpu 	 # Specify queue
#$ -l gpu_card=1
#s -pe smp 4 
#$ -N  run_April19thM	 # Specify job name


module load matlab
matlab -singleCompThread -nodisplay -nosplash < main.m
