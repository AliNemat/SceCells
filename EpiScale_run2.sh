#!/bin/csh

#$ -M anematba@nd.edu	 # Email address for job notification
#$ -m  abe		 # Send mail when job begins, ends and aborts
#$ -q  *@@acms_gpu 	 # Specify queue
#$ -l gpu_card=1 
#s -pe smp 4 
#$ -N  run_May12	 # Specify job name


module load matlab
module load gcc/4.9.2
module load cuda/7.0
module load bertini
setenv LD_LIBRARY_PATH /afs/crc.nd.edu/x86_64_linux/b/boost/1.58/gcc/4.9.2/build/lib:$LD_LIBRARY_PATH 
echo -n "It is currently: ";date
echo -n "I am logged on as ";who am i
echo -n "This computer is called ";hostname
echo -n "I am currently in the directory ";pwd
#setenv PATH /afs/crc.nd.edu/user/a/anematba/Public/2015/Oct/11th/SceCells/bin:$PATH
./bin/runDiscSimulation_M

