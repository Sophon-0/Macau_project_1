#!/usr/bin/env zsh

#BSUB -J mi_job
#BSUB -P jrc_combine
##BSUB -B -N -u mi.yang0586@gmail.com

### Maximum runtime for the job
#BSUB -W 100:00

### Number of processors you want to use. Adjust to whatever you use in your script.
#BSUB -n 1

### Memory to be reserved on the node(s)
#BSUB -M 1000



cd /home/my871390/MI_YANG/RWTH_Aachen/MACAU_PROJECT/script

export PYTHONPATH=/home/my871390/.local/lib/python2.7/site-packages:$PYTHONPATH
export PATH=/home/my871390/.local/bin:$PATH
export LD_LIBRARY_PATH=/home/my871390/.local/lib/python2.7/site-packages:$LD_LIBRARY_PATH

module switch intel gcc/6.3
module load python

python macau_GDSC.py

## "$tissue" 
