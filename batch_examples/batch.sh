#The script must include the specific command for job submission used in you batch system followed by $1.
#In the following, we provide examples for SLURM and Task Spooler.

#SLURM
#In this case here, nemo.sh is another script that contains the details of submission and takes as 
#argument the input file. This example assumes the the nemo.sh file is in the same directory as this script.

sbatch nemo.sh $1


#TASK SPOOLER
#If you use task spooler, all you need is the following line:
ts bash $1
