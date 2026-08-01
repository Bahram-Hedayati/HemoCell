#!/bin/bash --login
#SBATCH --job-name=stent_strut_ref_ch
#SBATCH --time 0-48:00:00 #runtime
#SBATCH --nodes=50
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --account=e723-hemocell
#SBATCH --partition=standard
#SBATCH --qos=long

output_dir=/work/e723/e723/cspieke/stent_strut/reference
echo $output_dir

cd $output_dir

# srun launches the parallel program based on the SBATCH options
srun ./stent_strut tmp/checkpoint/checkpoint.xml
