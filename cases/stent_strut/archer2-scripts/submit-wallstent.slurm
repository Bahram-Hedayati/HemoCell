#!/bin/bash --login
#SBATCH --job-name=stent_strut_wall
#SBATCH --time 0-48:00:00 #runtime
#SBATCH --nodes=18
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --account=e723-hemocell
#SBATCH --partition=standard
#SBATCH --qos=long

output_dir=/work/e723/e723/cspieke/stent_strut-crossing/wallstent
echo $output_dir

files=(crossing-wallstent.stl config_wallstent.xml RBC.xml PLT.xml cell_pos_files/wallstent/RBC.pos cell_pos_files/wallstent/PLT.pos stent_strut)

for file in "${files[@]}"; do
        echo $file
        cp $file $output_dir/
done

cd $output_dir

# srun launches the parallel program based on the SBATCH options
srun ./stent_strut config_wallstent.xml
