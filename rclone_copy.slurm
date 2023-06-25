#!/bin/bash

#SBATCH --job-name=rclone_copy
#SBATCH --mincpus=12
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64gb
#SBATCH --time=50:00:00
#SBATCH --partition=standard
#SBATCH --account=guangyao
###SBATCH --array 0-7
#SBATCH --output=%x.out

base_dir="/xdisk/guangyao/michellewei/PMC_OA_BIOGRID"

for i in {0..24};
do
    rclone copy ${base_dir}/ev_${i} aws:pmcbiogrid --progress -vv --transfers=100
done
