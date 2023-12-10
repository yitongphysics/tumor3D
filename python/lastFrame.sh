#!/bin/bash
#SBATCH -J my_python_program
#SBATCH --mem=10G
#SBATCH -t 4:00:00
#SBATCH -p scavenge
#SBATCH -o lastFrame.out

module load miniconda
conda activate py3_env
python lastFrame.py
