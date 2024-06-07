#!/bin/bash
#SBATCH -J py_gyc
#SBATCH --partition gpu
#SBATCH --gpus t4:1
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -o stdout.%j
#SBATCH -e stderr.%j       

module load conda
conda activate uni_core
python 8_read_inference_results.py