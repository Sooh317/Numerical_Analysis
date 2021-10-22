#!/bin/bash
#SBATCH --nodes 1 
#SBATCH --ntasks-per-node 1
#SBATCH -t 00:01:00
#SBATCH -o result.txt

./omp
./cuda
