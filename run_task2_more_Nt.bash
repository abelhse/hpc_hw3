#!/bin/bash


for nproc in 1 2 4 8 16 24; do
for Nx in 10000 25000 50000; do
for algo in 0 1 2; do
# srun -N nproc --ntasks-per-node=1 main_mpi T Nt Nx algo_num
sbatch --constraint="type_d" --ntasks=${nproc} --nodes=1 --wrap="srun main_mpi 0.001 501 ${Nx} ${algo}"
done;
done;
done;