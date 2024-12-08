#!/bin/bash

# srun -N nproc --ntasks-per-node=1 main_mpi T Nt Nx algo_num
srun -N 1 --ntasks-per-node=1 main_mpi 0.1 501 51 1
srun -N 5 --ntasks-per-node=1 main_mpi 0.1 501 51 1