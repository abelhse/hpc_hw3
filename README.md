```bash
module load openmpi/4.0.1
mpicxx -o main_mpi main_mpi.cpp
# srun -N7 --ntasks-per-node=1 main_mpi T Nt Nx
srun -N 5 --ntasks-per-node=1 main_mpi 0.1 501 51
```