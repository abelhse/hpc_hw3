```bash
module load openmpi/4.0.1
mpicxx -o main_mpi main_mpi.cpp
# srun -N nproc --ntasks-per-node=1 main_mpi T Nt Nx [isend_irecv_1 | isend_irecv | sendrecv]
srun -N 5 --ntasks-per-node=1 main_mpi 0.1 501 51 sendrecv
```