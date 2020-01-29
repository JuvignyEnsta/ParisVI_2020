from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.rank

if rank == 0 :
    u = np.ones(100, dtype=np.float64)
    comm.Send(u, dest=1, tag=101)
else:
    u = np.zeros(100, dtype=np.float64)
    comm.Recv(u, source=0, tag=101)

print(f"u = {u}")
