from mpi4py import MPI
import numpy as np
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

gros_tableau = None
taille_gros_tableau = 1000000

comm.barrier()
t_deb = time.time()
if rank == 0 :
    gros_tableau = np.empty(taille_gros_tableau, dtype=np.float64)
    gros_tableau[0::size] = rank
    comm.send(gros_tableau, dest=1, tag=101)
    gros_tableau = comm.recv(source = size-1, tag=101)
else:
    gros_tableau = comm.recv(source = rank-1, tag=101)
    gros_tableau[rank::size] = rank
    comm.send(gros_tableau, dest=(rank+1)%size, tag=101)
t_fin = time.time()

if rank == 0 :
    print(f"Temps de circulation dans l'anneau avec sérialisation : {t_fin-t_deb} secondes")
    print(f"Vérification des premières valeurs du tableau : {gros_tableau[:10]}")

comm.barrier()
t_deb = time.time()
gros_tableau = np.empty(taille_gros_tableau, dtype=np.float64)
if rank == 0 :
    gros_tableau[0::size] = rank
    comm.Send(gros_tableau, dest=1, tag=101)
    comm.Recv(gros_tableau, source = size-1, tag=101)
else:
    comm.Recv(gros_tableau, source = rank-1, tag=101)
    gros_tableau[rank::size] = rank
    comm.Send(gros_tableau, dest=(rank+1)%size, tag=101)
t_fin = time.time()

if rank == 0 :
    print(f"Temps de circulation dans l'anneau sans sérialisation : {t_fin-t_deb} secondes")
    print(f"Vérification des premières valeurs du tableau : {gros_tableau[:10]}")
