# Importation du module MPI
from mpi4py import MPI
# Pour mesurer le temps passer dans une partie du code
import time
# Duplication du communicateur :
comm = MPI.COMM_WORLD.Dup()
# Interrogation du contexte MPI :
rank = comm.rank
size = comm.size
# Ouverture d'un fichier nom unique mode simple
fileName = f"sortie{rank:03d}.txt"
print(f"filename 1 : {fileName}")
## Ouverture d'un fichier nom unique mode compliqué
###################################################
#nb_characters_for_rank = len(f"{size-1}")
#format_int = f":0{nb_characters_for_rank}d"
#fileName   = ("sortie{" + format_int + "}.txt").format(rank)
#print(f"filename 2 : {fileName}")
#
file = open(fileName, mode="w")
file.write(f"Rang du processus : {rank}\n")
file.write(f"Nombre de processus : {size}\n")

file.close()
