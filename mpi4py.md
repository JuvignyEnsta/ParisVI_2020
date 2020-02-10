# Documentation sur mpi4py

## Exécuter un script python avec MPI

En général, la plupart des programmes MPI s'exécutent avec la commande `mpiexec`. En pratique pour exécuter,
par exemple, un script python sur quatre processus, il faut lancer la commande :

    $mpiexec -n 4 python script.py

## Communicateurs

Les communicateurs sous python offrent les mêmes souplesses que ceux en C ou en Fortran.  Ils permettent une partition des processus, la définition de groupe, etc.

Par défaut, MPI propose un communicateur `MPI.COMM_WORLD`
englobant tous les processus lancés par la commande `mpiexec`. On construit généralement d'autres communicateurs
à partir de ce communicateur "global".

Remarque : Ecrire en Python la ligne suivante :
---------
    comm = MPI.COMM_WORLD
**ne copie pas** le communicateur global mais effectue une simple référence sur ce communicateur.

### Interrogation d'un communicateur à l'aide d'attributs dérivés

En python, afin de simplifier l'interface, l'interrogation d'un communicateur pour connaître le rang
ou le nombre de processus contenus dans le communicateur, peut se faire à l'aide d'un attribut dérivé,
c'est à dire une valeur vue comme attribut du communicateur *du point de vue utilisateur*  mais qui en réalité va appeler une fonction C donnant le rang ou le nombre de processus.

Exemple de code :
----------------
    from mpi4py import MPI

    comm = MPI.COMM_WORLD # référence au communicateur global
    numero_processus = comm.rank
    nombre_processus = comm.size

    print(f"Hello from {numero_processus}/{nombre_processus}")

### Interrogation d'un communicateur à l'aide de fonction à la "C"

Il est néanmoins possible d'appeler une fonction équivalente à la fonction C proposée par l'API MPI.
Cette fonction est en générale plus simple d'emploi que son homologue en C et appelle directement la
fonction C correspondante :

Ainsi, l'exemple donné ci-dessus peut aussi s'écrire "à la C" comme suit :
    
    from mpi4py import MPI

    comm = MPI.COMM_WORLD # référence au communicateur global
    numero_processus = comm.Get_rank() # Equivalent à comm.rank
    nombre_processus = comm.Get_size() # Equivalent à comm.size

    print(f"Hello from {numero_processus}/{nombre_processus}")

### Manipulation des communicateurs

Pour avoir une copie effective du communicateur ( et non un référence sur le communicateur en utilisant le
signe d'égalité ), on peut utiliser la méthode `Dup`:

    from mpi4py import MPI

    comm_ref = MPI.COMM_WORLD # On n'a pas copié le communicateur global, on fait une simple référence dessus
    comm_cpy = comm_ref.Dup() # Là, on a fait une vraie copie du communicateur dans un nouveau communicateur. 

Il est possible de définir des sous-communicateurs en utilisant la méthode `Split` de la même manière qu'en C :

    from mpi4py import MPI

    comm = MPI.COMM_WORLD # référence au communicateur global
    numero_processus = comm.Get_rank()

    # On fait une partition processus de rang pair/ processus de rang impair
    couleur = numero_processus % 2

    nouveau_comm = comm.Split(couleur, numero_processus)
    rank_oddeven = nouveau_comm.Get_rank()
    size_oddeven = nouveau_comm.Get_size()

## Communication point à point

Les fonctions de base permettant de transférer des données entre les différents processus, sous python,
sérialisent les objets à l'aide du module `pickle`. Cela permet une grande souplesse pour l'envoie des
données (on peut envoyer et recevoir n'importe quel objet sérialisable au sens Python) 
mais en contrepartie, il y aura un surcout supplémentaire dû à cette sérialisation.

### Envoie de données sérialisées

Voici un exemple commenté (remarque : le tag n'est pas obligatoire et sera mis à une valeur par défaut) :

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0 :
        dico = { 'alfred' : 42, 'victor' : 58, 'suzanne' : 37 }
        comm.send(dico, dest = 1, tag = 11)
    else:
        dico = comm.recv(source = 0, tag = 11)

### Envoie de données sérialisées non bloquante

Il est également possible d'envoyer des données sérialisées à l'aide d'opérations de communications non
bloquantes :

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        data = {'a': 7, 'b': 3.14}
        req = comm.isend(data, dest=1, tag=11)
        req.wait()
    elif rank == 1:
        req = comm.irecv(source=0, tag=11)
        data = req.wait()

## Envoie de données non sérialisées

Lorsqu'on doit uniquement envoyer des données basiques (entiers ou réels) de même type, afin d'optimiser le coût
en communication, il est préférable d'éviter une sérialisation et utiliser l'interface qui appelle directement
les fonctions de communication du C. Par exemple, si on doit envoyer des données de type entier ou réel stockées
dans un tableau numpy, la façon la plus rapide de la faire est la suivante :

    from mpi4py import MPI
    import numpy

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # passing MPI datatypes explicitly
    if rank == 0:
        data = numpy.arange(1000, dtype='i')
        comm.Send([data, MPI.INT], dest=1, tag=77)
    elif rank == 1:
        data = numpy.empty(1000, dtype='i')
        comm.Recv([data, MPI.INT], source=0, tag=77)

    # automatic MPI datatype discovery
    if rank == 0:
        data = numpy.arange(100, dtype=numpy.float64)
        comm.Send(data, dest=1, tag=13)
    elif rank == 1:
        data = numpy.empty(100, dtype=numpy.float64)
        comm.Recv(data, source=0, tag=13)

Noter deux façon d'appeler les routines, soit en précisant la nature des éléments contenus dans le tableau
numpy soit en laissant mpi4py "deviner" quels types de données sont stockées dans le tableau numpy (en appelant
des fonctions de numpy donnant le type de données stockées !).


## Communication collective

Pour faire une barrière de synchronisation :

    from mpi4py import MPI
    import time

    comm = MPI.COMM_WORLD

    # On fait une barrière de synchronisation pour mesurer le temps pris par les processus pour une
    # section parallèle.
    comm.Barrier()
    t1 = time.time()
    ... # Section parallèle dont on veut mesurer le temps
    t2 = time.time()
    print(f"Temps passé en parallèle dans la section mesurée : {t2-t1} secondes")

Exemple de diffusion d'un dictionnaire python :

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        data = {'key1' : [7, 2.72, 2+3j],
                'key2' : ( 'abc', 'xyz')}
    else:
        data = None
    data = comm.bcast(data, root=0)

Exemple de distribution d'un objet python :

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        data = [(i+1)**2 for i in range(size)]
    else:
        data = None
    data = comm.scatter(data, root=0)
    assert data == (rank+1)**2

Exemple de regroupement de données (gather) :

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    data = (rank+1)**2
    data = comm.gather(data, root=0)
    if rank == 0:
        for i in range(size):
            assert data[i] == (i+1)**2
    else:
        assert data is None

De même, on peut aussi optimiser ces communications collectives pour des données homogènes dans des tableaux numpy :

Diffusion d'un tableau numpy :

    from mpi4py import MPI
    import numpy as np

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        data = np.arange(100, dtype='i')
    else:
        data = np.empty(100, dtype='i')
    comm.Bcast(data, root=0)
    for i in range(100):
        assert data[i] == i

Distribution d'un tableau numpy :

    from mpi4py import MPI
    import numpy as np

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    sendbuf = None
    if rank == 0:
        sendbuf = np.empty([size, 100], dtype='i')
        sendbuf.T[:,:] = range(size)
    recvbuf = np.empty(100, dtype='i')
    comm.Scatter(sendbuf, recvbuf, root=0)
    assert np.allclose(recvbuf, rank)

Rassemblement de données dans un tableau numpy :

    from mpi4py import MPI
    import numpy as np

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    sendbuf = np.zeros(100, dtype='i') + rank
    recvbuf = None
    if rank == 0:
        recvbuf = np.empty([size, 100], dtype='i')
    comm.Gather(sendbuf, recvbuf, root=0)
    if rank == 0:
        for i in range(size):
            assert np.allclose(recvbuf[i,:], i)

Exemple de produit parallèle matrice--vecteur :

    from mpi4py import MPI
    import numpy

    def matvec(comm, A, x):
        m = A.shape[0] # local rows
        p = comm.Get_size()
        xg = numpy.zeros(m*p, dtype='d')
        comm.Allgather([x,  MPI.DOUBLE],
                       [xg, MPI.DOUBLE])
        y = numpy.dot(A, xg)
        return y



## MPI-IO

I/O collectives avec les tableaux numpy :

    from mpi4py import MPI
    import numpy as np

    amode = MPI.MODE_WRONLY|MPI.MODE_CREATE
    comm = MPI.COMM_WORLD
    fh = MPI.File.Open(comm, "./datafile.contig", amode)

    buffer = np.empty(10, dtype=np.int)
    buffer[:] = comm.Get_rank()

    offset = comm.Get_rank()*buffer.nbytes
    fh.Write_at_all(offset, buffer)

    fh.Close()

I/O collective non contigues avec des tableaux numpy et des datatypes :

    from mpi4py import MPI
    import numpy as np

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    amode = MPI.MODE_WRONLY|MPI.MODE_CREATE
    fh = MPI.File.Open(comm, "./datafile.noncontig", amode)

    item_count = 10

    buffer = np.empty(item_count, dtype='i')
    buffer[:] = rank

    filetype = MPI.INT.Create_vector(item_count, 1, size)
    filetype.Commit()

    displacement = MPI.INT.Get_size()*rank
    fh.Set_view(displacement, filetype=filetype)

    fh.Write_all(buffer)
    filetype.Free()
    fh.Close()

