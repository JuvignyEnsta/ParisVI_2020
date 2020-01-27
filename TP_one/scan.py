import numpy as np
import numpy.linalg as linalg

u1 = np.arange(0,10,0.01)
u2 = np.arange(0,10,0.01)

# Dependance avant : pas de soucis a
#                    vectoriser
for i in range(0,u1.shape[0]-1):
    u1[i] += u1[i+1]
# OK, vectorisation valide
u2[0:-1] += u2[1:]

# Difference entre la solution non vectorisee
# et la solution vectorisee( doit etre nulle).
diff = u2 - u1
errL2 = linalg.norm(diff)
print(f"Erreur L2 sur les deux calculs : {errL2}")


u1 = np.arange(0,10,0.01)
u2 = np.arange(0,10,0.01)

# Dependance arriere : Vectorisation
#                      fausse !
for i in range(1,u1.shape[0]):
    u1[i] += u1[i-1]
# vectorisation invalide !
u2[1:] += u2[0:-1]

# Difference entre la solution non vectorisee
# et la solution vectorisee (difference non
# nulle !)
diff = u2 - u1
errL2 = linalg.norm(diff)
print(f"Erreur L2 sur les deux calculs : {errL2}")

