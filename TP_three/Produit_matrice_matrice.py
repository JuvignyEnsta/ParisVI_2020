import math
import numpy as np

def compute_tensors(dim):
    """
    Calcul deux paires de tenseurs qui serviront pour définir deux matrices de rang 1
    
    :param      dim:  La dimension des matrices (carrées)
    :type       dim:  Entier
    """
    u1 = np.cos(np.array([1.67*i*math.pi/dim for i in range(dim)]))
    u2 = np.sin(np.array([2.03*i*math.pi/dim for i in range(dim)]))
    v1 = np.cos(np.array([1.23*i*i*math.pi/(7.5*dim) for i in range(dim)]))
    v2 = np.sin(np.array([0.675*i/(3.1*dim) for i in range(dim)]))
    return u1, u2, v1, v2

def verif_product( uA, vA, uB, vB, C):
    """
    Vérification du résultat produit matrice matrice grâce à la formule :
    A = u_A.v_A^{T}
    B = u_B.v_B^{T}
    => C= A.B = u_A.<v_A|u_B>v_B^{T} où < | > est un produit scalaire

    :param      uA:   Le vecteur u à gauche du produit tensoriel définissant A
    :type       uA:   Vecteur numpy
    :param      vA:   Le vecteur v à droite du produit tensoriel définissant A
    :type       vA:   Vecteur numpy
    :param      uB:   Le vecteur u à gauche du produit tensoriel définissant B
    :type       uB:   Vecteur numpy
    :param      vB:   Le vecteur v à droite du produit tensoriel définissant B
    :type       vB:   Vecteur numpy
    :param      C:    La matrice résultante du produit matrice-matrice
    :type       C:    Un vecteut deux dimensions de numpy
    """
    vA_dot_uB = np.dot(vA, uB)
    diff_values = np.abs(C - np.tensordot(uA, vA_dot_uB * vB, axes = 0))
    max_error = np.argmax(diff_values)
    if diff_values.flat[max_error] > 1.E-10:
        i = max_error//dim
        j = max_error%dim
        val = uA[i] * vA_dot_uB * vB[j]
        print(f"Erreur numerique : valeur attendue pour C({i},{j}) -> {val}")
        print(f"mais la valeur calculee est : {C[i,j]}")
        raise ValueError("Erreur dans le calcul du produit matrice matrice")

import sys
import time
dim = 1024
if len(sys.argv) > 1 : dim = int(sys.argv[1])

uA, vA, uB, vB = compute_tensors(dim)
A = np.tensordot(uA,vA,axes=0)
B = np.tensordot(uB,vB,axes=0)

start = time.time()
C = A.dot(B)
end   = time.time()

verif_product(uA, vA, uB, vB, C)
print("Test passe")
elapse_time = max(end-start,1.E-14)
print(f"Temps CPU produit matrice-matrice : {elapse_time} secondes")
print(f"MFlops pour le produit matrice matrice : {2*(dim**3)/elapse_time/1000000} MFlops")
