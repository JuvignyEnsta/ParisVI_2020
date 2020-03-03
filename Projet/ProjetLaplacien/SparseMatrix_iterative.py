#!/bin/env python
# -*- coding: utf-8 -*-
import mesh
import fem
import laplacian
import splitter
from math import cos,sin,pi,sqrt
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import VisuSplitMesh as VSM
import VisuSolution as VS
from scipy.sparse import linalg as sp_linalg

def g(x,y) :
    return cos(2*pi*x)+sin(2*pi*y)

m = mesh.read("CarreMedium.msh")
coords    = m[0]
elt2verts = m[1]
nbVerts = coords.shape[0]
nbElts  = elt2verts.shape[0]
print('nbVerts : {}'.format(nbVerts))
print('nbElts  : {}'.format(nbElts))
begVert2Elts, vert2elts = mesh.compvert2elts(elt2verts)

begRows, indCols = fem.comp_skel_csr_mat(elt2verts, (begVert2Elts, vert2elts) )
nz = begRows[-1]
print("Number of non zero in sparse matrix : {}".format(nz))

spCoefs = np.zeros( (nz,), np.double)

for iElt in range(nbElts):
    iVertices = elt2verts[iElt,:]
    crd1 = coords[iVertices[0],:]
    crd2 = coords[iVertices[1],:]
    crd3 = coords[iVertices[2],:]
    matElem = laplacian.comp_eltmat((crd1[0],crd1[1]), (crd2[0],crd2[1]),
                                    (crd3[0],crd3[1]))
    fem.add_elt_mat_to_csr_mat((begRows,indCols,spCoefs),
                               (iVertices, iVertices, matElem))

# Assemblage second membre :
f = np.zeros(nbVerts, np.double)
for iVert in range(nbVerts):
    if ( coords[iVert,3] > 0 ) :
        f[iVert] += g(coords[iVert,0],coords[iVert,1])
b = np.zeros(nbVerts, np.double)
for i in range(nbVerts) :
    for ptR in range(begRows[i],begRows[i+1]):
        b[i] -= spCoefs[ptR]*f[indCols[ptR]]
        
# Il faut maintenant tenir compte des conditions limites :
for iVert in range(nbVerts):
    if coords[iVert,3] > 0: # C'est une condition limite !
        # Suppression de la ligne avec 1 sur la diagonale :
        for i in range(begRows[iVert],begRows[iVert+1]):
            if indCols[i] != iVert :
                spCoefs[i] = 0.
            else :
                spCoefs[i] = 1.
        # Suppression des coefficients se trouvant sur la colonne iVert :
        for iRow in range(nbVerts):
            if iRow != iVert :
                for ptCol in range(begRows[iRow],begRows[iRow+1]):
                    if indCols[ptCol] == iVert :
                        spCoefs[ptCol] = 0.
                        
        b[iVert] = f[iVert]
# On definit ensuite la matrice :
spMatrix = sparse.csc_matrix((spCoefs, indCols, begRows),
                             shape=(nbVerts,nbVerts))
#print("Matrice creuse : {}".format(spMatrix))

# Visualisation second membre :
VS.view( coords, elt2verts, b, title = "second membre" )

# Résolution du problème séquentiel avec gradient conjugue :
iteration = 1
def display_iter(x_k):
    global iteration
    print("iteration {:03d} : ||A.xk-b||/||b|| = {}".format(iteration, np.linalg.norm(spMatrix * x_k - b)/np.linalg.norm(b)))
    iteration += 1

# Résoud à l'aide du gradient conjugué.
sol,info = sp_linalg.cg(spMatrix, b, x0=None, tol=1.E-14, M=None,callback=display_iter)
#print(sol)
d = spMatrix.dot(sol) - b
print("||A.x-b||/||b|| = {}".format(sqrt(d.dot(d)/b.dot(b))))
# Visualisation de la solution :
VS.view( coords, elt2verts, sol, title = "Solution" )
    
nbDoms = 8
ndsDomains = splitter.node( nbDoms, coords )
print("Nombre de noeuds par domaines : ")
i = 0
for a in ndsDomains :
    print("Domaine {} : {} sommets\n".format(i,a.shape[0]))
    i += 1

etsDomains = splitter.element( nbDoms, (elt2verts, coords) )
print("Nombre d'elements par domaines : ")
i = 0
for a in etsDomains :
    print("Domaine {} : {} elements\n".format(i,a.shape[0]))
    i += 1

elt2doms = np.zeros((nbElts,), np.double)
ia = 0.
for a in etsDomains :
    for e in a :
        elt2doms[e] = ia
    ia += 1

# Calcul l'interface :
ie = 0
mask = np.array([-1,]*nbVerts, np.short)
for e in elt2verts :
    d = elt2doms[ie]
    if mask[e[0]] == -1 :
        mask[e[0]] = d
    elif mask[e[0]] != d :
        mask[e[0]] = -2
    if mask[e[1]] == -1 :
        mask[e[1]] = d
    elif mask[e[1]] != d :
        mask[e[1]] = -2
    if mask[e[2]] == -1 :
        mask[e[2]] = d
    elif mask[e[2]] != d :
        mask[e[2]] = -2
    ie += 1

nbInterf = 0
for m in mask :
    if m == -2 :
        nbInterf += 1

interfNodes = np.empty(nbInterf, np.long)
nbInterf = 0
for im in range(mask.shape[0]):
    if mask[im] == -2 :
        interfNodes[nbInterf] = im
        nbInterf += 1

VSM.view( coords, elt2verts, nbDoms, elt2doms, indInterfNodes = interfNodes, title='Partition par elements')

vert2dom = np.zeros((nbVerts,), np.double)
ia = 0.
for a in ndsDomains :
    for v in a :
        vert2dom[v] = ia
    ia += 1

    
mask = np.zeros((nbVerts,), np.short)
for e in elt2verts :
    d1 = vert2dom[e[0]]
    d2 = vert2dom[e[1]]
    d3 = vert2dom[e[2]]
    if (d1 != d2) or (d1 != d3) or (d2 != d3) :
        mask[e[0]] = 1
        mask[e[1]] = 1
        mask[e[2]] = 1

nbInterf = 0
for m in mask :
    if m == 1 :
        nbInterf += 1
interfNodes = np.empty(nbInterf, np.long)
nbInterf = 0
for im in range(mask.shape[0]):
    if mask[im] == 1 :
        interfNodes[nbInterf] = im
        nbInterf += 1

VSM.view( coords, elt2verts, nbDoms, vert2dom, indInterfNodes = interfNodes, title='Partition par sommets')
