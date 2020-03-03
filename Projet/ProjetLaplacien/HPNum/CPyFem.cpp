# include <Python.h>
# include <math.h>
# include <list>
# include <algorithm>
# include <iostream>
# include <numpy/arrayobject.h>

static char compSkelSparseMatrix_doc[] = "\
Usage : compute_skeleton_SparseMatrix( pt_dof2elts, dof2elts )\n\
Calcul le graphe de la matrice creuse issue des éléments finis.\n\
elt2dofs    : Tableau donnant pour chaque élément les degrés de liberté associés. \n\
pt_dof2elts : pointeurs dans le tableau dof2elts donnant pour chaque degré de     \n\
              liberté les éléments associés.\n\
dof2elts    : tableau donnant pour chaque degré de liberté les éléments associés.\n";
static PyObject* py_compSkelSparseMatrix( PyObject* self, PyObject* args )
{
  PyArrayObject *PyElt2Dofs, *PyBegDof2Elts, *PyDof2Elts;

  if ( !PyArg_ParseTuple( args, "O!(O!O!)",
			  &PyArray_Type, &PyElt2Dofs,
			  &PyArray_Type, &PyBegDof2Elts,
			  &PyArray_Type, &PyDof2Elts ) ) return NULL;
  const long* elt2dofs    = (const long*)PyArray_DATA(PyElt2Dofs);
  const long* begDof2Elts = (const long*)PyArray_DATA(PyBegDof2Elts);
  const long* dof2elts    = (const long*)PyArray_DATA(PyDof2Elts);

  long nbDofs = long(PyArray_DIM(PyBegDof2Elts,0))-1;

  npy_intp nbRowsP1, nbNZ;
  nbRowsP1 = npy_intp(nbDofs) + 1;
  PyArrayObject* pyBegRows = (PyArrayObject*) PyArray_SimpleNew(1,&nbRowsP1,NPY_LONG);
  long* begRows = (long*)PyArray_DATA(pyBegRows);

  begRows[0] = 0;
  for ( int idof = 0; idof < nbDofs; ++idof ) {
    std::list<long> lst_neighbours;
    for ( long ptIElt = begDof2Elts[idof]; ptIElt < begDof2Elts[idof+1]; ++ptIElt ) {
      long iElt = dof2elts[ptIElt];
      for ( long ptDof = 3*iElt; ptDof < 3*(iElt+1); ptDof++ ) {
	lst_neighbours.push_back(elt2dofs[ptDof]);
      }
    }
    lst_neighbours.sort();
    lst_neighbours.unique();
    begRows[idof+1] = begRows[idof] + lst_neighbours.size();
  }
  
  nbNZ = begRows[nbDofs];
  PyArrayObject* pyIndCols = (PyArrayObject*) PyArray_SimpleNew(1,&nbNZ,NPY_LONG);
  long* indCols = (long*)PyArray_DATA(pyIndCols);

  for ( long idof = 0; idof < nbDofs; ++idof ) {
    std::list<long> lst_neighbours;
    for ( long ptIElt = begDof2Elts[idof]; ptIElt < begDof2Elts[idof+1]; ++ptIElt ) {
      long iElt = dof2elts[ptIElt];
      for ( long ptDof = 3*iElt; ptDof < 3*(iElt+1); ptDof++ ) {
	lst_neighbours.push_back(elt2dofs[ptDof]);
      }
    }
    lst_neighbours.sort();
    lst_neighbours.unique();
    long ind = 0;
    for ( std::list<long>::iterator itL = lst_neighbours.begin(); itL != lst_neighbours.end(); ++itL ) {
      indCols[begRows[idof]+ind] = (*itL);
      ind += 1;
    }
  }

  return Py_BuildValue("NN", pyBegRows, pyIndCols);
}
// ------------------------------------------------------------------------
static char addElementaryMatrixToCSRMatrix_doc[] = 
"Usage : add_elemMat_csrMatrix( (begRows, indCols, coefs), (indRows, indCols, elemMat) )"
"Rajoute la matrice élémentaire définie par le tuple (indRows, indCols,elemMat)"
"à la matrice creuse stockée CSR définie par le tuple (begRows, indCols, coefs).";
static PyObject*
py_addelementmatrix_csrmatrix( PyObject* self, PyObject* args )
{
  // Tableaux pour la matrice creuse :
  PyArrayObject *pysm_BegRows, *pysm_IndCols, *pysm_Coefs;
  // Tableaux pour la matrice élémentaire :
  PyArrayObject *pyem_IndRows, *pyem_IndCols, *pyem_Coefs;
  PyArrayObject *py_rowMask = NULL, *py_colMask = NULL;
  //
  if ( !PyArg_ParseTuple( args, "(O!O!O!)(O!O!O!)|O!O!",
			  &PyArray_Type, &pysm_BegRows,
			  &PyArray_Type, &pysm_IndCols,
			  &PyArray_Type, &pysm_Coefs,
			  &PyArray_Type, &pyem_IndRows,
			  &PyArray_Type, &pyem_IndCols,
			  &PyArray_Type, &pyem_Coefs,
			  &PyArray_Type, &py_rowMask,
			  &PyArray_Type, &py_colMask ) ) return NULL;
  const long* sm_begrows = (const long*)PyArray_DATA(pysm_BegRows);
  const long* sm_indcols = (const long*)PyArray_DATA(pysm_IndCols);
  double   * sm_coefs   = (   double*)PyArray_DATA(pysm_Coefs  );
  const long* sm_rowMask = NULL;
  if (py_rowMask != NULL) 
    sm_rowMask = (const long*)PyArray_DATA(py_rowMask  );
  const long* sm_colMask = NULL;
  if (py_colMask != NULL) 
    sm_colMask = (const long*)PyArray_DATA(py_colMask  );
  
  long nRowsMatElem = long(PyArray_DIM(pyem_IndRows,0));
  long nColsMatElem = long(PyArray_DIM(pyem_IndCols,0));
  const long* em_indrows = (const long*)PyArray_DATA(pyem_IndRows);
  const long* em_indcols = (const long*)PyArray_DATA(pyem_IndCols);

  if (sm_rowMask == NULL) {
    for (long iRow = 0; iRow < nRowsMatElem; ++iRow ) {
      long indRow = em_indrows[iRow];
      for ( long jCol = 0; jCol < nColsMatElem; ++jCol ) {
	long indCol = em_indcols[jCol];
	for ( long ptCol = sm_begrows[indRow]; ptCol<sm_begrows[indRow+1];++ptCol ) {
	  if ( sm_indcols[ptCol] == indCol ) {
	    sm_coefs[ptCol] += *(double*)PyArray_GETPTR2(pyem_Coefs,iRow,jCol);
	    break;
	  }
	}
      }
    }
  } else {
    for (long iRow = 0; iRow < nRowsMatElem; ++iRow ) {
      long indRow = em_indrows[iRow];
      if (sm_rowMask[indRow] == 1) {
	for ( long jCol = 0; jCol < nColsMatElem; ++jCol ) {
	  long indCol = em_indcols[jCol];
	  for ( long ptCol = sm_begrows[indRow]; ptCol<sm_begrows[indRow+1];
		++ptCol ) {
	    if ( sm_indcols[ptCol] == indCol ) {
	      sm_coefs[ptCol] += *(double*)PyArray_GETPTR2(pyem_Coefs,iRow,jCol);
	      break;
	    }
	  }
	}
      }
    }
  }
  Py_INCREF(Py_None);
  return Py_None;
}
// ========================================================================
static PyMethodDef Py_Methods[] =
  {
    {"comp_skel_csr_mat", py_compSkelSparseMatrix, METH_VARARGS, compSkelSparseMatrix_doc},
    {"add_elt_mat_to_csr_mat", py_addelementmatrix_csrmatrix, 
     METH_VARARGS, addElementaryMatrixToCSRMatrix_doc},
    {NULL, NULL} /* Guards */
  };
// ========================================================================
#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef femmoduledef = {
    PyModuleDef_HEAD_INIT,
    "fem",                      // Nom du module
    "Finite elements method",   // Documentation du module,
    -1                      ,   // Taille du module ?
    Py_Methods              ,   // Fonctions du module
    NULL,                       // m_reload ?
    NULL,                       // m_traverse ?
    NULL,                       // m_clear
    NULL,                       // m_free
  };
PyMODINIT_FUNC PyInit_fem(void)
{
  PyObject* m; 
  m = PyModule_Create(&femmoduledef);
  if (m==NULL) return NULL;
  /*  important : initialize numpy to use in subroutines !!!! */
  import_array();
  return m;
}
#else  
static char fem_doc[] = "Finite Element Method";
PyMODINIT_FUNC initfem()
{
  PyObject* m; 
  m = Py_InitModule4("fem",Py_Methods, fem_doc,
		     (PyObject*)NULL,PYTHON_API_VERSION);
  if (m==NULL) return;
  /*  important : initialize numpy to use in subroutines !!!! */
  import_array();
}
#endif
