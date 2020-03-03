# include <Python.h>
# include <math.h>
# include <list>
# include <iostream>
# include <algorithm>
# include <numpy/arrayobject.h>
# include <metis.h>

static char doc_splitNodeMesh[] =
  "Usage : metis.splitNodeMesh( nbDoms, (begRows, indCols) )"
  "where nbDoms is the number of subdomains"
  "and (begRows,indCols) are the graph of the matrix."
  "return set of nodes per domains\n";
static PyObject* py_splitNodeMesh( PyObject* self, PyObject* args )
{
  PyArrayObject *pyBegRows, *pyIndCols;
  int nbSubDomains;
  
  if ( !PyArg_ParseTuple( args, "i(O!O!)", &nbSubDomains,
			  &PyArray_Type, &pyBegRows,
			  &PyArray_Type, &pyIndCols ) ) return NULL;
  const long* begRows = (const long*)PyArray_DATA(pyBegRows);
  const long* indCols = (const long*)PyArray_DATA(pyIndCols);

  idx_t nbVerts = idx_t(PyArray_DIM(pyBegRows,0))-1;
  idx_t ncon    = 1;
  idx_t* xadj   = new idx_t[nbVerts+1];
  std::copy(begRows, begRows + nbVerts + 1, xadj );
  idx_t* adjncy = new idx_t[xadj[nbVerts]];
  std::copy(indCols, indCols + begRows[nbVerts], adjncy );
  idx_t nbDoms = nbSubDomains;
  idx_t objval;
  idx_t *part   = new idx_t[nbVerts];  
  int ok = METIS_PartGraphKway(&nbVerts, &ncon, xadj, adjncy, NULL, NULL, NULL,
			       &nbDoms, NULL, NULL, NULL, &objval, part );
  /* On compte maintenant le nombre de noeuds appartenant à chaque domaine :
   */
  long *nbNodesPerDomains = new long[nbSubDomains];
  std::fill_n(nbNodesPerDomains, nbSubDomains, 0);
  for ( int i = 0; i < nbVerts; ++i ) {
    nbNodesPerDomains[part[i]] += 1;
  }
  PyObject* pyLst = PyList_New(nbSubDomains);
  PyArrayObject* indVertSubDomains;
  long** ptIndVertices = new long*[nbSubDomains];
  for ( long i = 0; i < nbSubDomains; ++i ) {
    npy_intp nbVerts = nbNodesPerDomains[i];
    indVertSubDomains = (PyArrayObject*)PyArray_SimpleNew(1,&nbVerts,PyArray_LONG);
    PyList_SetItem(pyLst, i, (PyObject*)indVertSubDomains);
    ptIndVertices[i] = (long*)PyArray_DATA(indVertSubDomains);
  }
  std::fill_n(nbNodesPerDomains, nbSubDomains, 0);
  for ( long i = 0; i < nbVerts; ++i ) {
    ptIndVertices[part[i]][nbNodesPerDomains[part[i]]] = i;
    nbNodesPerDomains[part[i]] += 1;
  }
  delete [] ptIndVertices;
  delete [] nbNodesPerDomains; delete [] part;
  delete [] adjncy; delete [] xadj;
  return Py_BuildValue("N", pyLst);
}
// ------------------------------------------------------------------------
static char doc_splitEltMesh[] =
  "Usage : metis.splitEltMesh( nbDoms, nbVerts, elt2verts )"
  "where nbDoms is the number of subdomains, "
  "nbVerts the number of vertices in the mesh and "
  "elt2verts is the element to vertices connexion."
  "return a set of elements per domains\n";

static PyObject* py_splitEltMesh( PyObject* self, PyObject* args )
{
  int nbSubDoms, nbVerts;
  PyArrayObject *py_elt2verts;
  if (!PyArg_ParseTuple( args, "iiO!", &nbSubDoms, &nbVerts,
			 &PyArray_Type, &py_elt2verts ) )
    return NULL;
  const long* elt2verts = (const long*)PyArray_DATA(py_elt2verts);

  idx_t nn      = nbVerts;
  idx_t ne      = idx_t(PyArray_DIM(py_elt2verts,0));
  std::cout << "nn = " << nn << std::endl;
  std::cout << "ne = " << ne << std::endl;
  idx_t* eptr   = new idx_t[ne+1];
  eptr[0] = 0;
  for ( long i = 1; i <= ne; ++i ) eptr[i] = eptr[i-1]+3;
  idx_t* eind   = new idx_t[3*ne];
  std::copy(elt2verts, elt2verts + 3*ne, eind);
  idx_t ncommon = 2;
  idx_t nparts  = nbSubDoms;
  idx_t objval;
  idx_t* epart  = new idx_t[ne];
  idx_t* npart  = new idx_t[nn];
  int ok = METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL,
			      &ncommon, &nparts, NULL, NULL, &objval,
			      epart, npart);
  /* On compte le nombre d'éléments par domaine : */
  long* nbEltsPerDoms = new long[nbSubDoms];
  std::fill_n(nbEltsPerDoms, nbSubDoms, 0);
  for ( long i = 0; i < ne; ++i ) {
    nbEltsPerDoms[epart[i]] += 1;
  }
  PyObject* pyLst = PyList_New(nbSubDoms);
  PyArrayObject* indEltSubDoms;
  long** ptIndElts = new long*[nbSubDoms];
  for ( long i = 0; i < nbSubDoms; ++i ) {
    npy_intp nbElts = nbEltsPerDoms[i];
    indEltSubDoms = (PyArrayObject*)PyArray_SimpleNew(1,&nbElts,PyArray_LONG);
    PyList_SetItem(pyLst, i, (PyObject*)indEltSubDoms);
    ptIndElts[i] = (long*)PyArray_DATA(indEltSubDoms);
  }
  std::fill_n(nbEltsPerDoms, nbSubDoms, 0);
  for ( long i = 0; i < ne; ++i ) {
    ptIndElts[epart[i]][nbEltsPerDoms[epart[i]]] = i;
    nbEltsPerDoms[epart[i]] += 1;
  }
  
  delete [] nbEltsPerDoms;
  delete [] npart; delete [] epart;
  delete [] eind ; delete [] eptr ;
  return Py_BuildValue("N", pyLst);
}
// ========================================================================
static PyMethodDef Py_Methods[] =
  {
    {"splitNodeMesh", py_splitNodeMesh, METH_VARARGS, doc_splitNodeMesh},
    {"splitEltMesh" , py_splitEltMesh , METH_VARARGS, doc_splitEltMesh },
    {NULL, NULL} /* Guards */
  };
// ========================================================================
#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef splittermoduledef = {
    PyModuleDef_HEAD_INIT,
    "splitter",                      // Nom du module
    "Partitionner of mesh using metis",   // Documentation du module,
    -1                      ,   // Taille du module ?
    Py_Methods              ,   // Fonctions du module
    NULL,                       // m_reload ?
    NULL,                       // m_traverse ?
    NULL,                       // m_clear
    NULL,                       // m_free
  };
PyMODINIT_FUNC PyInit_splitter(void)
{
  PyObject* m; 
  m = PyModule_Create(&splittermoduledef);
  if (m==NULL) return NULL;
  /*  important : initialize numpy to use in subroutines !!!! */
  import_array();
  return m;
}
#else 
static char splitter_doc[] = "Partitionner of mesh using metis";
PyMODINIT_FUNC initsplitter()
{
  PyObject* m; 
  m = Py_InitModule4("splitter",Py_Methods, splitter_doc,
		     (PyObject*)NULL,PYTHON_API_VERSION);
  if (m==NULL) return;
  /*  important : initialize numpy to use in subroutines !!!! */
  import_array();
}
#endif
