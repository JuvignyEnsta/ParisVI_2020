# include <Python.h>
# include <math.h>
# include <algorithm>
# include <iostream>
# include <numpy/arrayobject.h>

static char compElemMatrix_doc[] =
  "Compute the elementary matrix of a laplacian in 2D.";
static PyObject* py_compElemMatrix( PyObject* self, PyObject* args )
{
  double v1x, v1y, v2x, v2y, v3x, v3y;
  npy_intp nd[2];

  if ( !PyArg_ParseTuple(args, "(dd)(dd)(dd)", &v1x, &v1y, &v2x, &v2y,
			 &v3x, &v3y ) ) return NULL;
  nd[0] = 3; nd[1] = 3;
  PyArrayObject* pyAElt = (PyArrayObject*)PyArray_SimpleNew(2,nd,PyArray_DOUBLE);
  if (pyAElt == NULL ) {
    return NULL;
  }
  double A1A3[2], A1A2[2], A2A3[2];
  A1A3[0] = v3x - v1x; A1A3[1] = v3y - v1y;
  A1A2[0] = v2x - v1x; A1A2[1] = v2y - v1y;
  A2A3[0] = v3x - v2x; A2A3[1] = v3y - v2y;

  double SqA1A3 = A1A3[0]*A1A3[0] + A1A3[1]*A1A3[1];
  double A1A2dotA1A3 = A1A2[0]*A1A3[0] + A1A2[1]*A1A3[1];
  double SqA2A3 = A2A3[0]*A2A3[0] + A2A3[1]*A2A3[1];
  double* AElt = (double*)PyArray_DATA(pyAElt);
  std::fill_n( AElt, 9, 0.);
  // Premiere ligne :
  AElt[0] =  0.5*(SqA1A3 -2.*A1A2dotA1A3 + SqA2A3 );
  AElt[1] = -0.5*( SqA1A3 - A1A2dotA1A3 ); 
  AElt[2] = -0.5*( SqA2A3 - A1A2dotA1A3 );
  // Deuxieme ligne :
  AElt[3] = AElt[1];
  AElt[4] =  0.5*SqA1A3; 
  AElt[5] = -0.5*A1A2dotA1A3;
  // Troisième ligne :
  AElt[6] = AElt[2];
  AElt[7] = AElt[5];
  AElt[8] =  0.5*SqA2A3;
  PyObject* ret =  Py_BuildValue("O",(PyObject*)pyAElt);
  Py_DECREF(pyAElt);
  return ret;
}
// ========================================================================
static PyMethodDef Py_Methods[] =
  {
    {"comp_eltmat", py_compElemMatrix, METH_VARARGS, compElemMatrix_doc},
    {NULL, NULL} /* Guards */
  };
// ========================================================================
#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef laplacianmoduledef = {
    PyModuleDef_HEAD_INIT,
    "laplacian",                      // Nom du module
    "Laplacian equation discretization",   // Documentation du module,
    -1                      ,   // Taille du module ?
    Py_Methods              ,   // Fonctions du module
    NULL,                       // m_reload ?
    NULL,                       // m_traverse ?
    NULL,                       // m_clear
    NULL,                       // m_free
  };
PyMODINIT_FUNC PyInit_laplacian(void)
{
  PyObject* m; 
  m = PyModule_Create(&laplacianmoduledef);
  if (m==NULL) return NULL;
  /*  important : initialize numpy to use in subroutines !!!! */
  import_array();
  return m;
}
#else  
static char laplacian_doc[] = "Finite Element Method";
PyMODINIT_FUNC initlaplacian()
{
  PyObject* m; 
  m = Py_InitModule4("laplacian",Py_Methods, laplacian_doc,
		     (PyObject*)NULL,PYTHON_API_VERSION);
  if (m==NULL) return;
  /*  important : initialize numpy to use in subroutines !!!! */
  import_array();
}
#endif
