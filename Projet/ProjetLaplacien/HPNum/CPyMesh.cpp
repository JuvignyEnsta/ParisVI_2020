# include <Python.h>
# include <math.h>
# include <algorithm>
# include <iostream>
# include <fstream>
# include <cstring>
# include <cassert>
# include <list>
# include <numpy/arrayobject.h>

static char compvert2elts_doc[] = 
  "Usage : comvert2elts( elt2verts )"
  "Calcul la connectivité du maillage des sommets vers les éléments."
  "elt2verts : Connectivité éléments vers sommets du maillage.\n";
static PyObject* py_compvert2elts( PyObject* self, PyObject* args )
{
  PyArrayObject *py_elt2verts, *py_begVert2elts, *py_vert2elts;
  PyArrayObject *py_glob2loc = NULL;
  const long* glob2loc = NULL;
  long nbVerts;
  if (!PyArg_ParseTuple( args, "O!|O!", &PyArray_Type, &py_elt2verts,
			 &PyArray_Type, &py_glob2loc))
    return NULL;
  if ( py_glob2loc != NULL ) glob2loc = (const long*)PyArray_DATA(py_glob2loc);
  long nbElts =  long(PyArray_DIM(py_elt2verts,0))-1;
  // Calcul du nombre de noeuds pris en compte par elt2verts :
  std::list<long> indVerts;
  for ( long i = 0; i < nbElts; ++i ) {
    const long* elt2verts = (const long*)PyArray_GETPTR2(py_elt2verts,i,0);
    long v1 = elt2verts[0];
    long v2 = elt2verts[1];
    long v3 = elt2verts[2];
    //std::cerr << "v1,v2, v3 : " << v1 << ", " << v2 << ", " << v3 << std::endl;
    assert(v1>=0); assert(v2 >= 0); assert(v3 >= 0);
    if ( glob2loc != NULL ) {
      if (glob2loc[v1] >= 0) indVerts.push_back(v1);
      if (glob2loc[v2] >= 0) indVerts.push_back(v2);
      if (glob2loc[v3] >= 0) indVerts.push_back(v3);
    } else {
      indVerts.push_back(elt2verts[0]);
      indVerts.push_back(elt2verts[1]);
      indVerts.push_back(elt2verts[2]);
    }
  }
  indVerts.sort();
  indVerts.unique();
  nbVerts = indVerts.size();
  
  npy_intp nd[2];
  nd[0] = npy_intp(nbVerts+1);
  py_begVert2elts = (PyArrayObject*)PyArray_SimpleNew(1,nd,PyArray_LONG);
  long* begVert2elts = (long*)PyArray_DATA(py_begVert2elts);
  long* nbEltsPerVert = new long[nbVerts];
  std::fill_n(nbEltsPerVert, nbVerts, 0);
  for ( long i = 0; i < nbElts; ++i ) {
    const long* elt2verts = (const long*)PyArray_GETPTR2(py_elt2verts,i,0);
    long v1 = elt2verts[0];
    long v2 = elt2verts[1];
    long v3 = elt2verts[2];
    assert(v1>=0); assert(v2 >= 0); assert(v3 >= 0);
    if ( glob2loc != NULL ) {
      if (glob2loc[v1] >= 0) nbEltsPerVert[glob2loc[v1]] += 1;
      if (glob2loc[v2] >= 0) nbEltsPerVert[glob2loc[v2]] += 1;
      if (glob2loc[v3] >= 0) nbEltsPerVert[glob2loc[v3]] += 1;
    } else {
      nbEltsPerVert[v1] += 1;
      nbEltsPerVert[v2] += 1;
      nbEltsPerVert[v3] += 1;      
    }
  }
  begVert2elts[0] = 0;
  for ( long iV = 0; iV < nbVerts; ++iV ) {
    begVert2elts[iV+1] = begVert2elts[iV] + nbEltsPerVert[iV];
  }
  nd[0] = begVert2elts[nbVerts];
  py_vert2elts = (PyArrayObject*)PyArray_SimpleNew(1,nd,PyArray_LONG);
  long* vert2elts = (long*)PyArray_DATA(py_vert2elts);
  
  std::fill_n(nbEltsPerVert, nbVerts, 0);
  for ( long i = 0; i < nbElts; ++i ) {
    const long* elt2verts = (const long*)PyArray_GETPTR2(py_elt2verts,i,0);
    long iV = (glob2loc == NULL ? elt2verts[0] : glob2loc[elt2verts[0]]);
    if (iV >= 0) {
      vert2elts[begVert2elts[iV] + nbEltsPerVert[iV]] = i;
      nbEltsPerVert[iV] += 1;
    }
    iV = (glob2loc == NULL ? elt2verts[1] : glob2loc[elt2verts[1]]);
    if (iV >= 0) {
      vert2elts[begVert2elts[iV] + nbEltsPerVert[iV]] = i;
      nbEltsPerVert[iV] += 1;
    }
    iV = (glob2loc == NULL ? elt2verts[2] : glob2loc[elt2verts[2]]);
    if (iV >= 0) {
      vert2elts[begVert2elts[iV] + nbEltsPerVert[iV]] = i;
      nbEltsPerVert[iV] += 1;
    }
  }
  delete [] nbEltsPerVert;
  return Py_BuildValue("NN", py_begVert2elts, py_vert2elts);
}
// ------------------------------------------------------------------------
static char readMesh[] =
  "Usage : read_fmtmesh(fileName)"
  "Lit un maillage à partir du fichier formaté donné par fileName"
  "Retourne un tuple (vertices, elt2verts) où :"
  " - vertices contient les coordonnées du sommet du maillage plus une coloration ( 0 par défaut, non nul pour Dirichlet )"
  " - elt2verts contient les indices des sommets définissant chaque élément du maillage\n";
static PyObject* py_read_fmtmesh( PyObject* self, PyObject* args )
{
  char buffer[4096];
  char* filename;
  double version;
  int file_type,  data_size;
  PyArrayObject *pyVertices = NULL, *pyElt2Verts = NULL;
  double* vertices = NULL;
  long* elt2verts = NULL;
  int nbBC = 0;
  
  if ( !PyArg_ParseTuple(args, "s", &filename) ) return NULL;

  std::ifstream input(filename);
  if (!input) {
    PyErr_SetString(PyExc_IOError, "Failed to open the file !");
    return NULL;
  }
  input.getline(buffer, 1024, '\n');
  input >> version >> file_type >> data_size;
  if ( (file_type != 0) || (data_size != 8) ) {
    PyErr_SetString(PyExc_IOError, "Wrong formated file format !");
    return NULL;
  }
  input.getline(buffer, 1024, '\n');
  input.getline(buffer, 1024, '\n');
  while (input) {
    input.getline(buffer, 1024, '\n');
    if (strncmp(buffer, "$Nodes",6) == 0) {
      double* vertices;
      long nbNodes, iNode;
      double x,y,z;
      input >> nbNodes;
      npy_intp dim[2];
      dim[0] = nbNodes;
      dim[1] = 4;
      pyVertices = (PyArrayObject*) PyArray_SimpleNew(2,dim,NPY_DOUBLE);
      vertices = (double*)PyArray_DATA(pyVertices);
    
      for ( long i = 0; i < nbNodes; ++i ) {
	input >> iNode >> x >> y >> z;
	iNode -= 1;
	assert(iNode < nbNodes);
	vertices[4*iNode  ] = x;
	vertices[4*iNode+1] = y;
	vertices[4*iNode+2] = z;
	vertices[4*iNode+3] = 0;
      }
      input.getline(buffer, 1024, '\n');
      input.getline(buffer, 1024, '\n');
      assert(strcmp(buffer, "$EndNodes") == 0);
    }// End if (Nodes
    else if (strncmp(buffer, "$Elements",9) == 0) {
      long nbElts, iElt, elt_type, nbTags, tag;
      input >> nbElts;
      pyElt2Verts = NULL;
      assert(pyVertices != NULL);
      vertices = (double*)PyArray_DATA(pyVertices);
      for ( long i = 0; i < nbElts; ++i ) {
	input >> iElt >> elt_type >> nbTags;
	iElt -= 1;
	for ( long t = 0; t < nbTags; ++t ) {
	  input >> tag;
	}
	if ( elt_type == 1 ) // C'est une ligne "physique", donc une condition limite :
	  {
	    long nd1, nd2;
	    input >> nd1 >> nd2;
	    nd1 -= 1; nd2 -= 1;
	    vertices[4*nd1+3] = tag;
	    vertices[4*nd2+3] = tag;
	    nbBC += 1;
	  }
	else if ( elt_type == 2 ) // C'est un triangle :
	  {
	    iElt -= nbBC;
	    if ( pyElt2Verts == NULL ) {
	      npy_intp dim[2];
	      dim[0] = nbElts-nbBC;
	      dim[1] = 3;
	      pyElt2Verts = (PyArrayObject*) PyArray_SimpleNew(2,dim,NPY_LONG);
	      elt2verts = (long*)PyArray_DATA(pyElt2Verts);
	    }
	    input >> elt2verts[3*iElt+0] >> elt2verts[3*iElt+1] >> elt2verts[3*iElt+2];
	    elt2verts[3*iElt+0] -= 1; elt2verts[3*iElt+1] -= 1; elt2verts[3*iElt+2] -= 1;
	  }
      }
      input.getline(buffer, 1024, '\n');
      input.getline(buffer, 1024, '\n');
    }
    else { // Autre, on va "skipper"
      std::string flag(buffer);
      if ( flag.size()  > 0 ) {
      std::string endBuffer = "$End" + flag.substr(1);
      do {
	input.getline(buffer, 1024, '\n');
      } while ( ( std::string(buffer).substr(0,4) != endBuffer.substr(0,4) ) && (input) );
      }
    }
  }
  return Py_BuildValue("NN", pyVertices, pyElt2Verts);
}
// ========================================================================
static PyMethodDef Py_Methods[] =
  {
    {"compvert2elts", py_compvert2elts, METH_VARARGS, compvert2elts_doc},
    {"read", py_read_fmtmesh, METH_VARARGS, readMesh},
    {NULL, NULL} /* Guards */
  };
// ========================================================================
#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef meshmoduledef = {
    PyModuleDef_HEAD_INIT,
    "mesh",                      // Nom du module
    "Mesh utilities",   // Documentation du module,
    -1                      ,   // Taille du module ?
    Py_Methods              ,   // Fonctions du module
    NULL,                       // m_reload ?
    NULL,                       // m_traverse ?
    NULL,                       // m_clear
    NULL,                       // m_free
  };
PyMODINIT_FUNC PyInit_mesh(void)
{
  PyObject* m; 
  m = PyModule_Create(&meshmoduledef);
  if (m==NULL) return NULL;
  /*  important : initialize numpy to use in subroutines !!!! */
  import_array();
  return m;
}
#else  
static char mesh_doc[] = "Mesh utilities";
PyMODINIT_FUNC initmesh()
{
  PyObject* m; 
  m = Py_InitModule4("mesh",Py_Methods, mesh_doc,
		     (PyObject*)NULL,PYTHON_API_VERSION);
  if (m==NULL) return;
  /*  important : initialize numpy to use in subroutines !!!! */
  import_array();
}
#endif
