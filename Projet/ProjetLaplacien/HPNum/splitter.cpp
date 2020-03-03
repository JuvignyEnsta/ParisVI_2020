# include <Python.h>
# include <cmath>
# include <vector>
# include <iostream>
# include <algorithm>
# include <numpy/arrayobject.h>
#include "hilbert_curve.hpp"

struct point2D
{
    double x,y;
    point2D(double xc, double yc) : x(xc), y(yc)
    {}

};

inline point2D mid_point( const point2D& p1, const point2D& p2 )
{
    return {0.5*(p1.x+p2.x), 0.5*(p1.y+p2.y)};
}

static std::pair<point2D,point2D> compute_bounding_box( const unsigned long nb_verts, const double* coords, int stride=4)
{
    point2D cmin(coords[0], coords[1]);
    point2D cmax(coords[0], coords[1]);
    for ( unsigned long i = 1; i < nb_verts; ++i )
    {
        cmin.x = std::min(cmin.x, coords[stride*i+0]);
        cmin.y = std::min(cmin.y, coords[stride*i+1]);
        cmax.x = std::max(cmax.x, coords[stride*i+0]);
        cmax.y = std::max(cmax.y, coords[stride*i+1]);
    }
    return {cmin, cmax};
}

using morton_number = std::pair<unsigned long long, unsigned long long>;


static std::vector<morton_number>
compute_morton_ordering( const unsigned long nb_verts, const double* coords,
                         const std::pair<point2D,point2D>& bbox, int N, int stride=4 )
{
    std::vector<morton_number> numbering;
    numbering.reserve(nb_verts);
    auto hilbert_curve = HilbertCurve(N, 2);
    unsigned long long pN = (1<<N)-1;
    point2D lgth{ bbox.second.x - bbox.first.x, bbox.second.y - bbox.first.y };
    for ( unsigned long iVert = 0; iVert < nb_verts; ++iVert )
    {
        std::pair<point2D,point2D> cur_bbox = bbox;
        point2D center_box = mid_point(cur_bbox.first, cur_bbox.second);
        morton_number index;
        index.first = iVert;
        double xc = coords[stride*iVert+0];
        double yc = coords[stride*iVert+1];
        double xn = (xc-bbox.first.x)/lgth.x;
        double yn = (yc-bbox.first.y)/lgth.y;
        unsigned long long ix = (unsigned long long)(pN*xn);
        unsigned long long iy = (unsigned long long)(pN*yn);
        assert(ix>=0);
        assert(iy>=0);
        assert(ix<=pN);
        assert(iy<=pN);

        index.second= hilbert_curve.distance_from_coordinates({ix,iy});
/*
        for ( int i = 0; i < N; ++i )
        {
            int msk = index.second&3;
            int ry = 0, rx = 0;
            if ((msk == 3)  or (msk == 2)) rx = 1;
            if ((msk == 1) or (msk == 2)) ry = 1;
            index.second <<= 2;
            bool is_right = (xc > center_box.x);
            bool is_up    = (yc > center_box.y);
            int inc_x = 0, inc_y = 0;
            // Enregistre le quadrant suivant l'ordonnancement de Morton
            if (is_right) 
            {
                cur_bbox.first.x = center_box.x;
                inc_x = 1;
            }
            else
            {
                cur_bbox.second.x= center_box.x;
            }

            if (is_up) 
            {
                cur_bbox.first.y = center_box.y;
                inc_y = 1;
            }
            else
                cur_bbox.second.y= center_box.y;
            center_box = mid_point(cur_bbox.first, cur_bbox.second);
            rot(2, inc_x, inc_y, rx, ry);
            if (inc_x==1)
                if (inc_y==0)
                    index.second += 3;
                else
                    index.second += 2;
            else if ( inc_y == 1 ) index.second += 1;
        }*/
        numbering.emplace_back(index);
    }
    return numbering;
}

static char doc_splitNodeMesh[] =
  "Usage : splitter.NodeMesh( nbDoms, vertices_coords )"
  "where nbDoms is the number of subdomains"
  "and vertices_coords are the coordinates of the vertice of the mesh."
  "return set of nodes per domains\n";
static PyObject* py_split_nodeMesh( PyObject* self, PyObject* args )
{
    PyArrayObject *py_coords;
    PyObject *splitted_vertices;
    int nb_subdomains;
    if ( !PyArg_ParseTuple( args, "iO!", &nb_subdomains, &PyArray_Type, &py_coords) ) return NULL;
    const double* coords = (const double*)PyArray_DATA(py_coords);

    npy_intp nbVerts = PyArray_DIM(py_coords,0);
    auto bbox = compute_bounding_box(nbVerts, coords);
    // Calcul de l'index de morton pour chaque sommet :
    auto morton_numbering = compute_morton_ordering(nbVerts, coords, bbox, 20);
    std::sort(morton_numbering.begin(), morton_numbering.end(), [](const morton_number& n1,
                                                                   const morton_number& n2){
        return n1.second < n2.second;
    });
    unsigned long nb_loc_verts = nbVerts/unsigned(nb_subdomains);
    unsigned long suppl_verts  = nbVerts%unsigned(nb_subdomains);
    splitted_vertices = PyList_New(nb_subdomains);
    unsigned long long start_indices = 0;
    for ( int iDom = 0; iDom < nb_subdomains; ++iDom )
    {
        npy_intp nbLocVerts = nb_loc_verts;
        if (iDom < suppl_verts) nbLocVerts += 1;
        PyArrayObject* loc_indices = (PyArrayObject*)PyArray_EMPTY(1, &nbLocVerts, NPY_UINT64, 0);
        unsigned long long *data_indices = (unsigned long long*)PyArray_DATA(loc_indices);
        for ( unsigned long long i_ind = 0; i_ind < nbLocVerts; ++i_ind )
        {
            data_indices[i_ind] = morton_numbering[start_indices + i_ind].first;
        }
        start_indices += nbLocVerts;
        PyList_SetItem(splitted_vertices, iDom, (PyObject*)loc_indices);
    }
    return splitted_vertices;
}
// 
static char doc_splitEltMesh[] =
  "Usage : splitter.EltMesh( nbDoms, (elt2vert,vertices_coords) )"
  "where nbDoms is the number of subdomains"
  "and vertices_coords are the coordinates of the vertice of the mesh."
  "return set of nodes per domains\n";
static PyObject* py_split_eltMesh( PyObject* self, PyObject* args )
{
    PyArrayObject *py_coords, *py_elt2vert;
    PyObject *splitted_elements;
    int nb_subdomains;
    if ( !PyArg_ParseTuple( args, "i(O!O!)", &nb_subdomains, &PyArray_Type, &py_elt2vert, &PyArray_Type, &py_coords) ) 
        return NULL;
    const long* elt2vert = (const long*)PyArray_DATA(py_elt2vert);
    const double* coords = (const double*)PyArray_DATA(py_coords);

    npy_intp nbVerts = PyArray_DIM(py_coords,0);
    npy_intp nbElts  = PyArray_DIM(py_elt2vert,0);
    auto bbox = compute_bounding_box(nbVerts, coords);
    // Calcul du barycentre de chaque élément du maillage :
    std::vector<double> bary_coords(2*nbElts);
    for ( unsigned long long iElt = 0; iElt < nbElts; ++iElt )
    {
        long nd_1 = elt2vert[3*iElt+0], nd_2 = elt2vert[3*iElt+1] , nd_3 = elt2vert[3*iElt+2];
        bary_coords[2*iElt+0] = (coords[4*nd_1+0]+coords[4*nd_2+0]+coords[4*nd_3+0])/3.;
        bary_coords[2*iElt+1] = (coords[4*nd_1+1]+coords[4*nd_2+1]+coords[4*nd_3+1])/3.;
    }
    // Calcul de l'index de morton pour chaque sommet :
    auto morton_numbering = compute_morton_ordering(nbElts, bary_coords.data(), bbox, 20, 2);
    // A partir de la renumerotation des noeuds, on reordonne les éléments :
    std::sort(morton_numbering.begin(), morton_numbering.end(), [](const morton_number& n1,
                                                                   const morton_number& n2){
        return n1.second < n2.second;
    });
    unsigned long nb_loc_elts = nbElts/unsigned(nb_subdomains);
    unsigned long suppl_elts  = nbElts%unsigned(nb_subdomains);
    splitted_elements = PyList_New(nb_subdomains);
    unsigned long long start_indices = 0;
    for ( int iDom = 0; iDom < nb_subdomains; ++iDom )
    {
        npy_intp nbLocElts = nb_loc_elts;
        if (iDom < suppl_elts) nbLocElts += 1;
        PyArrayObject* loc_indices = (PyArrayObject*)PyArray_EMPTY(1, &nbLocElts, NPY_UINT64, 0);
        unsigned long long *data_indices = (unsigned long long*)PyArray_DATA(loc_indices);
        for ( unsigned long long i_ind = 0; i_ind < nbLocElts; ++i_ind )
        {
            data_indices[i_ind] = morton_numbering[start_indices + i_ind].first;
        }
        start_indices += nbLocElts;
        PyList_SetItem(splitted_elements, iDom, (PyObject*)loc_indices);
    }
    return splitted_elements;
}
// ========================================================================
static PyMethodDef Py_Methods[] =
  {
    {"node"    , py_split_nodeMesh, METH_VARARGS, doc_splitNodeMesh},
    {"element" , py_split_eltMesh , METH_VARARGS, doc_splitEltMesh },
    {NULL, NULL} /* Guards */
  };
// ========================================================================
  static struct PyModuleDef splittermoduledef = {
    PyModuleDef_HEAD_INIT,
    "splitter",                      // Nom du module
    "Partitionner of mesh",   // Documentation du module,
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
