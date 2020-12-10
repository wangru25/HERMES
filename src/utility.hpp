#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <limits>
#include <random>
#include <chrono>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point3;
typedef Kernel::Vector_3 Vector3;

typedef CGAL::Triangulation_vertex_base_3<Kernel> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<Kernel> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb> Tds;

typedef CGAL::Delaunay_triangulation_3<Kernel,Tds,CGAL::Fast_location> Delaunay;
typedef Delaunay::Vertex_iterator VertexIterator;
typedef Delaunay::Edge_iterator EdgeIterator;
typedef Delaunay::Facet_iterator FacetIterator;
typedef Delaunay::Cell_iterator CellIterator;
typedef Delaunay::Facet_circulator FacetCirculator;
typedef Delaunay::Cell_circulator CellCirculator;
typedef Delaunay::Facet Facet;
typedef Delaunay::Edge Edge;
typedef CGAL::Handle_hash_function CGALHash;

typedef CGALHash VertexIteratorHash;
class EdgeIteratorHash{
public:
    size_t operator()(const EdgeIterator& ei) const{
        return CGALHash()(ei->first) ^ std::hash<int>()(ei->second) ^ std::hash<int>()(ei->third);
    }
};
class FacetIteratorHash{
public:
    size_t operator()(const FacetIterator& fi) const{
        return CGALHash()(fi->first) ^ std::hash<int>()(fi->second);
    }
};
typedef CGALHash CellIteratorHash;

class EdgeHash{
public:
    size_t operator()(const Edge& e) const{
        return CGALHash()(e.first) ^ std::hash<int>()(e.second) ^ std::hash<int>()(e.third);
    }
};
class FacetHash{
public:
    size_t operator()(const Facet& f) const{
        return CGALHash()(f.first) ^ std::hash<int>()(f.second);
    }
};


typedef std::pair<VertexIterator, int> VertexAsBoundary;
typedef std::pair<EdgeIterator, int> EdgeAsBoundary;
typedef std::pair<FacetIterator, int> FacetAsBoundary;

class VertexAsBoundaryHash{
public:
    size_t operator()(const VertexAsBoundary& v) const{
        return CGALHash()(v.first) ^ std::hash<int>()(v.second);
    }
};
class EdgeAsBoundaryHash{
public:
    size_t operator()(const EdgeAsBoundary& e) const{
        return CGALHash()(e.first) ^ std::hash<int>()(e.second);
    }
};
class FacetAsBoundaryHash{
public:
    size_t operator()(const FacetAsBoundary& f) const{
        return CGALHash()(f.first) ^ std::hash<int>()(f.second);
    }
};

typedef std::unordered_set<VertexAsBoundary, VertexAsBoundaryHash> BoundaryVertexSet;
typedef std::unordered_set<EdgeAsBoundary, EdgeAsBoundaryHash> BoundaryEdgeSet;
typedef std::unordered_set<FacetAsBoundary, FacetAsBoundaryHash> BoundaryFacetSet;

typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrix;
typedef Eigen::MatrixXd DenseMatrix;
typedef Eigen::VectorXd ColumnVector;

const std::string CONSOLE_RED("\033[0;31m");
const std::string CONSOLE_GREEN("\033[1;32m");
const std::string CONSOLE_YELLOW("\033[1;33m");
const std::string CONSOLE_CYAN("\033[0;36m");
const std::string CONSOLE_MAGENTA("\033[0;35m");
const std::string CONSOLE_WHITE("\033[0m");

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::duration<float> TimeSecond;

#endif