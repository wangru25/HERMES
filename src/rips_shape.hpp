#ifndef _RIPS_H_
#define _RIPS_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <numeric>

class SimpleVertex{
public:
    int index;
    int *adj_vertices; //adjacency list in increasing distance
    
    SimpleVertex(int ind) : index{ind} {};
    SimpleVertex() = default;
};
class VertexCompare{
public:
    bool operator()(const SimpleVertex& v1, const SimpleVertex& v2) const{
        return v1.index < v2.index;
    }
};
class SimpleEdge{
public:
    int index;
    int v1;
    int v2;
    double length;
    
    SimpleEdge(int ind, int i, int j, double l)
    : index{ind}, v1{i}, v2{j}, length{l} {};
    SimpleEdge() = default;
    
    bool operator==(SimpleEdge other_edge) const { return this->index == other_edge.index; }
};
class EdgeCompare{
public:
    bool operator()(const SimpleEdge& e1, const SimpleEdge& e2) const{
        return e1.index < e2.index;
    }
};
class SimpleFacet{
public:
    int index;
    int v1;
    int v2;
    int v3;
    double maxdistance;
    
    SimpleFacet(int ind, int x, int y, int z, double max)
    : index{ind}, v1{x}, v2{y}, v3{z}, maxdistance{max} {};
    SimpleFacet() = default;
    
    bool operator==(SimpleFacet other_facet) const { return this->index == other_facet.index; }
};
class FacetCompare{
public:
    bool operator()(const SimpleFacet& f1, const SimpleFacet& f2) const{
        return f1.index < f2.index;
    }
};
class SimpleCell{
public:
    int index;
    int v1;
    int v2;
    int v3;
    int v4;
    double maxdistance;
    
    SimpleCell(int ind, int x, int y, int z, int w, double max)
    : index{ind}, v1{x}, v2{y}, v3{z}, v4{w}, maxdistance{max} {};
    SimpleCell() = default;
};
class CellCompare{
public:
    bool operator()(const SimpleCell& c1, const SimpleCell& c2) const{
        return c1.index < c2.index;
    }
};


class RipsShape{
private:
    std::vector<std::vector<double>> m_point_cloud;
    std::vector<std::vector<double>> distance;
    const double distance_threshold = 11.; //TODO: Set by user input

    // basic
    std::map<SimpleVertex, int, VertexCompare> m_vertex_idx;
    std::map<SimpleEdge, int, EdgeCompare> m_edge_idx;
    std::map<SimpleFacet, int, FacetCompare> m_facet_idx;
    std::map<SimpleCell, int, CellCompare> m_cell_idx;

    std::map<int, SimpleVertex> m_idx_to_vertex;
    std::map<int, SimpleEdge> m_idx_to_edge;
    std::map<int, SimpleFacet> m_idx_to_facet;
    std::map<int, SimpleCell> m_idx_to_cell;

    int m_num_vertices;
    int m_num_edges;
    int m_num_facets;
    int m_num_cells;
 
    std::map<std::pair<int, int>, SimpleEdge> m_vertices_to_edge;
    std::map<std::pair<int, std::pair<int, int>>, SimpleFacet> m_vertices_to_facet;

public:
    RipsShape();
    ~RipsShape();

    // accessor
    int vertexIdx(const SimpleVertex& v) const;
    int edgeIdx(const SimpleEdge& e) const;
    int facetIdx(const SimpleFacet& f) const;
    int cellIdx(const SimpleCell& c) const;

    SimpleVertex idxToVertex(const int& idx) const;
    SimpleEdge idxToEdge(const int& idx) const;
    SimpleFacet idxToFacet(const int& idx) const;
    SimpleCell idxToCell(const int& idx) const;

    int numVertices() const;
    int numEdges() const;
    int numFacets() const;
    int numCells() const;

    // main functions
    void readPointCloud(std::string filename);
    void readDistanceMatrix(std::string filename);
    void convertToDistanceMatrix();
    void preprocess();
    void debug();
    
    std::vector<SimpleEdge> m_simple_edges;
    std::vector<SimpleFacet> m_simple_facets;
    std::vector<SimpleCell> m_simple_cells;
    
    std::map<int, int> positiveEdgeMap;
    std::map<int, int> negativeEdgeMap;
    std::map<int, std::pair<int, int>> positiveFacetMap;
    std::map<int, int> negativeFacetMap;
    std::map<int, std::pair<int, int>> positiveCellMap;
    std::map<int, std::pair<int, int>> negativeCellMap;

};

#endif
