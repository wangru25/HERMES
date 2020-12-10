#include "utility.hpp"

class AlphaShape{
private:
    std::vector<Point3> m_point_cloud;
    Delaunay m_triangulation;

    // basic
    std::unordered_map<VertexIterator, bool, VertexIteratorHash> m_valid_vertex;
    std::unordered_map<EdgeIterator, bool, EdgeIteratorHash> m_valid_edge;
    std::unordered_map<FacetIterator, bool, FacetIteratorHash> m_valid_facet;
    std::unordered_map<CellIterator, bool, CellIteratorHash> m_valid_cell;

    std::unordered_map<VertexIterator, int, VertexIteratorHash> m_vertex_idx;
    std::unordered_map<EdgeIterator, int, EdgeIteratorHash> m_edge_idx;
    std::unordered_map<FacetIterator, int, FacetIteratorHash> m_facet_idx;
    std::unordered_map<CellIterator, int, CellIteratorHash> m_cell_idx;
    
    std::unordered_map<int, VertexIterator> m_idx_to_vertex;
    std::unordered_map<int, EdgeIterator> m_idx_to_edge;
    std::unordered_map<int, FacetIterator> m_idx_to_facet;
    std::unordered_map<int, CellIterator> m_idx_to_cell;

    int m_num_vertices;
    int m_num_edges;
    int m_num_facets;
    int m_num_cells;
 
    std::unordered_map<Edge, EdgeIterator, EdgeHash> m_cell_edge_map;
    std::unordered_map<Facet, FacetIterator, FacetHash> m_cell_facet_map;

    std::unordered_map<EdgeIterator, double, EdgeIteratorHash> m_edge_length;

    std::unordered_map<EdgeIterator, Point3, EdgeIteratorHash> m_edge_circumcenter;
    std::unordered_map<FacetIterator, Point3, FacetIteratorHash> m_facet_circumcenter;
    std::unordered_map<CellIterator, Point3, CellIteratorHash> m_cell_circumcenter;

    std::unordered_map<EdgeIterator, double, EdgeIteratorHash> m_edge_radius;
    std::unordered_map<FacetIterator, double, FacetIteratorHash> m_facet_radius;
    std::unordered_map<CellIterator, double, CellIteratorHash> m_cell_radius;

    std::unordered_map<EdgeIterator, double, EdgeIteratorHash> m_edge_alpha;
    std::unordered_map<FacetIterator, double, FacetIteratorHash> m_facet_alpha;
    std::unordered_map<CellIterator, double, CellIteratorHash> m_cell_alpha;    

    // boundary
    std::unordered_map<EdgeIterator, BoundaryVertexSet, EdgeIteratorHash> m_edge_boundary;
    std::unordered_map<FacetIterator, BoundaryEdgeSet, FacetIteratorHash> m_facet_boundary;
    std::unordered_map<CellIterator, BoundaryFacetSet, CellIteratorHash> m_cell_boundary;

public:
    AlphaShape();
    ~AlphaShape();

    // accessor
    bool validVertex(const VertexIterator& vi) const;
    bool validEdge(const EdgeIterator& ei) const;
    bool validFacet(const FacetIterator& fi) const;
    bool validCell(const CellIterator& ci) const;

    int vertexIdx(const VertexIterator& vi) const;
    int edgeIdx(const EdgeIterator& ei) const;
    int facetIdx(const FacetIterator& fi) const;
    int cellIdx(const CellIterator& ci) const;

    VertexIterator idxToVertex(const int& idx) const;
    EdgeIterator idxToEdge(const int& idx) const;
    FacetIterator idxToFacet(const int& idx) const;
    CellIterator idxToCell(const int& idx) const;

    int numVertices() const;
    int numEdges() const;
    int numFacets() const;
    int numCells() const;

    EdgeIterator cellEdgeMap(const Edge& e) const;
    FacetIterator cellFacetMap(const Facet& f) const;

    double edgeLength(const EdgeIterator& ei) const;

    double edgeAlpha(const EdgeIterator& ei) const;
    double facetAlpha(const FacetIterator& fi) const;
    double cellAlpha(const CellIterator& ci) const;

    BoundaryVertexSet edgeBoundary(const EdgeIterator& ei) const;
    BoundaryEdgeSet facetBoundary(const FacetIterator& fi) const;
    BoundaryFacetSet cellBoundary(const CellIterator& ci) const;

    // main functions
    void readPointCloud(std::string filename);
    void buildTriangulation();
    void preprocess();
    void debug();

protected:
    void classifyElements(); // valid or invalid (have infinite vertex or not)
    void indexElements();
    void buildCellToEdgeFacetMap();
    void computeEdgeLength();
    void computeCircum();
    void computeAlpha();
    void computeBoundary();

protected:
    void classifyVertices();
    void classifyEdges();
    void classifyFacets();
    void classifyCells();

    void indexVertices();
    void indexEdges();
    void indexFacets();
    void indexCells();

    void buildCellToEdgeMap();
    void buildCellToFacetMap();

    void computeEdgeCircum();
    void computeFacetCircum();
    void computeCellCircum();

    void initAlpha();
    void processCellAlpha();
    void processFacetAlpha();
    void processEdgeAlpha();

    void computeEdgeBoundary();
    void computeFacetBoundary();
    void computeCellBoundary();

private:
    // debug functions
    void writeVertices();
    void writeEdges();
    void writeFacets();
    void writeCells();

    void writeSortedAlpha();
};


