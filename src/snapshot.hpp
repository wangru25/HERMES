#include "alpha_shape.hpp"

class EdgeRadiusCompare{
public:
    AlphaShape* m_alpha;

    EdgeRadiusCompare(AlphaShape* alpha){
        m_alpha = alpha;
    }

    bool operator()(const EdgeIterator& ei1, const EdgeIterator& ei2) const{
        return m_alpha->edgeAlpha(ei1) < m_alpha->edgeAlpha(ei2);
    }
};

class FacetRadiusCompare{
public:
    AlphaShape* m_alpha;

    FacetRadiusCompare(AlphaShape* alpha){
        m_alpha = alpha;
    }

    bool operator()(const FacetIterator& fi1, const FacetIterator& fi2) const{
        return m_alpha->facetAlpha(fi1) < m_alpha->facetAlpha(fi2);
    }
};

class CellRadiusCompare{
public:
    AlphaShape* m_alpha;

    CellRadiusCompare(AlphaShape* alpha){
        m_alpha = alpha;
    }

    bool operator()(const CellIterator& ci1, const CellIterator& ci2) const{
        return m_alpha->cellAlpha(ci1) < m_alpha->cellAlpha(ci2);
    }
};

typedef std::priority_queue<EdgeIterator, std::vector<EdgeIterator>, EdgeRadiusCompare> EdgePriorityQueue;
typedef std::priority_queue<FacetIterator, std::vector<FacetIterator>, FacetRadiusCompare> FacetPriorityQueue;
typedef std::priority_queue<CellIterator, std::vector<CellIterator>, CellRadiusCompare> CellPriorityQueue;

class Snapshot{
protected:
    AlphaShape m_alpha;
    
    std::unordered_map<EdgeIterator, double, EdgeIteratorHash> m_edge_radius;
    std::unordered_map<FacetIterator, double, FacetIteratorHash> m_facet_radius;
    std::unordered_map<CellIterator, double, CellIteratorHash> m_cell_radius;

    double m_edge_length_min;
    double m_edge_length_max;

    std::vector<EdgeIterator> m_sorted_edges;
    std::vector<FacetIterator> m_sorted_facets;
    std::vector<CellIterator> m_sorted_cells;

    std::unordered_map<EdgeIterator, int, EdgeIteratorHash> m_sorted_edge_idx;
    std::unordered_map<FacetIterator, int, FacetIteratorHash> m_sorted_facet_idx;
    std::unordered_map<CellIterator, int, CellIteratorHash> m_sorted_cell_idx;

    // min heap
    EdgePriorityQueue* m_edge_heap;
    FacetPriorityQueue* m_facet_heap;
    CellPriorityQueue* m_cell_heap;

    std::vector<double> m_filtration;
    int m_num_eigenvalues;
    double m_p = 0.; // p-persistence
    char m_complex = 'a'; // 'a' for alpha complex and 'r' for rips complex
    // TODO: Czech complex

    std::unique_ptr<matlab::engine::MATLABEngine> m_matlab_engine;

    std::vector<std::vector<double>> m_vertex_snapshots;
    std::vector<std::vector<double>> m_edge_snapshots;
    std::vector<std::vector<double>> m_facet_snapshots;

public:
    // read a point cloud file to build the alpha shape
    void initializeParameters(int num_eigenvalues);
    void initializeParameters(int num_eigenvalues, double p);
    void initializeParameters(int num_eigenvalues, double p, char complex);
    void readFiltration(const std::string& filename);
    void buildAlphaShape(const std::string& filename);
    void preprocess();
    void takeSnapshots();
    void write();
    void debug();

protected:
    void computeMinMaxEdgeLength();
    void computeRadius();
    void buildPriorityQueue();
    void assembleSortedElements();
    void indexSortedElements();

protected:
    void computeEdgeRadius();
    void computeFacetRadius();
    void computeCellRadius();

    void buildEdgePriorityQueue();
    void buildFacetPriorityQueue();
    void buildCellPrioiryQueue();

    void assembleSortedEdges();
    void assembleSortedFacets();
    void assembleSortedCells();

    void indexSortedEdges();
    void indexSortedFacets();
    void indexSortedCells();

    void writeVertexSnapshots();
    void writeEdgeSnapshots();
    void writeFacetSnapshots();

protected:
    void determineNextEdgeSize(const int& current_size, int& next_size, double filtration);
    void determineNextFacetSize(const int& current_size, int& next_size, double filtration);
    void determineNextCellSize(const int& current_size, int& next_size, double filtration);

    // build exteior derivative operator (transpose of boundary operator)
    void fillED0T(SparseMatrix& ED0T, const int& current_size, const int& next_size); 
    void fillED1T(SparseMatrix& ED1T, const int& current_size, const int& next_size);
    void fillED2T(SparseMatrix& ED2T, const int& current_size, const int& next_size);

    void buildLaplacian0(SparseMatrix& L0, const SparseMatrix& ED0T, const int& edge_size);
    void buildLaplacian1(SparseMatrix& L1, const SparseMatrix& ED0T, const SparseMatrix& ED1T, const int& edge_size, const int& facet_size);
    void buildLaplacian2(SparseMatrix& L2, const SparseMatrix& ED1T, const SparseMatrix& ED2T, const int& edge_size, const int& facet_size, const int& cell_size);
    
    void buildPersistentLaplacian0(SparseMatrix& L0, const SparseMatrix& ED0T, const SparseMatrix& ED1T, const int& edge_size, const double filtration, double p);
	void buildPersistentLaplacian1(SparseMatrix& L1, const SparseMatrix& ED0T, const SparseMatrix& ED1T, const int& edge_size, const int& facet_size, const double filtration, double p);
    void buildPersistentLaplacian2(SparseMatrix& L2, const SparseMatrix& ED1T, const SparseMatrix& ED2T, const int& edge_size, const int& facet_size, const int& cell_size, const double filtration, double p);

    void takeVertexSnapshots(const SparseMatrix& L0, const int matrix_size); // beta_0
    void takeEdgeSnapshots(const SparseMatrix& L1, const int matrix_size); // beta_1
    void takeFacetSnapshots(const SparseMatrix& L2, const int matrix_size); // beta_2

protected:
    void startMatlab();
    void matlabEIGS(std::vector<double>& eval, std::vector<ColumnVector>& evec, int matrix_size,
        std::vector<double>& row, std::vector<double>& col, std::vector<double>& val, int m);

private:
    //debug
    void writeVertices();

    void testEdgeFiltration();
    void testFacetFiltration();
    void testCellFiltration();
};
