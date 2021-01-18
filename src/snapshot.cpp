#include "snapshot.hpp"

void Snapshot::initializeParameters(int num_eigenvalues){
    m_num_eigenvalues = num_eigenvalues;
}

void Snapshot::initializeParameters(int num_eigenvalues, double p){
    m_num_eigenvalues = num_eigenvalues;
    m_p = p;
}

void Snapshot::initializeParameters(int num_eigenvalues, double p, char c){
    m_num_eigenvalues = num_eigenvalues;
    m_p = p;
    m_complex = c;
}

void Snapshot::readFiltration(const std::string& filename){
    std::ifstream file(filename);
    double val = 0;
    while(file>>val){
        m_filtration.push_back(val);
        //std::cout<<val<<std::endl;
    }
    file.close();
}

void Snapshot::buildAlphaShape(const std::string& filename){
    m_alpha.readPointCloud(filename);
    m_alpha.buildTriangulation();
    m_alpha.preprocess();
    m_alpha.debug();
}

void Snapshot::preprocess(){
    computeMinMaxEdgeLength();
    computeRadius();
    assembleSortedElements();
    indexSortedElements();
}

void Snapshot::takeSnapshots(){
    startMatlab();

    // for each new snapshot filtration value
    // we will first determine the new elements to be added
    // represented by elements in [current_element_size, next_element_size)
    int current_edge_size = 0;
    int current_facet_size = 0;
    int current_cell_size = 0;
    int next_edge_size;
    int next_facet_size;
    int next_cell_size;

    // Take care!!
    // 1.
    // We are building the transpose of exterior derivative operator
    // 2.
    // We can reserve number of non-zero elements for each column
    // Thus we can avoid reallocation when dynamically fill the matrices
    // We reserve 2 slots for the column of ED0T because each edge has 2 boundary vertices
    // Similar for others
    SparseMatrix ED0T(m_alpha.numVertices(), m_alpha.numEdges());
    SparseMatrix ED1T(m_alpha.numEdges(), m_alpha.numFacets());
    SparseMatrix ED2T(m_alpha.numFacets(), m_alpha.numCells());
    ED0T.reserve(Eigen::VectorXi::Constant(m_alpha.numEdges(), 2));
    ED1T.reserve(Eigen::VectorXi::Constant(m_alpha.numFacets(), 3));
    ED2T.reserve(Eigen::VectorXi::Constant(m_alpha.numCells(), 4));
    
    fillED0T(ED0T, 0, m_alpha.numEdges());
    fillED1T(ED1T, 0, m_alpha.numFacets());
    fillED2T(ED2T, 0, m_alpha.numCells());

    for(int i=0; i<m_filtration.size(); ++i){
        std::cout<<CONSOLE_GREEN<<"Filtration: "<<i<<CONSOLE_WHITE<<std::endl;

        std::cout<<CONSOLE_GREEN<<"Determine next size ... "<<CONSOLE_WHITE<<std::endl;
        determineNextEdgeSize(current_edge_size, next_edge_size, m_filtration[i]);
        determineNextFacetSize(current_facet_size, next_facet_size, m_filtration[i]);
        determineNextCellSize(current_cell_size, next_cell_size, m_filtration[i]);

        std::cout<<CONSOLE_GREEN<<"Fill boundary operators ... "<<CONSOLE_WHITE<<std::endl;
//        fillED0T(ED0T, current_edge_size, next_edge_size);
//        fillED1T(ED1T, current_facet_size, next_facet_size);
//        fillED2T(ED2T, current_cell_size, next_cell_size);

        std::cout<<CONSOLE_GREEN<<"Build Laplacians ... "<<CONSOLE_WHITE<<std::endl;
        SparseMatrix L0, L1, L2;
        // buildLaplacian0(L0, ED0T, next_edge_size);
        // buildLaplacian1(L1, ED0T, ED1T, next_edge_size, next_facet_size);
        // buildLaplacian2(L2, ED1T, ED2T, next_edge_size, next_facet_size, next_cell_size);
        buildPersistentLaplacian0(L0, ED0T, ED1T, next_edge_size, m_filtration[i], m_p);
        buildPersistentLaplacian1(L1, ED0T, ED1T, next_edge_size, next_facet_size, m_filtration[i], m_p);
        buildPersistentLaplacian2(L2, ED1T, ED2T, next_edge_size, next_facet_size, next_cell_size, m_filtration[i], m_p);

        std::cout<<CONSOLE_GREEN<<"Take snapshots ... "<<CONSOLE_WHITE<<std::endl;
        takeVertexSnapshots(L0, m_alpha.numVertices());
        takeEdgeSnapshots(L1, next_edge_size);
        takeFacetSnapshots(L2, next_facet_size);

        current_edge_size = next_edge_size;
        current_facet_size = next_facet_size;
        current_cell_size = next_cell_size;
    }
    
    //

    /*
    std::vector<EdgeIterator> edges;
    std::vector<FacetIterator> facets;
    std::vector<CellIterator> cells;
    edges.reserve(m_alpha.numEdges());
    facets.reserve(m_alpha.numFacets());
    cells.reserve(m_alpha.numCells());

    // Take care!!
    // 1.
    // We are building the transpose of exterior derivative operator
    // 2.
    // We can reserve number of non-zero elements for each column
    // Thus we can avoid reallocation when dynamically fill the matrices
    // We reserve 2 slots for the column of ED0T because each edge has 2 boundary vertices
    // Similar for others
    SparseMatrix ED0T(m_alpha.numVertices(), m_alpha.numEdges());
    SparseMatrix ED1T(m_alpha.numEdges(), m_alpha.numFacets());
    SparseMatrix ED2T(m_alpha.numFacets(), m_alpha.numCells());
    ED0T.reserve(Eigen::VectorXi::Constant(m_alpha.numEdges(), 2));
    ED1T.reserve(Eigen::VectorXi::Constant(m_alpha.numFacets(), 3));
    ED2T.reserve(Eigen::VectorXi::Constant(m_alpha.numCells(), 4));

    std::cout<<CONSOLE_RED<<"Begin filtration ..."<<CONSOLE_WHITE<<std::endl;
    for(int i=0; i<m_filtration.size(); ++i){
        std::cout<<i<<std::endl;
        std::cout<<CONSOLE_GREEN<<"Get additional edges ..."<<CONSOLE_WHITE<<std::endl;
        edges.clear();
        facets.clear();
        cells.clear();
        getAdditionalEdgeAtFiltration(edges, m_filtration[i]);
        getAdditionalFacetAtFiltration(facets, m_filtration[i]);
        getAdditionalCellAtFiltration(cells, m_filtration[i]);

        std::cout<<CONSOLE_GREEN<<"Build exterior derivative ..."<<CONSOLE_WHITE<<std::endl;
        std::cout<<"Build ED0T ..."<<std::endl;
        buildED0T(ED0T, edges);
        std::cout<<"Build ED1T ..."<<std::endl;
        buildED1T(ED1T, facets);
        buildED2T(ED2T, cells);

        std::cout<<CONSOLE_GREEN<<"build laplacian ..."<<CONSOLE_WHITE<<std::endl;
        SparseMatrix L0, L1, L2;
        buildLaplacian0(L0, ED0T);
        buildLaplacian1(L1, ED0T, ED1T);
        buildLaplacian2(L2, ED1T, ED2T);

        std::cout<<CONSOLE_GREEN<<"Compute eigenvalues ..."<<CONSOLE_WHITE<<std::endl;
        takeVertexSnapshots(L0);
        takeEdgeSnapshots(L1);
        takeFacetSnapshots(L2);
    }
    */
}

void Snapshot::write(){
    writeVertexSnapshots();
    writeEdgeSnapshots();
    writeFacetSnapshots();
}

void Snapshot::debug(){
    writeVertices();
    
    testEdgeFiltration();
    testFacetFiltration();
    testCellFiltration();


}

void Snapshot::computeMinMaxEdgeLength(){
    m_edge_length_min = std::numeric_limits<double>::max();
    m_edge_length_max = std::numeric_limits<double>::min();
    for(int i=0; i<m_alpha.numEdges(); ++i){
        auto ei = m_alpha.idxToEdge(i);
        m_edge_length_min = std::min(m_edge_length_min, m_alpha.edgeLength(ei));
        m_edge_length_max = std::max(m_edge_length_max, m_alpha.edgeLength(ei));
    }
    //std::cout<<m_edge_length_min<<" "<<m_edge_length_max<<std::endl;
}

void Snapshot::computeRadius(){
    computeEdgeRadius();
    computeFacetRadius();
    computeCellRadius();
}

void Snapshot::buildPriorityQueue(){
    buildEdgePriorityQueue();
    buildFacetPriorityQueue();
    buildCellPrioiryQueue();
}

void Snapshot::assembleSortedElements(){
    assembleSortedEdges();
    assembleSortedFacets();
    assembleSortedCells();
}

void Snapshot::indexSortedElements(){
    indexSortedEdges();
    indexSortedFacets();
    indexSortedCells();
}

void Snapshot::computeEdgeRadius(){
    if (m_complex == 'r'){
        for(int i=0; i<m_alpha.numEdges(); ++i){
            auto ei = m_alpha.idxToEdge(i);
            m_edge_radius[ei] = m_alpha.edgeLength(ei);
        }
    }
    else{
        for(int i=0; i<m_alpha.numEdges(); ++i){
            auto ei = m_alpha.idxToEdge(i);
            Point3 p1 = ei->first->vertex(ei->second)->point();
            Point3 p2 = ei->first->vertex(ei->third)->point();
            Point3 c = CGAL::circumcenter(p1, p2);
            m_edge_radius[ei] = (p1-c).squared_length();
            //std::cout<<m_edge_radius[ei]<<std::endl;
        }
    }
}

void Snapshot::computeFacetRadius(){
    if (m_complex == 'r'){
        for(int i=0; i<m_alpha.numFacets(); ++i){
            auto fi = m_alpha.idxToFacet(i);
            double radius = std::numeric_limits<double>::min();
            
            std::vector<int> idx;
            idx.reserve(3);
            for(int j=0; j<4; ++j){
                if(j==fi->second) continue;
                idx.push_back(j);
            }
            
            for(int j=0; j<idx.size(); ++j){
                int vid1 = idx[j];
                int vid2 = idx[(j+1)%idx.size()];
                if(vid1>vid2) std::swap(vid1, vid2);
                
                auto ei = m_alpha.cellEdgeMap(Edge(fi->first, vid1, vid2));
                
                radius = std::max(radius, m_alpha.edgeLength(ei));
            }
            m_facet_radius[fi] = radius;
        }
    }
    else {
        for(int i=0; i<m_alpha.numFacets(); ++i){
            auto fi = m_alpha.idxToFacet(i);
            std::vector<Point3> vp;
            vp.reserve(3);
            for(int j=0; j<4; ++j){
                if(j==fi->second) continue;
                vp.push_back(fi->first->vertex(j)->point());
            }
            Point3 c = CGAL::circumcenter(vp[0], vp[1], vp[2]);
            m_facet_radius[fi] = (vp[0]-c).squared_length();
        }
    }
}

void Snapshot::computeCellRadius(){
    if (m_complex == 'r'){
        for(int i=0; i<m_alpha.numCells(); ++i){
            auto ci = m_alpha.idxToCell(i);
            double radius = std::numeric_limits<double>::min();

            for(int j=0; j<4; ++j){
                for(int k=j+1; k<4; ++k){
                    auto ei = m_alpha.cellEdgeMap(Edge(ci, j, k));
                    radius = std::max(radius, m_alpha.edgeLength(ei));
                }
            }
            m_cell_radius[ci] = radius;
        }
    }
    else {
        for(int i=0; i<m_alpha.numCells(); ++i){
            auto ci = m_alpha.idxToCell(i);
            std::vector<Point3> vp;
            vp.reserve(4);
            for(int j=0; j<4; ++j){
                vp.push_back(ci->vertex(j)->point());
            }
            Point3 c = CGAL::circumcenter(vp[0], vp[1], vp[2], vp[3]);
            m_cell_radius[ci] = (vp[0]-c).squared_length();
        }
    }
}

void Snapshot::buildEdgePriorityQueue(){
    /*
    EdgeRadiusCompare comp(&m_edge_radius);
    m_edge_heap = new std::priority_queue<EdgeIterator, std::vector<EdgeIterator>, EdgeRadiusCompare>(comp);

    for(int i=0; i<m_alpha.numEdges(); ++i){
        auto ei = m_alpha.idxToEdge(i);
        m_edge_heap->push(ei);
    }
    */
}

void Snapshot::buildFacetPriorityQueue(){
    /*
    FacetRadiusCompare comp(&m_facet_radius);
    m_facet_heap = new std::priority_queue<FacetIterator, std::vector<FacetIterator>, FacetRadiusCompare>(comp);

    for(int i=0; i<m_alpha.numFacets(); ++i){
        auto fi = m_alpha.idxToFacet(i);
        m_facet_heap->push(fi);
    }
    */
}

void Snapshot::buildCellPrioiryQueue(){
    /*
    CellRadiusCompare comp(&m_cell_radius);
    m_cell_heap = new std::priority_queue<CellIterator, std::vector<CellIterator>, CellRadiusCompare>(comp);

    for(int i=0; i<m_alpha.numCells(); ++i){
        auto ci = m_alpha.idxToCell(i);
        m_cell_heap->push(ci);
    }
    */
}

void Snapshot::assembleSortedEdges(){
    m_sorted_edges.reserve(m_alpha.numEdges());
    for(int i=0; i<m_alpha.numEdges(); ++i){
        auto ei = m_alpha.idxToEdge(i);
        m_sorted_edges.push_back(ei);
    }

    std::sort(m_sorted_edges.begin(), m_sorted_edges.end(), EdgeRadiusCompare(&m_alpha));
}

void Snapshot::assembleSortedFacets(){
    m_sorted_facets.reserve(m_alpha.numFacets());
    for(int i=0; i<m_alpha.numFacets(); ++i){
        auto fi = m_alpha.idxToFacet(i);
        m_sorted_facets.push_back(fi);
    }

    std::sort(m_sorted_facets.begin(), m_sorted_facets.end(), FacetRadiusCompare(&m_alpha));
}

void Snapshot::assembleSortedCells(){
    m_sorted_cells.reserve(m_alpha.numCells());
    for(int i=0; i<m_alpha.numCells(); ++i){
        auto ci = m_alpha.idxToCell(i);
        m_sorted_cells.push_back(ci);
    }

    std::sort(m_sorted_cells.begin(), m_sorted_cells.end(), CellRadiusCompare(&m_alpha));
}

void Snapshot::indexSortedEdges(){
    for(int i=0; i<m_sorted_edges.size(); ++i){
        m_sorted_edge_idx[m_sorted_edges[i]] = i;
    }
}

void Snapshot::indexSortedFacets(){
    for(int i=0; i<m_sorted_facets.size(); ++i){
        m_sorted_facet_idx[m_sorted_facets[i]] = i;
    }
}

void Snapshot::indexSortedCells(){
    for(int i=0; i<m_sorted_cells.size(); ++i){
        m_sorted_cell_idx[m_sorted_cells[i]] = i;
    }
}

void Snapshot::writeVertexSnapshots(){
    std::ofstream out("snapshots_vertex.txt");
    //out<<std::scientific;
    for(int i=0; i<m_vertex_snapshots.size(); ++i){
        for(int j=0; j<m_vertex_snapshots[i].size(); ++j){
            out<<m_vertex_snapshots[i][j]<<" ";
        }
        out<<std::endl;
    }
    out.close();
}

void Snapshot::writeEdgeSnapshots(){
    std::ofstream out("snapshots_edge.txt");
    //out<<std::scientific;
    for(int i=0; i<m_edge_snapshots.size(); ++i){
        for(int j=0; j<m_edge_snapshots[i].size(); ++j){
            out<<m_edge_snapshots[i][j]<<" ";
        }
        out<<std::endl;
    }
    out.close();
}

void Snapshot::writeFacetSnapshots(){
    std::ofstream out("snapshots_facet.txt");
    //out<<std::scientific;
    for(int i=0; i<m_facet_snapshots.size(); ++i){
        for(int j=0; j<m_facet_snapshots[i].size(); ++j){
            out<<m_facet_snapshots[i][j]<<" ";
        }
        out<<std::endl;
    }
    out.close();
}

void Snapshot::determineNextEdgeSize(const int& current_size, int& next_size, double filtration){
    next_size = current_size;
    while(next_size < m_sorted_edges.size() && m_alpha.edgeAlpha(m_sorted_edges[next_size]) < filtration){
        ++next_size;
    }
}

void Snapshot::determineNextFacetSize(const int& current_size, int& next_size, double filtration){
    next_size = current_size;
    while(next_size < m_facet_radius.size() && m_alpha.facetAlpha(m_sorted_facets[next_size]) < filtration){
        ++next_size;
    }
}

void Snapshot::determineNextCellSize(const int& current_size, int& next_size, double filtration){
    next_size = current_size;
    while(next_size < m_cell_radius.size() && m_alpha.cellAlpha(m_sorted_cells[next_size]) < filtration){
        ++next_size;
    }
}

/*
void Snapshot::getAdditionalEdgeAtFiltration(std::vector<EdgeIterator>& edges, double filtration){
    while(m_edge_radius[m_edge_heap->top()] < filtration){
        edges.push_back(m_edge_heap->top());
        m_edge_heap->pop();
    }
}

void Snapshot::getAdditionalFacetAtFiltration(std::vector<FacetIterator>& facets, double filtration){
    while(m_facet_radius[m_facet_heap->top()] < filtration){
        facets.push_back(m_facet_heap->top());
        m_facet_heap->pop();
    }
}

void Snapshot::getAdditionalCellAtFiltration(std::vector<CellIterator>& cells, double filtration){
    while(m_cell_radius[m_cell_heap->top()] < filtration){
        cells.push_back(m_cell_heap->top());
        m_cell_heap->pop();
    }
}
*/

void Snapshot::fillED0T(SparseMatrix& ED0T, const int& current_size, const int& next_size){
    for(int i=current_size; i<next_size; ++i){
        // boundary vertex set
        auto bvs = m_alpha.edgeBoundary(m_sorted_edges[i]);
        for(auto it = bvs.begin(); it != bvs.end(); ++it){
            ED0T.insert(m_alpha.vertexIdx(it->first), i) = it->second;
        }
    }
    
    /*
    for(int i=0; i<edges.size(); ++i){
        // boundary vertex set
        auto bvs = m_alpha.edgeBoundary(edges[i]);
        for(auto it = bvs.begin(); it != bvs.end(); ++it){
            ED0T.insert(m_alpha.vertexIdx(it->first), m_alpha.edgeIdx(edges[i])) = it->second;
        }
    }
    */
}

void Snapshot::fillED1T(SparseMatrix& ED1T, const int& current_size, const int& next_size){
    for(int i=current_size; i<next_size; ++i){
        // boundary edge set
        auto bes = m_alpha.facetBoundary(m_sorted_facets[i]);
        for(auto it = bes.begin(); it != bes.end(); ++it){
            ED1T.insert(m_sorted_edge_idx[it->first], i) = it->second;
        }
    }

    /*
    for(int i=0; i<facets.size(); ++i){
        // boundary edge set
        auto bes = m_alpha.facetBoundary(facets[i]);
        for(auto it = bes.begin(); it != bes.end(); ++it){
            ED1T.insert(m_alpha.edgeIdx(it->first), m_alpha.facetIdx(facets[i])) = it->second;
        }
    }
    */
}

void Snapshot::fillED2T(SparseMatrix& ED2T, const int& current_size, const int& next_size){
    for(int i=current_size; i<next_size; ++i){
        // boundary facet set
        auto bfs = m_alpha.cellBoundary(m_sorted_cells[i]);
        for(auto it = bfs.begin(); it != bfs.end(); ++it){
            ED2T.insert(m_sorted_facet_idx[it->first], i) = it->second;
        }
    }

    /*
    for(int i=0; i<cells.size(); ++i){
        // boundary facet set
        auto bfs = m_alpha.cellBoundary(cells[i]);
        for(auto it = bfs.begin(); it != bfs.end(); ++it){
            ED2T.insert(m_alpha.facetIdx(it->first), m_alpha.cellIdx(cells[i])) = it->second;
        }
    }
    */
}

void Snapshot::buildLaplacian0(SparseMatrix& L0, const SparseMatrix& ED0T, const int& edge_size){
    const int vertex_size = m_alpha.numVertices();
    SparseMatrix ED0T_block = ED0T.block(0,0,vertex_size,edge_size);
    L0 = ED0T_block*ED0T_block.transpose();
    //L0 = ED0T*ED0T.transpose();
}

void Snapshot::buildLaplacian1(SparseMatrix& L1, const SparseMatrix& ED0T, const SparseMatrix& ED1T, const int& edge_size, const int& facet_size){
    const int vertex_size = m_alpha.numVertices();
    SparseMatrix ED0T_block = ED0T.block(0,0,vertex_size,edge_size);
    SparseMatrix ED1T_block = ED1T.block(0,0,edge_size,facet_size);
    L1 = ED0T_block.transpose()*ED0T_block + ED1T_block*ED1T_block.transpose();
    //L1 = ED0T.transpose()*ED0T + ED1T*ED1T.transpose();
}

void Snapshot::buildLaplacian2(SparseMatrix& L2, const SparseMatrix& ED1T, const SparseMatrix& ED2T, const int& edge_size, const int& facet_size, const int& cell_size) {
	SparseMatrix ED1T_block = ED1T.block(0, 0, edge_size, facet_size);
	SparseMatrix ED2T_block = ED2T.block(0, 0, facet_size, cell_size);
	L2 = ED1T_block.transpose()*ED1T_block + ED2T_block * ED2T_block.transpose();
	//L2 = ED1T.transpose()*ED1T + ED2T*ED2T.transpose();
}

void Snapshot::buildPersistentLaplacian0(SparseMatrix& L0, const SparseMatrix& ED0T, const SparseMatrix& ED1T, const int& edge_size, const double filtration, double p) {
		p = pow(sqrt(filtration) + p, 2.) - filtration;

		const int vertex_size = m_alpha.numVertices();
		int p_vertex = 0; //always 0 for alpha-complex;

		int next_size = edge_size;

		while (next_size < m_sorted_edges.size() && m_alpha.edgeAlpha(m_sorted_edges[next_size]) < filtration + p) {
			++next_size;
		}
		int p_edge = next_size - edge_size;

		SparseMatrix ED0T_block = ED0T.block(0, 0, vertex_size, edge_size + p_edge);

		//For alpha-complex or Rips-complex, all 0-simplices are present in the first subcomplex of the filtration
		L0 = ED0T_block * ED0T_block.transpose();
}


void Snapshot::buildPersistentLaplacian1(SparseMatrix& L1, const SparseMatrix& ED0T, const SparseMatrix& ED1T, const int& edge_size, const int& facet_size, const double filtration, double p) {
    
    p = pow(sqrt(filtration) + p, 2.) - filtration;
    
    const int vertex_size = m_alpha.numVertices();
    int p_vertex = 0; //always 0 for alpha-complex;
    
    int next_size = edge_size;
    
    while(next_size < m_sorted_edges.size() && m_alpha.edgeAlpha(m_sorted_edges[next_size]) < filtration + p){
        ++next_size;
    }
    int p_edge = next_size - edge_size;
    
    next_size = facet_size;
    while(next_size < m_sorted_facets.size() && m_alpha.facetAlpha(m_sorted_facets[next_size]) < filtration + p){
        ++next_size;
    }
    int p_face = next_size - facet_size;
    
    SparseMatrix ED0T_block = ED0T.block(0, 0, vertex_size, edge_size);
    SparseMatrix ED1T_block = ED1T.block(0, 0, edge_size, facet_size + p_face);
    // size (edge_size, facet_size + p_face) starting at (0, 0)
    
    if (p_edge != 0 && p_face != 0 && p != 0) {
        SparseMatrix L1diff(p_edge, p_edge);
        SparseMatrix Diff = ED1T.block(edge_size, 0, p_edge, facet_size + p_face);
        
        SparseMatrix eye(facet_size + p_face, facet_size + p_face);
        eye.setIdentity();
        
        SparseMatrix gauge = ED0T.block(0, edge_size, vertex_size, p_edge);

		SparseMatrix gaugeTranspose = gauge.transpose();
        
        // get rid of all cols corresponding to a vertex in gaugeTranspose
        // setting dependence on the vertices touched by alpha complex to 0
        for (int k = 0; k < edge_size; ++k) {
            for (SparseMatrix::InnerIterator it(ED0T, k); it; ++it) {
                if (fabs(it.value()) > 0.5 && it.index()< vertex_size) {
                    for (SparseMatrix::InnerIterator it2(gaugeTranspose, it.index()); it2; ++it2) {
                        it2.valueRef()=0.;
                    }
                }
            }
        }
        ////L1diff = gaugeTranspose * gauge + Diff * Diff.transpose(); //Need to have deficiency fixing if the different complex has more than 1 boundary pieces properly
        SparseMatrix eye1(p_edge, p_edge);
        //SparseMatrix eye1(edge_size + p_edge, edge_size + p_edge);
        eye1.setIdentity();
        L1diff = gaugeTranspose*gauge + Diff*Diff.transpose() + 1e-12*eye1; //This is just adding a bit of damping
        Eigen::ConjugateGradient<SparseMatrix> solver;
        solver.compute(L1diff);
        SparseMatrix I(p_edge, p_edge);
        //SparseMatrix I(edge_size + p_edge, edge_size + p_edge);
        I.setIdentity();
        SparseMatrix L1diff_inv = solver.solve(I);
        ////L1 = ED0T_block.transpose()*ED0T_block + ED1T_block*(eye-Diff.transpose()*L1diff_inv*Diff)* ED1T_block.transpose();
        SparseMatrix L1_t1 = L1diff*L1diff_inv;
        SparseMatrix L1_t = L1diff_inv*L1_t1;
        SparseMatrix L1_temp = eye - Diff.transpose()*L1_t*Diff;
        SparseMatrix L1_temp1 = ED1T_block*L1_temp*ED1T_block.transpose();
        //Eigen::SparseMatrix<double> L1_temp1 = ED1T_block*ED1T_block.transpose();
        L1  = ED0T_block.transpose()*ED0T_block + L1_temp1;
        //L1 = ED0T_block.transpose()*ED0T_block + ED1T_block*(eye-Diff.transpose()*(L1diff_inv*L1diff*L1diff_inv)*Diff)* ED1T_block.transpose();

    } else {
        L1 = ED0T_block.transpose()*ED0T_block + ED1T_block*ED1T_block.transpose();
    }
}

void Snapshot::buildPersistentLaplacian2(SparseMatrix& L2, const SparseMatrix& ED1T, const SparseMatrix& ED2T, const int& edge_size, const int& facet_size, const int & cell_size, const double filtration, double p) {

	p = pow(sqrt(filtration) + p, 2.) - filtration;

	int next_size = edge_size;

	while (next_size < m_sorted_edges.size() && m_alpha.edgeAlpha(m_sorted_edges[next_size]) < filtration + p) {
		++next_size;
	}
	int p_edge = next_size - edge_size;

	next_size = facet_size;
	while (next_size < m_sorted_facets.size() && m_alpha.facetAlpha(m_sorted_facets[next_size]) < filtration + p) {
		++next_size;
	}
	int p_face = next_size - facet_size;

	next_size = cell_size;
	while (next_size < m_sorted_cells.size() && m_alpha.cellAlpha(m_sorted_cells[next_size]) < filtration + p) {
		++next_size;
	}
	int p_cell = next_size - cell_size;

	SparseMatrix ED1T_block = ED1T.block(0, 0, edge_size, facet_size);
	SparseMatrix ED2T_block = ED2T.block(0, 0, facet_size, cell_size + p_cell);

	if (p_face != 0 && p_cell != 0 && p != 0) {
		SparseMatrix L2diff(p_face, p_face);
		SparseMatrix Diff = ED2T.block(facet_size, 0, p_face, cell_size + p_cell);

		SparseMatrix eye(cell_size + p_cell, cell_size + p_cell);
		eye.setIdentity();

		SparseMatrix gauge = ED1T.block(0, facet_size, edge_size, p_face);

		SparseMatrix gaugeTranspose = gauge.transpose();

		// get rid of all cols corresponding to a vertex in gaugeTranspose
		// setting dependence on the vertices touched by alpha complex to 0
		for (int k = 0; k < facet_size; ++k) {
			for (SparseMatrix::InnerIterator it(ED1T, k); it; ++it) {
				if (fabs(it.value()) > 0.5 && it.index() < edge_size) {
					for (SparseMatrix::InnerIterator it2(gaugeTranspose, it.index()); it2; ++it2) {
						it2.valueRef() = 0.;
					}
				}
			}
		}
		////L2diff = gaugeTranspose * gauge + Diff * Diff.transpose(); //Need to have deficiency fixing if the different complex has more than 1 boundary pieces properly
		SparseMatrix eye1(p_face, p_face);
		eye1.setIdentity();
		
		L2diff = gaugeTranspose * gauge + Diff * Diff.transpose() + 1e-12*eye1; //This is just adding a bit of damping
		Eigen::ConjugateGradient<SparseMatrix> solver;
		solver.compute(L2diff);
		SparseMatrix I(p_face, p_face);
		I.setIdentity();
		SparseMatrix L2diff_inv = solver.solve(I);
		////L1 = ED0T_block.transpose()*ED0T_block + ED1T_block*(eye-Diff.transpose()*L1diff_inv*Diff)* ED1T_block.transpose();
		SparseMatrix L2_t1 = L2diff * L2diff_inv;
		SparseMatrix L2_t = L2diff_inv * L2_t1;
		SparseMatrix L2_temp = eye - Diff.transpose()*L2_t*Diff;
		SparseMatrix L2_temp1 = ED2T_block * L2_temp*ED2T_block.transpose();

		L2 = ED1T_block.transpose()*ED1T_block + L2_temp1;
		//L2 = ED1T_block.transpose()*ED1T_block + ED2T_block*(eye-Diff.transpose()*(L2diff_inv.transpose()*L2diff*L2diff_inv)*Diff)* ED2T_block.transpose();

	}
	else {
		L2 = ED1T_block.transpose()*ED1T_block + ED2T_block * ED2T_block.transpose();
	}
}

void Snapshot::takeVertexSnapshots(const SparseMatrix& L0, const int matrix_size){
    int size = L0.nonZeros();
    std::vector<double> row, col, val;
    row.reserve(size);
    col.reserve(size);
    val.reserve(size);

    for(int i=0; i<L0.outerSize(); ++i){
        for(SparseMatrix::InnerIterator iter(L0, i); iter; ++iter){
            row.push_back(static_cast<double>(iter.row()+1));
            col.push_back(static_cast<double>(iter.col()+1));
            val.push_back(static_cast<double>(iter.value()));
        }
    }

    std::vector<double> eigenvalues;
    std::vector<ColumnVector> eigenvectors;
    
    std::cout<<"Matrix size ["<< L0.innerSize() <<", " << L0.outerSize() <<"], ";
    auto time_start = Clock::now();
    matlabEIGS(eigenvalues, eigenvectors, m_alpha.numVertices(), row, col, val, m_num_eigenvalues);
    auto time_end = Clock::now();
    auto secs = std::chrono::duration_cast<TimeSecond>(time_end-time_start);
    std::cout<<"Elapsed time: "<<secs.count()<<std::endl;

    if(eigenvalues.empty()){
        int n = std::min(m_alpha.numVertices(), m_num_eigenvalues);
        eigenvalues.resize(n, 0);
    }
    m_vertex_snapshots.push_back(eigenvalues);
    
}

void Snapshot::takeEdgeSnapshots(const SparseMatrix& L1, const int matrix_size){
    int size = L1.nonZeros();
    std::vector<double> row, col, val;
    row.reserve(size);
    col.reserve(size);
    val.reserve(size);

    for(int i=0; i<L1.outerSize(); ++i){
        for(SparseMatrix::InnerIterator iter(L1, i); iter; ++iter){
            row.push_back(static_cast<double>(iter.row()+1));
            col.push_back(static_cast<double>(iter.col()+1));
            val.push_back(static_cast<double>(iter.value()));
        }
    }

    std::vector<double> eigenvalues;
    std::vector<ColumnVector> eigenvectors;
    
    std::cout<<"Matrix size ["<< L1.innerSize() <<", " << L1.outerSize() <<"], ";
    auto time_start = Clock::now();
    matlabEIGS(eigenvalues, eigenvectors, matrix_size, row, col, val, m_num_eigenvalues);
    auto time_end = Clock::now();
    auto secs = std::chrono::duration_cast<TimeSecond>(time_end-time_start);
    std::cout<<"Elapsed time: "<<secs.count()<<std::endl;

    m_edge_snapshots.push_back(eigenvalues);
}

void Snapshot::takeFacetSnapshots(const SparseMatrix& L2, const int matrix_size){
    int size = L2.nonZeros();
    std::vector<double> row, col, val;
    row.reserve(size);
    col.reserve(size);
    val.reserve(size);

    for(int i=0; i<L2.outerSize(); ++i){
        for(SparseMatrix::InnerIterator iter(L2, i); iter; ++iter){
            row.push_back(static_cast<double>(iter.row()+1));
            col.push_back(static_cast<double>(iter.col()+1));
            val.push_back(static_cast<double>(iter.value()));
        }
    }

    std::vector<double> eigenvalues;
    std::vector<ColumnVector> eigenvectors;
    
    std::cout<<"Matrix size ["<< L2.innerSize() <<", " << L2.outerSize() <<"], ";
    auto time_start = Clock::now();
    matlabEIGS(eigenvalues, eigenvectors, matrix_size, row, col, val, m_num_eigenvalues);
    auto time_end = Clock::now();
    auto secs = std::chrono::duration_cast<TimeSecond>(time_end-time_start);
    std::cout<<"Elapsed time: "<<secs.count()<<std::endl;
    
    m_facet_snapshots.push_back(eigenvalues);
}

void Snapshot::startMatlab(){
    m_matlab_engine = matlab::engine::startMATLAB();
}

void Snapshot::matlabEIGS(std::vector<double>& eval, std::vector<ColumnVector>& evec, int matrix_size,
    std::vector<double>& row, std::vector<double>& col, std::vector<double>& val, int es){
    std::vector<double> vms = {static_cast<double>(matrix_size)};
    std::vector<double> ves = {static_cast<double>(es)};
    using namespace matlab;

    data::ArrayFactory factory;
    
    data::TypedArray<double> mms = factory.createArray<std::vector<double>::iterator>({ 1, 1 }, vms.begin(), vms.end());
    data::TypedArray<double> mrow = factory.createArray<std::vector<double>::iterator>({ row.size(), 1 }, row.begin(), row.end());
    data::TypedArray<double> mcol = factory.createArray<std::vector<double>::iterator>({ col.size(), 1 }, col.begin(), col.end());
    data::TypedArray<double> mval = factory.createArray<std::vector<double>::iterator>({ val.size(), 1 }, val.begin(), val.end());
    data::TypedArray<double> mes = factory.createArray<std::vector<double>::iterator>({ 1, 1 }, ves.begin(), ves.end());
    
    m_matlab_engine->setVariable(u"size", std::move(mms));
    m_matlab_engine->setVariable(u"row", std::move(mrow));
    m_matlab_engine->setVariable(u"col", std::move(mcol));
    m_matlab_engine->setVariable(u"val", std::move(mval));
    m_matlab_engine->setVariable(u"es", std::move(mes));
    m_matlab_engine->eval(u"A=sparse(row, col, val, size, size);");
    m_matlab_engine->eval(u"A=A+speye(size)*1e-12;");
    m_matlab_engine->eval(u"if size > es \n num = es; \n else \n num = size; \n end \n");
    m_matlab_engine->eval(u"[V,D]=eigs((A+A')/2, num, 'smallestabs');");
    data::TypedArray<double> dd = m_matlab_engine->getVariable(u"D");
    data::TypedArray<double> vv = m_matlab_engine->getVariable(u"V");

    eval.reserve(dd.getDimensions()[0]);
    for(int i=0; i<dd.getDimensions()[0]; ++i){
        eval.push_back(dd[i][i]);
    }
}

void Snapshot::writeVertices(){
    std::ofstream out("triangulation_vertices.vtk");
    
    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Triangulation Vertices"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    out<<"POINTS "<<m_alpha.numVertices()<<" double"<<std::endl;
    for(int i=0; i<m_alpha.numVertices(); ++i){
        auto vi = m_alpha.idxToVertex(i);
        out<<vi->point().x()<<" "<<vi->point().y()<<" "<<vi->point().z()<<std::endl;
    }

    out<<"CELLS "<<m_alpha.numVertices()<<" "<<2*m_alpha.numVertices()<<std::endl;
    for(int i=0; i<m_alpha.numVertices(); ++i){
        out<<"1 "<<i<<std::endl;
    }

    out<<"CELL_TYPES "<<m_alpha.numVertices()<<std::endl;
    for(size_t i=0; i<m_alpha.numVertices(); ++i){
        out<<"1"<<std::endl;
    }
}

void Snapshot::testEdgeFiltration(){
    // get edges under a threshold
    double threshold = 10000;
    std::vector<EdgeIterator> edges;
    /*
    while(!m_edge_heap->empty() && m_edge_radius[m_edge_heap->top()]<threshold){
        edges.push_back(m_edge_heap->top());
        m_edge_heap->pop();
    }
    */
    for(int i=0; i<m_sorted_edges.size(); ++i){
        if(m_edge_radius[m_sorted_edges[i]] > threshold) break;
        edges.push_back(m_sorted_edges[i]);
    }

    std::ofstream out("triangulation_edges.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Triangulation Edges"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    out<<"POINTS "<<m_alpha.numVertices()<<" double"<<std::endl;
    for(int i=0; i<m_alpha.numVertices(); ++i){
        auto vi = m_alpha.idxToVertex(i);
        out<<vi->point().x()<<" "<<vi->point().y()<<" "<<vi->point().z()<<std::endl;
    }

    out<<"CELLS "<<edges.size()<<" "<<3*edges.size()<<std::endl;
    for(int i=0; i<edges.size(); ++i){
        const auto& ei = edges[i];
        auto vi1 = ei->first->vertex(ei->second);
        auto vi2 = ei->first->vertex(ei->third);
        out<<"2 "<<m_alpha.vertexIdx(vi1)<<" "<<m_alpha.vertexIdx(vi2)<<std::endl;
    }

    out<<"CELL_TYPES "<<edges.size()<<std::endl;
    for(size_t i=0; i<edges.size(); ++i){
        out<<"3"<<std::endl;
    }

    out.close();
}

void Snapshot::testFacetFiltration(){
    // get facets under a threshold
    double threshold = 10000;
    std::vector<FacetIterator> facets;
    /*
    while(!m_facet_heap->empty() && m_facet_radius[m_facet_heap->top()]<threshold){
        facets.push_back(m_facet_heap->top());
        m_facet_heap->pop();
    }
    */
    for(int i=0; i<m_sorted_facets.size(); ++i){
        if(m_facet_radius[m_sorted_facets[i]] > threshold) break;
        facets.push_back(m_sorted_facets[i]);
    }

    std::ofstream out("triangulation_facets.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Triangulation Facets"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    out<<"POINTS "<<m_alpha.numVertices()<<" double"<<std::endl;
    for(int i=0; i<m_alpha.numVertices(); ++i){
        auto vi = m_alpha.idxToVertex(i);
        out<<vi->point().x()<<" "<<vi->point().y()<<" "<<vi->point().z()<<std::endl;
    }

    out<<"CELLS "<<facets.size()<<" "<<4*facets.size()<<std::endl;
    for(int i=0; i<facets.size(); ++i){
        const auto& fi = facets[i];
        out<<"3 ";
        for(int j=0; j<4; ++j){
            if(j==fi->second) continue;
            out<<m_alpha.vertexIdx(fi->first->vertex(j))<<" ";
        }
        out<<std::endl;
    }

    out<<"CELL_TYPES "<<facets.size()<<std::endl;
    for(size_t i=0; i<facets.size(); ++i){
        out<<"5"<<std::endl;
    }

    out.close();
}

void Snapshot::testCellFiltration(){
    // get cells under a threshold
    double threshold = 10000;
    std::vector<CellIterator> cells;
    /*
    while(!m_cell_heap->empty() && m_cell_radius[m_cell_heap->top()] < threshold){
        cells.push_back(m_cell_heap->top());
        m_cell_heap->pop();
    }
    */
    for(int i=0; i<m_sorted_cells.size(); ++i){
        if(m_cell_radius[m_sorted_cells[i]] > threshold) break;
        cells.push_back(m_sorted_cells[i]);
    }

    std::ofstream out("triangulation_cells.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Triangulation Cells"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    out<<"POINTS "<<m_alpha.numVertices()<<" double"<<std::endl;
    for(int i=0; i<m_alpha.numVertices(); ++i){
        auto vi = m_alpha.idxToVertex(i);
        out<<vi->point().x()<<" "<<vi->point().y()<<" "<<vi->point().z()<<std::endl;
    }

    out<<"CELLS "<<cells.size()<<" "<<5*cells.size()<<std::endl;
    for(int i=0; i<cells.size(); ++i){
        const auto&ci = cells[i];
        out<<"4 ";
        for(int j=0; j<4; ++j){
            out<<m_alpha.vertexIdx(ci->vertex(j))<<" ";
        }
        out<<std::endl;
    }

    out<<"CELL_TYPES "<<cells.size()<<std::endl;
    for(size_t i=0; i<cells.size(); ++i){
        out<<"10"<<std::endl;
    }

    out.close();
}
