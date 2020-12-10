#include "alpha_shape.hpp"

AlphaShape::AlphaShape(){

}

AlphaShape::~AlphaShape(){

}

bool AlphaShape::validVertex(const VertexIterator& vi) const{
    return m_valid_vertex.at(vi);
}

bool AlphaShape::validEdge(const EdgeIterator& ei) const{
    return m_valid_edge.at(ei);
}

bool AlphaShape::validFacet(const FacetIterator& fi) const{
    return m_valid_facet.at(fi);
}

bool AlphaShape::validCell(const CellIterator& ci) const{
    return m_valid_cell.at(ci);
}

int AlphaShape::vertexIdx(const VertexIterator& vi) const{
    return m_vertex_idx.at(vi);
}

int AlphaShape::edgeIdx(const EdgeIterator& ei) const{
    return m_edge_idx.at(ei);
}

int AlphaShape::facetIdx(const FacetIterator& fi) const{
    return m_facet_idx.at(fi);
}

int AlphaShape::cellIdx(const CellIterator& ci) const{
    return m_cell_idx.at(ci);
}

VertexIterator AlphaShape::idxToVertex(const int& idx) const{
    return m_idx_to_vertex.at(idx);
}

EdgeIterator AlphaShape::idxToEdge(const int& idx) const{
    return m_idx_to_edge.at(idx);
}

FacetIterator AlphaShape::idxToFacet(const int& idx) const{
    return m_idx_to_facet.at(idx);
}

CellIterator AlphaShape::idxToCell(const int& idx) const{
    return m_idx_to_cell.at(idx);
}

int AlphaShape::numVertices() const{
    return m_num_vertices;
}

int AlphaShape::numEdges() const{
    return m_num_edges;
}

int AlphaShape::numFacets() const{
    return m_num_facets;
}

int AlphaShape::numCells() const{
    return m_num_cells;
}

EdgeIterator AlphaShape::cellEdgeMap(const Edge& e) const{
    return m_cell_edge_map.at(e);
}

FacetIterator AlphaShape::cellFacetMap(const Facet& f) const{
    return m_cell_facet_map.at(f);
}

double AlphaShape::edgeLength(const EdgeIterator& ei) const{
    return m_edge_length.at(ei);
}

double AlphaShape::edgeAlpha(const EdgeIterator& ei) const{
    return m_edge_alpha.at(ei);
}

double AlphaShape::facetAlpha(const FacetIterator& fi) const{
    return m_facet_alpha.at(fi);
}

double AlphaShape::cellAlpha(const CellIterator& ci) const{
    return m_cell_alpha.at(ci);
}

BoundaryVertexSet AlphaShape::edgeBoundary(const EdgeIterator& ei) const{
    return m_edge_boundary.at(ei);
}

BoundaryEdgeSet AlphaShape::facetBoundary(const FacetIterator& fi) const{
    return m_facet_boundary.at(fi);
}

BoundaryFacetSet AlphaShape::cellBoundary(const CellIterator& ci) const{
    return m_cell_boundary.at(ci);
}

void AlphaShape::readPointCloud(std::string filename){
    std::ifstream file(filename);
    
    // add perturbation
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    std::string line;
    while(std::getline(file, line)){
        std::stringstream ss(line);
        double x, y, z;
        ss>>x>>y>>z;
        double x_pert = 0;
        double y_pert = 0;
        double z_pert = 0;
        //double x_pert = distribution(generator)*1e-15;
        //double y_pert = distribution(generator)*1e-15;
        //double z_pert = distribution(generator)*1e-15;
        //std::cout<<x+x_pert<<" "<<y+y_pert<<" "<<z+z_pert<<std::endl;
        m_point_cloud.push_back(Point3(x+x_pert,y+y_pert,z+z_pert));
    }
    file.close();
}

void AlphaShape::buildTriangulation(){
    m_triangulation = Delaunay(m_point_cloud.begin(), m_point_cloud.end());
    
    //std::ofstream out("output");
    //out<<m_triangulation;
    //out.close();
}   

void AlphaShape::preprocess(){
    classifyElements();
    indexElements();
    buildCellToEdgeFacetMap();
    computeEdgeLength();
    computeCircum();
    computeAlpha();
    computeBoundary();
}

void AlphaShape::debug(){
    writeVertices();
    writeEdges();
    writeFacets();
    writeCells();

    writeSortedAlpha();
}

void AlphaShape::classifyElements(){
    classifyVertices();
    classifyEdges();
    classifyFacets();
    classifyCells();
}

void AlphaShape::indexElements(){
    indexVertices();
    indexEdges();
    indexFacets();
    indexCells();

    std::cout<<"Vertex Num: "<<m_num_vertices<<std::endl;
    std::cout<<"Edge Num: "<<m_num_edges<<std::endl;
    std::cout<<"Facet Num: "<<m_num_facets<<std::endl;
    std::cout<<"Cell Num: "<<m_num_cells<<std::endl;
}

void AlphaShape::buildCellToEdgeFacetMap(){
    buildCellToEdgeMap();
    buildCellToFacetMap();
}

void AlphaShape::computeEdgeLength(){
    for(auto ei = m_triangulation.edges_begin(); ei != m_triangulation.edges_end(); ++ei){
        if(!validEdge(ei)) continue;
        VertexIterator vi1 = ei->first->vertex(ei->second);
        VertexIterator vi2 = ei->first->vertex(ei->third);
        m_edge_length[ei] = sqrt((vi1->point()-vi2->point()).squared_length());
    }
}

void AlphaShape::computeCircum(){
    computeEdgeCircum();
    computeFacetCircum();
    computeCellCircum();
}

void AlphaShape::computeAlpha(){
    /*
        Simplex processing order must obey: Cell->Facet->Edge
        Because there is alpha value propagation
    */

    initAlpha();
    processCellAlpha();
    processFacetAlpha();
    processEdgeAlpha();
}

void AlphaShape::computeBoundary(){
    /*
        Natural orientation: increasing indices
            - Edge: direction from smaller index to larger index
            - Face: following edge direction, using right hand rule
    */

    computeEdgeBoundary();
    computeFacetBoundary();
    computeCellBoundary();
}

void AlphaShape::classifyVertices(){
    for(auto vi = m_triangulation.vertices_begin(); vi != m_triangulation.vertices_end(); ++vi){
        m_valid_vertex[vi] =  !m_triangulation.is_infinite(vi);
    }
}

void AlphaShape::classifyEdges(){
    for(auto ei = m_triangulation.edges_begin(); ei != m_triangulation.edges_end(); ++ei){
        VertexIterator vi1 = ei->first->vertex(ei->second);
        VertexIterator vi2 = ei->first->vertex(ei->third);
        m_valid_edge[ei] = validVertex(vi1) && validVertex(vi2);
    }
}

void AlphaShape::classifyFacets(){
    for(auto fi = m_triangulation.facets_begin(); fi != m_triangulation.facets_end(); ++fi){
        bool ans = true;
        for(int i=0; i<4; ++i){
            if(i == fi->second) continue;
            ans = ans && validVertex(fi->first->vertex(i));
        }
        m_valid_facet[fi] = ans;
    }
}

void AlphaShape::classifyCells(){
    for(auto ci = m_triangulation.cells_begin(); ci != m_triangulation.cells_end(); ++ci){
        bool ans = true;
        for(int i=0; i<4; ++i){
            ans = ans && validVertex(ci->vertex(i));
        }
        m_valid_cell[ci] = ans;
    }
}

void AlphaShape::indexVertices(){
    int k = 0;
    for(auto vi = m_triangulation.vertices_begin(); vi != m_triangulation.vertices_end(); ++vi){
        if(!validVertex(vi)) continue;
        m_vertex_idx[vi] = k;
        m_idx_to_vertex[k] = vi;
        ++k;
    }
    m_num_vertices = k;
}

void AlphaShape::indexEdges(){
    int k = 0;
    for(auto ei = m_triangulation.edges_begin(); ei != m_triangulation.edges_end(); ++ei){
        //std::cout<<validEdge(ei)<<std::endl;
        if(!validEdge(ei)) continue;
        m_edge_idx[ei] = k;
        m_idx_to_edge[k] = ei;
        ++k;
    }
    m_num_edges = k;
}

void AlphaShape::indexFacets(){
    int k = 0;
    for(auto fi = m_triangulation.facets_begin(); fi != m_triangulation.facets_end(); ++fi){
        if(!validFacet(fi)) continue;
        m_facet_idx[fi] = k;
        m_idx_to_facet[k] = fi;
        ++k;
    }
    m_num_facets = k;
}

void AlphaShape::indexCells(){
    int k = 0;
    for(auto ci = m_triangulation.cells_begin(); ci != m_triangulation.cells_end(); ++ci){
        if(!validCell(ci)) continue;
        m_cell_idx[ci] = k;
        m_idx_to_cell[k] = ci;
        ++k;
    }
    m_num_cells = k;
}

void AlphaShape::buildCellToEdgeMap(){
    for(auto ei = m_triangulation.edges_begin(); ei != m_triangulation.edges_end(); ++ei){
        if(!validEdge(ei)) continue;

        VertexIterator vi1 = ei->first->vertex(ei->second);
        VertexIterator vi2 = ei->first->vertex(ei->third);
        CellCirculator cc = m_triangulation.incident_cells(*ei);
        CellCirculator end = cc;
        
        do{
            // We should consider invalid cell
            // Because the representative cell of the edge maybe invalid cell
            int i,j;
            cc->has_vertex(vi1, i);
            cc->has_vertex(vi2, j);
            if(i>j) std::swap(i,j);

            Edge e(cc, i, j);
            m_cell_edge_map[e] = ei;
            ++cc;
        }while(cc != end);
    }
}

void AlphaShape::buildCellToFacetMap(){
    for(auto fi = m_triangulation.facets_begin(); fi != m_triangulation.facets_end(); ++fi){
        if(!validFacet(fi)) continue;

        // We should consider invalid cell
        // Because the representative cell of the facet maybe invalid cell
        Facet f(fi->first, fi->second);
        m_cell_facet_map[f] = fi;

        Facet mirror_facet = m_triangulation.mirror_facet(*fi);
        Facet of(mirror_facet.first, mirror_facet.second);
        m_cell_facet_map[of] = fi;
    }
}

void AlphaShape::computeEdgeCircum(){
    for(int i=0; i<m_num_edges; ++i){
        auto ei = m_idx_to_edge[i];
        Point3 p1 = ei->first->vertex(ei->second)->point();
        Point3 p2 = ei->first->vertex(ei->third)->point();
        Point3 c = CGAL::circumcenter(p1, p2);
        m_edge_circumcenter[ei] = c;
        m_edge_radius[ei] = std::sqrt((p1-c).squared_length());
    }
}

void AlphaShape::computeFacetCircum(){
    for(int i=0; i<m_num_facets; ++i){
        auto fi = m_idx_to_facet[i];
        std::vector<Point3> vp;
        vp.reserve(3);
        for(int j=0; j<4; ++j){
            if(j==fi->second) continue;
            vp.push_back(fi->first->vertex(j)->point());
        }
        Point3 c = CGAL::circumcenter(vp[0], vp[1], vp[2]);
        m_facet_circumcenter[fi] = c;
        m_facet_radius[fi] = std::sqrt((vp[0]-c).squared_length());
    }
}

void AlphaShape::computeCellCircum(){
    for(int i=0; i<m_num_cells; ++i){
        auto ci = m_idx_to_cell[i];
        std::vector<Point3> vp;
        vp.reserve(4);
        for(int j=0; j<4; ++j){
            vp.push_back(ci->vertex(j)->point());
        }
        Point3 c = CGAL::circumcenter(vp[0], vp[1], vp[2], vp[3]);
        m_cell_circumcenter[ci] = c;
        m_cell_radius[ci] = std::sqrt((vp[0]-c).squared_length());
    }
}

void AlphaShape::initAlpha(){
    // init cell/facet/edge alpha value to positive infinity
    for(int i=0; i<m_num_cells; ++i){
        auto ci = m_idx_to_cell[i];
        m_cell_alpha[ci] = std::numeric_limits<double>::max();
    }

    for(int i=0; i<m_num_facets; ++i){
        auto fi = m_idx_to_facet[i];
        m_facet_alpha[fi] = std::numeric_limits<double>::max();
    }

    for(int i=0; i<m_num_edges; ++i){
        auto ei = m_idx_to_edge[i];
        m_edge_alpha[ei] = std::numeric_limits<double>::max();
    }
}

void AlphaShape::processCellAlpha(){
    for(int i=0; i<m_num_cells; ++i){
        auto ci = m_idx_to_cell[i];
        if(m_cell_alpha[ci] != std::numeric_limits<double>::max()) continue;
        m_cell_alpha[ci] = std::pow(m_cell_radius[ci],2);
    }

    // propagation to facet
    for(int i=0; i<m_num_cells; ++i){
        auto ci = m_idx_to_cell[i];
        for(int j=0; j<4; ++j){
            Facet f(ci, j);
            auto fi = m_cell_facet_map[f];
            if(m_facet_alpha[fi] == std::numeric_limits<double>::max()){
                Vector3 d = ci->vertex(j)->point() - m_facet_circumcenter[fi];
                double r = std::sqrt(d.squared_length());
                if(r>m_facet_radius[fi]) continue;
                m_facet_alpha[fi] = m_cell_alpha[ci];
            }else{
                m_facet_alpha[fi] = std::min(m_facet_alpha[fi], m_cell_alpha[ci]);
            }
        }
    }
}

void AlphaShape::processFacetAlpha(){
    for(int i=0; i<m_num_facets; ++i){
        auto fi = m_idx_to_facet[i];
        if(m_facet_alpha[fi] != std::numeric_limits<double>::max()) continue;
        m_facet_alpha[fi] = std::pow(m_facet_radius[fi],2);
    }

    // propagation to edge
    for(int i=0; i<m_num_facets; ++i){
        auto fi = m_idx_to_facet[i];
        std::vector<int> lid; // local idx for vertices in this facet
        lid.reserve(3);
        for(int j=0; j<4; ++j){
            if(j == fi->second) continue;
            lid.push_back(j);
        }

        for(int j=0; j<lid.size(); ++j){
            int vid1 = lid[j];
            int vid2 = lid[(j+1)%lid.size()];
            if(vid1>vid2) std::swap(vid1, vid2);
            Edge e(fi->first, vid1, vid2);
            auto ei = m_cell_edge_map[e];
            if(m_edge_alpha[ei] == std::numeric_limits<double>::max()){
                Vector3 d = fi->first->vertex(lid[(j+2)%lid.size()])->point() - m_edge_circumcenter[ei];
                double r = std::sqrt(d.squared_length());
                if(r>m_edge_radius[ei]) continue;
                m_edge_alpha[ei] = m_facet_alpha[fi];
            }else{
                m_edge_alpha[ei] = std::min(m_edge_alpha[ei], m_facet_alpha[fi]);
            }
        }
    }
}

void AlphaShape::processEdgeAlpha(){
    for(int i=0; i<m_num_edges; ++i){
        auto ei = m_idx_to_edge[i];
        if(m_edge_alpha[ei] != std::numeric_limits<double>::max()) continue;
        m_edge_alpha[ei] = std::pow(m_edge_radius[ei],2);
    }
}

void AlphaShape::computeEdgeBoundary(){
    for(auto ei = m_triangulation.edges_begin(); ei != m_triangulation.edges_end(); ++ei){
        if(!validEdge(ei)) continue;
        VertexIterator vi1 = ei->first->vertex(ei->second);
        VertexIterator vi2 = ei->first->vertex(ei->third);
        if(vertexIdx(vi1) > vertexIdx(vi2)) std::swap(vi1, vi2);
        m_edge_boundary[ei].insert(VertexAsBoundary(vi1, -1));
        m_edge_boundary[ei].insert(VertexAsBoundary(vi2, 1));
    }
}

void AlphaShape::computeFacetBoundary(){
    for(auto fi = m_triangulation.facets_begin(); fi != m_triangulation.facets_end(); ++fi){
        if(!validFacet(fi)) continue;
        typedef std::pair<int,int> Pair; // global index, local index
        std::vector<Pair> idx;
        idx.reserve(3);
        for(int i=0; i<4; ++i){
            if(i==fi->second) continue;
            idx.push_back(Pair(vertexIdx(fi->first->vertex(i)), i));
        }
        
        std::sort(idx.begin(), idx.end(), [](const Pair& p1, const Pair& p2) -> bool{
            return p1.first < p2.first;
        });
        
        for(int i=0; i<idx.size(); ++i){
            int lid1 = idx[i].second; // local index
            int lid2 = idx[(i+1)%idx.size()].second; // local index
            if(lid1>lid2) std::swap(lid1, lid2);
            int val = i==idx.size()-1?-1:1;
            m_facet_boundary[fi].insert(EdgeAsBoundary(cellEdgeMap(Edge(fi->first, lid1, lid2)),val));
        }
    }
}

void AlphaShape::computeCellBoundary(){
    // TODO: need to make sure the orientation is correct
    //       i.e. it really computes the divergence (total outward flux)
    for(auto ci = m_triangulation.cells_begin(); ci != m_triangulation.cells_end(); ++ci){
        if(!validCell(ci)) continue;
        for(int i=0; i<4; ++i){
            int even;
            if(i%2==0) even = 1;
            else even = -1;

            std::vector<int> idx;
            idx.reserve(3);
            for(int j=0; j<4; ++j){
                if(j==i) continue;
                idx.push_back(vertexIdx(ci->vertex(j)));
            }

            for(int j=0; j<idx.size(); ++j){
                for(int k=j+1; k<idx.size(); ++k){
                    if(idx[j] < idx[k]) continue;
                    std::swap(idx[j], idx[k]);
                    even = -even;
                }
            }

            m_cell_boundary[ci].insert(FacetAsBoundary(cellFacetMap(Facet(ci, i)),even));
        }
    }
}



// debug
void AlphaShape::writeVertices(){
    std::ofstream out("triangulation_vertices.vtk");
    
    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Triangulation Vertices"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    out<<"POINTS "<<m_num_vertices<<" double"<<std::endl;
    for(auto vi = m_triangulation.vertices_begin(); vi != m_triangulation.vertices_end(); ++vi){
        if(!validVertex(vi)) continue;
        out<<vi->point().x()<<" "<<vi->point().y()<<" "<<vi->point().z()<<std::endl;
    }

    out<<"CELLS "<<m_num_vertices<<" "<<2*m_num_vertices<<std::endl;
    for(int i=0; i<m_num_vertices; ++i){
        out<<"1 "<<i<<std::endl;
    }

    out<<"CELL_TYPES "<<m_num_vertices<<std::endl;
    for(size_t i=0; i<m_num_vertices; ++i){
        out<<"1"<<std::endl;
    }
}

void AlphaShape::writeEdges(){
    std::ofstream out("triangulation_edges.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Triangulation Edges"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    out<<"POINTS "<<m_num_vertices<<" double"<<std::endl;
    for(auto vi = m_triangulation.vertices_begin(); vi != m_triangulation.vertices_end(); ++vi){
        if(!validVertex(vi)) continue;
        out<<vi->point().x()<<" "<<vi->point().y()<<" "<<vi->point().z()<<std::endl;
    }

    out<<"CELLS "<<m_num_edges<<" "<<3*m_num_edges<<std::endl;
    for(auto ei = m_triangulation.edges_begin(); ei != m_triangulation.edges_end(); ++ei){
        if(!validEdge(ei)) continue;
        VertexIterator vi1 = ei->first->vertex(ei->second);
        VertexIterator vi2 = ei->first->vertex(ei->third);
        out<<"2 "<<vertexIdx(vi1)<<" "<<vertexIdx(vi2)<<std::endl;
    }

    out<<"CELL_TYPES "<<m_num_edges<<std::endl;
    for(size_t i=0; i<m_num_edges; ++i){
        out<<"3"<<std::endl;
    }
}

void AlphaShape::writeFacets(){
    std::ofstream out("triangulation_facets.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Triangulation Facets"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    out<<"POINTS "<<m_num_vertices<<" double"<<std::endl;
    for(auto vi = m_triangulation.vertices_begin(); vi != m_triangulation.vertices_end(); ++vi){
        if(!validVertex(vi)) continue;
        out<<vi->point().x()<<" "<<vi->point().y()<<" "<<vi->point().z()<<std::endl;
    }

    out<<"CELLS "<<m_num_facets<<" "<<4*m_num_facets<<std::endl;
    for(auto fi = m_triangulation.facets_begin(); fi != m_triangulation.facets_end(); ++fi){
        if(!validFacet(fi)) continue;
        out<<"3 ";
        for(int i=0; i<4; ++i){
            if(i==fi->second) continue;
            out<<vertexIdx(fi->first->vertex(i))<<" ";
        }
        out<<std::endl;
    }

    out<<"CELL_TYPES "<<m_num_facets<<std::endl;
    for(size_t i=0; i<m_num_facets; ++i){
        out<<"5"<<std::endl;
    }
}

void AlphaShape::writeCells(){
    std::ofstream out("triangulation_cells.vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Triangulation Cells"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

    out<<"POINTS "<<m_num_vertices<<" double"<<std::endl;
    for(auto vi = m_triangulation.vertices_begin(); vi != m_triangulation.vertices_end(); ++vi){
        if(!validVertex(vi)) continue;
        out<<vi->point().x()<<" "<<vi->point().y()<<" "<<vi->point().z()<<std::endl;
    }

    out<<"CELLS "<<m_num_cells<<" "<<5*m_num_cells<<std::endl;
    for(auto ci = m_triangulation.cells_begin(); ci != m_triangulation.cells_end(); ++ci){
        if(!validCell(ci)) continue;
        out<<"4 ";
        for(int i=0; i<4; ++i){
            out<<vertexIdx(ci->vertex(i))<<" ";
        }
        out<<std::endl;
    }

    out<<"CELL_TYPES "<<m_num_cells<<std::endl;
    for(size_t i=0; i<m_num_cells; ++i){
        out<<"10"<<std::endl;
    }
}

void AlphaShape::writeSortedAlpha(){
    std::ofstream out("sorted_alpha.txt");

    std::vector<double> alpha;
    alpha.reserve(m_num_vertices+m_num_edges+m_num_facets+m_num_cells);
    
    for(int i=0; i<m_num_vertices; ++i){
        alpha.push_back(0);
    }

    for(int i=0; i<m_num_edges; ++i){
        auto ei = m_idx_to_edge[i];
        alpha.push_back(m_edge_alpha[ei]);
    }

    for(int i=0; i<m_num_facets; ++i){
        auto fi = m_idx_to_facet[i];
        alpha.push_back(m_facet_alpha[fi]);
    }

    for(int i=0; i<m_num_cells; ++i){
        auto ci = m_idx_to_cell[i];
        alpha.push_back(m_cell_alpha[ci]);
    }

    std::sort(alpha.begin(), alpha.end());

    for(int i=0; i<alpha.size(); ++i){
        out<<alpha[i]<<std::endl;
    }
    out.close();
}
