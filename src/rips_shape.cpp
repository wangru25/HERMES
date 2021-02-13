#include "rips_shape.hpp"

RipsShape::RipsShape(){

}

RipsShape::~RipsShape(){

}

int RipsShape::numVertices() const{
    return m_num_vertices;
}

int RipsShape::numEdges() const{
    return m_num_edges;
}

int RipsShape::numFacets() const{
    return m_num_facets;
}

int RipsShape::numCells() const{
    return m_num_cells;
}

void RipsShape::readPointCloud(std::string filename){
    std::ifstream file(filename);
    
    std::string line;
    while(std::getline(file, line)){
        std::stringstream ss(line);
        double x, y, z;
        ss>>x>>y>>z;
        std::vector<double> tmp;
        tmp.push_back(x);
        tmp.push_back(y);
        tmp.push_back(z);
        m_point_cloud.push_back(tmp);
    }
    file.close();
}

void RipsShape::convertToDistanceMatrix(){
    //init
    for(int i = 0; i < m_point_cloud.size(); i++){
        std::vector<double> tmp(m_point_cloud.size(), 0.);
        distance.push_back(tmp);
    }
    
    for(int i = 0; i < m_point_cloud.size(); i++){
        for(int j = i+1; j < m_point_cloud.size(); j++){
            double dist = std::sqrt(std::inner_product(m_point_cloud[i].begin(), m_point_cloud[i].end(),
                m_point_cloud[j].begin(), 0., std::plus<double>(), [](double u, double v)
                { return (u - v) * (u - v); }));

            distance[i][j] = dist;
            distance[j][i] = dist;
        }
    }
}

void RipsShape::readDistanceMatrix(std::string filename){
    std::ifstream file(filename);
    
    std::string line;
    int col_index = 0;
    int row_index = 0;
    while(std::getline(file, line)){
        std::stringstream ss(line);
        std::vector<double> temp;
        double x;
        while (ss>>x){
            temp.push_back(x);
            col_index++;
        }
        distance.push_back(temp);
        col_index = 0;
        row_index++;
    }
}

void RipsShape::preprocess(){
    m_num_vertices = distance.size();
    
    int e_index = 0;
    int f_index = 0;
    int c_index = 0;
    for (int i = 0; i < m_num_vertices; i++){
        SimpleVertex v(i);
        m_vertex_idx[v] = i;
        m_idx_to_vertex[i] = v;
        
//        // build adjacency lists
//        int adjacency_list[m_num_vertices];
//
//        std::map<double, int> distance_ij;
//        for(int j = 0; j < m_num_vertices; j++){
//            // leave out vertex farther than threshold?
//            if (i != j) //&& distance[i][j] <= distance_threshold)
//                distance_ij.insert(std::make_pair(distance[i][j], j));
//        }
//
//        auto it = distance_ij.begin();
//        for(int j = 0; j < m_num_vertices; j++){
//            if (i != j)
//                adjacency_list[j] = it++->second;
//        }
//        v.adj_vertices = adjacency_list;
        
        // TODO: loop through adj list of i until reach threshold distance
        for (int j = i+1; j < m_num_vertices; j++){
            
            SimpleEdge e(e_index, i, j, distance[i][j]);
            m_idx_to_edge[e_index] = e;
            m_edge_idx[e] = e_index++;
            m_vertices_to_edge[std::make_pair(e.v1, e.v2)] = e;
            m_simple_edges.push_back(e);
            
            // TODO: loop through adj list of j until reach threshold distance
            for (int k = j+1; k < m_num_vertices; k++){
                double dist_ij, dist_jk, dist_ik;
                dist_ij = distance[i][j];
                dist_jk = distance[j][k];
                dist_ik = distance[i][k];
                double maxdistance = std::max(std::max(dist_ij, dist_jk), dist_ik);
                
                SimpleFacet f(f_index, i, j, k, maxdistance);
                m_idx_to_facet[f_index] = f;
                m_facet_idx[f] = f_index++;
                m_vertices_to_facet[std::make_pair(i, std::make_pair(j, k))] = f;
                m_simple_facets.push_back(f);
                
                double max_d_face = maxdistance;
                // TODO: loop through adj list of k until reach threshold distance
                for (int l = k+1; l < m_num_vertices; l++){
                    double dist_il, dist_jl, dist_kl;
                    dist_il = distance[i][l];
                    dist_jl = distance[j][l];
                    dist_kl = distance[k][l];
                    maxdistance = std::max(std::max(std::max(max_d_face, dist_il), dist_jl), dist_kl);

                    SimpleCell c(c_index, i, j, k, l, maxdistance);
                    m_idx_to_cell[c_index] = c;
                    m_cell_idx[c] = c_index++;
                    m_simple_cells.push_back(c);
                }
            }
        }
    }
    m_num_edges = e_index;
    m_num_facets = f_index;
    m_num_cells = c_index;

    std::sort(m_simple_edges.begin(), m_simple_edges.end(), [](SimpleEdge a, SimpleEdge b) { return a.length < b.length; });
    const double t = distance_threshold;
    auto ite = std::find_if(std::begin(m_simple_edges), std::end(m_simple_edges), [t](SimpleEdge e) { return t < e.length; });
    // remove all edges starting at it
    m_simple_edges.erase(ite, m_simple_edges.end());
    
    for(int i = 0; i < m_simple_edges.size(); i++){
        SimpleEdge e = m_simple_edges[i];
        
        positiveEdgeMap[i] = e.v1;
        negativeEdgeMap[i] = e.v2;
    }
    
    std::sort(m_simple_facets.begin(), m_simple_facets.end(), [](SimpleFacet a, SimpleFacet b) {
        return a.maxdistance < b.maxdistance;});
    auto itf = std::find_if(std::begin(m_simple_facets), std::end(m_simple_facets), [t](SimpleFacet f) { return t < f.maxdistance; });
    // remove all facets starting at it
    m_simple_facets.erase(itf, m_simple_facets.end());
    
    for(int i = 0; i < m_simple_facets.size(); i++){
        SimpleFacet f = m_simple_facets[i];

        //exclude v3, +1 (1-based indexing)
        SimpleEdge e = m_vertices_to_edge[std::make_pair(f.v1, f.v2)];
        auto it = std::find(std::begin(m_simple_edges), std::end(m_simple_edges), e);
        int edge_index3 = it - m_simple_edges.begin();
        
        //exclude v1, +1
        e = m_vertices_to_edge[std::make_pair(f.v2, f.v3)];
        it = std::find(std::begin(m_simple_edges), std::end(m_simple_edges), e);
        int edge_index1 = it - m_simple_edges.begin();
        
        //exclude v2, -1
        e = m_vertices_to_edge[std::make_pair(f.v1, f.v3)];
        it = std::find(std::begin(m_simple_edges), std::end(m_simple_edges), e);
        int edge_index2 = it - m_simple_edges.begin();
        
        positiveFacetMap[i] = std::make_pair(edge_index3, edge_index1);
        negativeFacetMap[i] = edge_index2;
    }
    
    std::sort(m_simple_cells.begin(), m_simple_cells.end(), [](SimpleCell a, SimpleCell b) {
        return a.maxdistance < b.maxdistance;});
    auto itc = std::find_if(std::begin(m_simple_cells), std::end(m_simple_cells), [t](SimpleCell c) { return t < c.maxdistance; });
    // remove all facets starting at it
    m_simple_cells.erase(itc, m_simple_cells.end());
    
    for(int i = 0; i < m_simple_cells.size(); i++){
        SimpleCell c = m_simple_cells[i];

        //exclude v4, -1 (1-based indexing)
        SimpleFacet f = m_vertices_to_facet[std::make_pair(c.v1, std::make_pair(c.v2, c.v3))];
        auto it = std::find(std::begin(m_simple_facets), std::end(m_simple_facets), f);
        int facet_index4 = it - m_simple_facets.begin();
        
        //exclude v3, +1
        f = m_vertices_to_facet[std::make_pair(c.v1, std::make_pair(c.v2, c.v4))];
        it = std::find(std::begin(m_simple_facets), std::end(m_simple_facets), f);
        int facet_index3 = it - m_simple_facets.begin();
        
        //exclude v2, -1
        f = m_vertices_to_facet[std::make_pair(c.v1, std::make_pair(c.v3, c.v4))];
        it = std::find(std::begin(m_simple_facets), std::end(m_simple_facets), f);
        int facet_index2 = it - m_simple_facets.begin();
        
        //exclude v1, +1
        f = m_vertices_to_facet[std::make_pair(c.v2, std::make_pair(c.v3, c.v4))];
        it = std::find(std::begin(m_simple_facets), std::end(m_simple_facets), f);
        int facet_index1 = it - m_simple_facets.begin();
        
        //add to maps to read when building coboundary matrix
        positiveCellMap[i] = std::make_pair(facet_index3, facet_index1);
        negativeCellMap[i] = std::make_pair(facet_index4, facet_index2);
    }
    
}
