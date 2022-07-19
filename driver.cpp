//
// Created by usr on 18/07/2022.
//

#include <PyGEL/graph_functions.h>
#include "GEL/Geometry/graph_io.h"
#include "GEL/Geometry/Graph.h"

#include "GEL/HMesh/Manifold.h"
#include "chrono"
#include "GEL/Geometry/graph_skeletonize.h"

using Graph = Geometry::AMGraph3D;


void skeletonize(Graph &g, Graph* skel_ptr, bool samplping){
    std::cout << "Generating skeleton" << std::endl;

    auto t0 = std::chrono::high_resolution_clock::now();

    Geometry::local_separators(g,samplping);

    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Time to generate: "<<(t1-t0).count() * 1e-9 << std::endl;
}

void graph_from_mesh(HMesh::Manifold* m, Geometry::AMGraph3D* g) {
    HMesh::VertexAttributeVector<Geometry::AMGraph::NodeID> v2n;

    for(auto v : m->vertices())
        v2n[v] = g->add_node(m->pos(v));
    for(auto h: m->halfedges()) {
        HMesh::Walker w = m->walker(h);
        if(h<w.opp().halfedge())
            g->connect_nodes(v2n[w.opp().vertex()], v2n[w.vertex()]);
    }
}

int main(int argc, char* argv[]){
    bool success;
    auto g = Graph();

    std::string path = "../package/fertility.off";

    bool graph_is_mesh = false;

    HMesh::Manifold m;
    // If source is mesh, load and convert to graph
    if(path.compare(path.size()-6,6,".graph")){
        m = HMesh::Manifold();
        success = HMesh::load(path, m);
        if (!success) {
            std::cout << "ERROR : Could not load mesh!";
            return 0;
        }
        graph_from_mesh(&m,&g);
    } else {
        g = Geometry::graph_load(path);
    }

    success = !g.empty();
    if (!success) std::cout << "ERROR : Graph is empty!" << std::endl;

    // Make skeleton
    auto skel = Geometry::AMGraph3D(); // Target graph

    skeletonize(g,&skel,true);

}
