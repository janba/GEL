/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/HMesh/refine_edges.h>

#include <vector>
#include <iterator>

#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/AttributeVector.h>

namespace HMesh
{
    using namespace std;

    float average_edge_length(const Manifold& m)
    {
        float lsum = 0;
        for(auto h : m.halfedges())
            lsum += length(m, h);
        return lsum / m.no_halfedges();
    }
    
    float median_edge_length(const Manifold& m)
    {
        vector<double> lengths;
        for(auto h : m.halfedges()) {
            auto w = m.walker(h);
            if(h == w.hmin())
                lengths.push_back(sqr_length(m.pos(w.vertex())-m.pos(w.opp().vertex())));
        }
        if (lengths.size()>0) {
            nth_element(begin(lengths), begin(lengths)+lengths.size()/2, end(lengths));
            return sqrt(lengths[lengths.size()/2]);
        }
        return 0;
    }


    int refine_edges(Manifold& m, float t)
    {
        vector<HalfEdgeID> hedges;
        hedges.reserve(m.no_halfedges());
        copy(m.halfedges_begin(), m.halfedges_end(), back_inserter(hedges));

        HalfEdgeAttributeVector<int> touched(m.allocated_halfedges(), 0);

        int work = 0;
        for(HalfEdgeID h : hedges){
            Walker w = m.walker(h);

            if(!m.in_use(h) || w.face() == InvalidFaceID || length(m, h) < t || touched[h])
                continue;
            touched[w.opp().halfedge()] = 1;
            m.split_edge(h);
            ++work;
        }
        return work;
    }

}
