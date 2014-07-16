/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "ConnectivityKernel.h"

namespace HMesh
{
    using namespace std;

    void ConnectivityKernel::cleanup(IDRemap& map)
    {
        VertexID::IndexType vid = 0;
        for(VertexID v = vertices.index_begin(); v != vertices.index_end(); v = vertices.index_next(v), ++vid)
            map.vmap[v] = VertexID(vid);

        FaceID::IndexType fid = 0;
        for(FaceID f = faces.index_begin(); f != faces.index_end(); f = faces.index_next(f), ++fid)
            map.fmap[f] = FaceID(fid);

        HalfEdgeID::IndexType hid = 0;
        for(HalfEdgeID h = halfedges.index_begin(); h != halfedges.index_end(); h = halfedges.index_next(h), ++hid){
            map.hmap[h] = HalfEdgeID(hid);
        }

        //2. update the connectivity kernel connectivity with the new locations
        for(VertexID v = vertices.index_begin(); v != vertices.index_end(); v = vertices.index_next(v))
            set_out(v, map.hmap[out(v)]);

        for(FaceID f = faces.index_begin(); f != faces.index_end(); f = faces.index_next(f))
            set_last(f, map.hmap[last(f)]);

        for(HalfEdgeID h = halfedges.index_begin(); h != halfedges.index_end(); h = halfedges.index_next(h)){ 
            // do not update holes
            if(face(h) != InvalidFaceID)
                set_face(h, map.fmap[face(h)]);
            set_next(h, map.hmap[next(h)]);
            set_prev(h, map.hmap[prev(h)]);
            set_opp(h, map.hmap[opp(h)]);
            set_vert(h, map.vmap[vert(h)]);
        }

        vertices.cleanup();
        faces.cleanup();
        halfedges.cleanup();
    }
}