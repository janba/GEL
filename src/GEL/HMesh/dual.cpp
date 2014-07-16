/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "dual.h"

#include <vector>
#include "../CGLA/Vec3d.h"

#include "Manifold.h"
#include "AttributeVector.h"

namespace HMesh
{
    using namespace std;
    using namespace CGLA;
    
    void dual(Manifold& m)
    {
        // Create new vertices. Each face becomes a vertex whose position
        // is the centre of the face
        int i = 0;
        FaceAttributeVector<int> ftouched;
        vector<Vec3d> vertices;
        vertices.resize(m.no_faces());
        for(auto f : m.faces())
            vertices[ftouched[f] = i++] = centre(m, f);
        
        // Create new faces. Each vertex is a new face with N=valency of vertex
        // edges.
        vector<int> faces;
        vector<int> indices;
        for(auto v : m.vertices())
            if(valency(m, v) > 2 && !(boundary(m, v)))
            {
				int N = circulate_vertex_ccw(m, v, (std::function<void(FaceID)>)[&](FaceID fid) {
                    indices.push_back(ftouched[fid]);
                });
                // Insert face valency in the face vector.
                faces.push_back(N);
            }
        
        // Clear the manifold before new geometry is inserted.
        m.clear();
        
        // And build
        m.build(    vertices.size(),
                reinterpret_cast<double*>(&vertices[0]),
                faces.size(),
                &faces[0],
                &indices[0]);
    }
}