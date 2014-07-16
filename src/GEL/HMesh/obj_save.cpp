/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "obj_save.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../CGLA/Vec3f.h"

#include "Manifold.h"
#include "AttributeVector.h"

namespace HMesh
{
    using namespace std;
    using namespace CGLA;

    bool obj_save(const string& filename, Manifold& m)
    {
        ofstream os(filename.data());
        if(os.bad())
            return false;

        VertexAttributeVector<int> vmap;
        int k = 0;
        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
            Vec3d p = m.pos(*v);
            os << "v "<< p[0] << " " << p[1] << " " << p[2] << "\n";
            vmap[*v] = k++;
        }

        for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){        
            vector<int> verts;
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                int idx = vmap[w.vertex()];			
                assert(static_cast<size_t>(idx) < m.no_vertices());
                // move subscript range from 0..size-1 to 1..size according to OBJ standards
                verts.push_back(idx + 1);
            }
			os << "f ";
            for(size_t i = 0; i < verts.size() ; ++i){
                os << verts[i];
                if (i+1==verts.size())
                    break;
                os << " ";
            }
			os<<endl;
        }

        return true;
    }

}
