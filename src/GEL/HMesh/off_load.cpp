/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "off_load.h"

#include <fstream>

#include "Manifold.h"

using namespace std;
using namespace CGLA;
using namespace HMesh;
using namespace Geometry;

namespace HMesh
{
    bool off_load(const std::string& filename, HMesh::Manifold& m)
    {
        ifstream ifs(filename.c_str());
        string str;
        if(ifs.good()){
            ifs >> str;
        }
        else return false;
        if(str != "OFF") {
            return false;
        }

        size_t NF, NV, NE;

        ifs >> NV >> NF >> NE;

        vector<Vec3f> vertices(NV);
        for(size_t i = 0; i < NV; ++i){
            ifs >> vertices[i];
        }

        vector<int> faces(NF);
        vector<int> indices;
        for(size_t i = 0; i < NF; ++i){
            int verts_in_face;
            ifs >> verts_in_face;
            faces[i] = verts_in_face;
            for(int j = 0; j < verts_in_face; ++j){
                int idx;
                ifs >> idx;
                indices.push_back(idx);
            }
        }
        m.build(NV, reinterpret_cast<float *>(&vertices[0]), NF, &faces[0], &indices[0]);
        return true;
    }
}

