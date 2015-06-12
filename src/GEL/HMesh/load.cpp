/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */



#include "load.h"

#include "Manifold.h"

#include "ply_load.h"
#include "x3d_load.h"
#include "obj_load.h"
#include "off_load.h"
#include "cleanup.h"

using namespace std;
using namespace CGLA;

namespace HMesh
{
    using std::string;
    bool load(const string& file_name, Manifold& mani)
    {
        if(file_name.length()<5){
            return false;
        }
        if(file_name.substr(file_name.length()-4,file_name.length())==".obj"){
            return obj_load(file_name, mani);
        }
        else if(file_name.substr(file_name.length()-4,file_name.length())==".x3d"){
            return x3d_load(file_name, mani);
        }
        else if(file_name.substr(file_name.length()-4,file_name.length())==".ply"){
            return ply_load(file_name, mani);
        }
        else if(file_name.substr(file_name.length()-4,file_name.length())==".off"){
            return off_load(file_name, mani);
        }
        return false;
    }
    
    
     void safe_build(Manifold& m, size_t no_vertices,
                              const double* vertvec,
                              size_t no_faces,
                              const int* facevec,
                              const int* indices)
    {
        int k=0;
        VertexAttributeVector<int> cluster_id;
        for(int i=0;i<no_faces;++i) {
            vector<Vec3d> pts(facevec[i]);
            for(int j=0;j<facevec[i]; ++j) {
                const double* v = &vertvec[3*indices[j+k]];
                pts[j] = Vec3d(v[0],v[1],v[2]);
            }
            FaceID f = m.add_face(pts);
            int j=0;
            circulate_face_ccw(m, f, (std::function<void(VertexID)>)[&](VertexID v){
                cluster_id[v] = indices[j+k];
                ++j;
            });
            k += facevec[i];
        }
        stitch_mesh(m, cluster_id);
    }

}