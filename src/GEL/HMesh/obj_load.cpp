/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <fstream>
#include "obj_load.h"
#include "Manifold.h"
#include "cleanup.h"

using namespace std;
using namespace CGLA;

namespace HMesh
{
    using std::string;
    
    bool obj_load(const string& filename, Manifold& m)
    {
        ifstream obj_file(filename.data());
        
        if(obj_file)
        {
            string buf;
            vector<Vec3d> vertices;
            vector<int> faces;
            vector<int> indices;
            
            while(obj_file >> buf)
            {
                switch(buf[0])
                {
                    case 'v': // v, vn, vt
                    {
                        if(buf == "v")
                        {
                            Vec3d v_geo;
                            obj_file >> v_geo;
                            vertices.push_back(v_geo);
                        }
                        else if(buf == "vn")
                        {
                            Vec3d v_norm;
                            obj_file >> v_norm;
                        }
                        else if(buf == "vt")
                        {
                            Vec3d v_tex;
                            obj_file >> v_tex[0] >> v_tex[1];
                            v_tex[2]=1;
                        }
                    }
                        break;
                    case 'f':
                    {
                        char line[1024];
                        
                        vector<int> v_indices;
                        vector<int> t_indices;
                        vector<int> n_indices;
                        
                        obj_file.getline(line, 1022);
                        char* pch = strtok(line, " \t");
                        int ctr = 0;
                        while(pch != 0)
                        {
                            int v,t,n;
                            if(sscanf(pch, "%d/%d/%d", &v, &t, &n)==3)
                            {
                                indices.push_back(v>0?v-1:vertices.size()+v);
                                ++ctr;
                            }
                            else if(sscanf(pch, "%d//%d", &v, &n)==2)
                            {
                                indices.push_back(v>0?v-1:vertices.size()+v);
                                ++ctr;
                            }
                            else if(sscanf(pch, "%d/%d", &v, &t)==2)
                            {
                                indices.push_back(v>0?v-1:vertices.size()+v);
                                ++ctr;
                            }
                            else if(sscanf(pch, "%d", &v)==1)
                            {
                                indices.push_back(v>0?v-1:vertices.size()+v);
                                ++ctr;
                            }
                            pch = strtok(0, " \t");
                        }
                        faces.push_back(ctr);
                    }
                        break;
                    case 'm':
                        obj_file >> buf;
                        break;
                    case 'u':
                        obj_file >> buf;
                        break;
                    case '#':
                    default:
                        obj_file.ignore(1024, '\n');
                        break;
                }
            }
            cout << "Loaded " << vertices.size() << " vertices and " << faces.size() << " faces"<< endl;
#if 0 // old school loading
            m.clear();
            m.build(vertices.size(),
                    reinterpret_cast<double*>(&vertices[0]),
                    faces.size(),
                    &faces[0],
                    &indices[0]);
#else // robust loading
            m.clear();
            int k=0;
            VertexAttributeVector<int> cluster_id;
            for(int i=0;i<faces.size();++i) {
                vector<Vec3d> pts(faces[i]);
                for(int j=0;j<faces[i]; ++j)
                    pts[j] = vertices[indices[j+k]];
                FaceID f = m.add_face(pts);
                int j=0;
                circulate_face_ccw(m, f, [&](VertexID v){
                    cluster_id[v] = indices[j+k];
                    ++j;
                });
                k += faces[i];
            }
            stitch_mesh(m, cluster_id);
#endif
            return true;
        }
        return false;
        
    }
}