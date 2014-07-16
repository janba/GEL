/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <fstream>
#include "obj_load.h"
#include "Manifold.h"

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
            vector<Vec3f> vertices;
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
                            Vec3f v_geo;
                            obj_file >> v_geo;
                            vertices.push_back(v_geo);
                        }
                        else if(buf == "vn")
                        {
                            Vec3f v_norm;
                            obj_file >> v_norm;
                        }
                        else if(buf == "vt")
                        {
                            Vec3f v_tex;
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
            m.clear();
            m.build(vertices.size(),
                    reinterpret_cast<float*>(&vertices[0]),
                    faces.size(),
                    &faces[0],
                    &indices[0]);
            
            return true;
        }
        return false;
        
    }
}