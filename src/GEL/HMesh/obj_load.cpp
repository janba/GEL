/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <sstream>
#include <fstream>
#include <GEL/HMesh/load.h>
#include <GEL/HMesh/obj_load.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/cleanup.h>

using namespace std;
using namespace CGLA;

namespace HMesh
{
    using std::string;
    
    void right_trim(std::string& str) {
        if (!str.empty())
            while(isspace(str.back()))
                str.pop_back();
    }
    
    istream& get_multi_line(istream& is, string& buf) {
        getline(is,buf);
        right_trim(buf);
        if (!buf.empty())
            while(buf.back() == '\\') {
                buf.pop_back();
                string continuation;
                getline(is, continuation);
                right_trim(continuation);
                buf += continuation;
            }
        return is;
    }
    
    bool obj_load(const std::string& filename, Manifold& m, VertexAttributeVector<int>& orig_vertex_indices)
    {
        ifstream obj_file(filename.data());
        
        if(obj_file)
        {
            string buf;
            vector<Vec3d> vertices;
            vector<int> faces;
            vector<int> indices;
            while(get_multi_line(obj_file,buf))
            {
                istringstream iss(move(buf));
                string code;
                if(iss>>code) {
                    if(code == "v") {
                        Vec3d v;
                        iss >> v;
                        vertices.push_back(v);
                    }
                    else if (code == "f") {
                        string str;
                        int ctr = 0;
                        while(iss>>str) {
                            int v,t,n;
                            if(sscanf(str.c_str(), "%d/%d/%d", &v, &t, &n)==3 ||
                               sscanf(str.c_str(), "%d//%d", &v, &n)==2 ||
                               sscanf(str.c_str(), "%d/%d", &v, &t)==2 ||
                               sscanf(str.c_str(), "%d", &v)==1) {
                                indices.push_back(v>0?v-1:vertices.size()+v);
                                ++ctr;
                            }
                        }
                        faces.push_back(ctr);
                    }
                }
            }
//            cout << "Loaded " << filename << " : " << vertices.size() << " vertices and "
//            << faces.size() << " faces"<< endl;
            m.clear();
            
            orig_vertex_indices = build(m, vertices.size(),
                  reinterpret_cast<double*>(&vertices[0]),
                  faces.size(),
                  &faces[0],
                  &indices[0]);
            
            return true;
        }
        return false;
        
    }
    
    bool obj_load(const string& filename, Manifold& m) {
        VertexAttributeVector<int> orig_vertex_indices;
        return obj_load(filename, m, orig_vertex_indices);
    }

}
