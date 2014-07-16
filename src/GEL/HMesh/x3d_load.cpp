/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "x3d_load.h"

#include <iostream>
#include <vector>
#include <string>

#include "../CGLA/Vec3f.h"
#include "../Util/Timer.h"
#include "../Util/Parse.h"
#include "../Util/XmlParser.h"

#include "Manifold.h"

namespace HMesh
{
    using namespace std;
    using namespace CGLA;
    using namespace Util;

    namespace
    {
        vector<int> faces;
        vector<int> indices;
        vector<float> vertices;

        void coord_index_to_face_vec(const vector<int>& coord_index, vector<int>& faces, vector<int>& indices)
        {
            int k = 0;
            for(size_t i = 0; i < coord_index.size(); ++i){
                int idx = coord_index[i];
                if (idx == -1) {
                    faces.push_back(k);
                    k=0;
                }
                else{
                    indices.push_back(idx);
                    ++k;
                }
            }
        }
    }

    void handle_Shape(XmlElement& elem)
    {
        cout << "Found Shape" << endl;
        elem.process_elements();
        cout << "Shape ends" << endl;				
    }

    void handle_IndexedFaceSet(XmlElement& elem)
    {
        cout << "Found IndexedFaceSet" << endl;
        vector<int> coord_index;
        parse(elem.atts["coordIndex"].c_str(), coord_index);
        coord_index_to_face_vec(coord_index, faces, indices);
        elem.process_elements();
        cout << "IndexedFaceSet ends" << endl;
    }

    void handle_Coordinate(XmlElement& elem)
    {
        cout << "Found Coordinate" << endl;
        parse(elem.atts["point"].c_str(), vertices);
        cout << "Coordinate ends" << endl;
    }

    int find_last_of(const string& F, const string& C)
    {
        int pos = F.find_last_of(C);
        if (pos == string::npos) 
            return -1;
        return pos;
    }

    bool x3d_load(const string& filename, Manifold& m) 
    {
        faces.clear();
        indices.clear();
        vertices.clear();

        Timer tim;
        tim.start();

        string baseurl;
        int idx = max(find_last_of(filename, "\\"), 
            find_last_of(filename, "/"));

        if(idx != -1)
            baseurl = string(filename.substr(0, idx+1));

        XmlDoc x3d_doc(filename.c_str());
        
        if(!x3d_doc.is_valid())
            return false;
        
        x3d_doc.add_handler("Shape", handle_Shape);    
        x3d_doc.add_handler("IndexedFaceSet", handle_IndexedFaceSet);
        x3d_doc.add_handler("Coordinate", handle_Coordinate);
        x3d_doc.process_elements();
        x3d_doc.close();
        
        cout << "vertices " << vertices.size() << endl;

        m.build(vertices.size()/3, 
            reinterpret_cast<float*>(&vertices[0]), 
            faces.size(), 
            &faces[0], 
            &indices[0]);
        
        cout << " Loading took " << tim.get_secs() << endl;
        return true;
    }
}
