/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "x3d_save.h"

#include <fstream>
#include <vector>
#include "../CGLA/Vec3f.h"

#include "Manifold.h"
#include "AttributeVector.h"


namespace HMesh
{
    using namespace std;
    using namespace CGLA;


    namespace
    {
        const string X3D_BEGIN = 
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
            "<!DOCTYPE X3D PUBLIC \"http://www.web3D.org/TaskGroups/x3d/translation/x3d-compact.dtd\"\n"
            "\"/www.web3d.org/TaskGroups/x3d/translation/x3d-compact.dtd\">\n"
            "<X3D>\n"
            "  <head>\n"
            "  </head>\n"
            "  <Scene>\n"
            "		<Viewpoint description=\"Pyramid\" orientation=\"0 1 0 .2\" position=\"4 0 25\"/>\n"
            "    <NavigationInfo type=\"EXAMINE ANY\"/>\n";

        const string X3D_END = 
            "  </Scene>\n"
            "</X3D>\n";

    }

    bool x3d_save(const std::string& filename, Manifold& m)
    {
        ofstream os(filename.data());
        if(os.bad()){ 
            return false;
        }
        os << X3D_BEGIN << "\n";

        os << "<Shape>\n"
            << "<Appearance>\n"
            << "<Material diffuseColor=\"0.8 0.8 0.2\" specularColor=\"0 0 0.5\"/>\n"
            << "</Appearance>\n";

        os <<  "<Coordinate point=\"" << "\n";

        VertexAttributeVector<int> vmap(m.allocated_vertices());

        int k = 0;
        for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
            Vec3d p = m.pos(*v);
            os << p[0] << " " << p[1] << " " << p[2] << "\n";
            vmap[*v] = k++;
        }
        os << "\"/>" << "\n";

        os << "<IndexedFaceSet coordIndex=\"" << endl;
        for(FaceIDIterator f =  m.faces_begin(); f != m.faces_end(); ++f){
            vector<int> verts;
            for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
                int idx = vmap[w.vertex()];
                assert(static_cast<size_t>(idx) < m.no_vertices());
                verts.push_back(idx);
            }
            //assert(verts.size()==3);
            for(size_t i = 0; i < verts.size(); ++i){
                os << verts[i] << " ";
            }
            os << " -1\n";
        }
        os << "\">" << "\n";
        os << "</IndexedFaceSet>\n"
            << "</Shape>\n";
        os << X3D_END << "\n";

        return true;
    }
}
