/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/HMesh/stl_save.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <GEL/CGLA/Vec3f.h>

#include <GEL/HMesh/AttributeVector.h>
#include <GEL/HMesh/Manifold.h>

namespace HMesh {
    using namespace std;
    using namespace CGLA;

    bool stl_save_ascii(const string& filename, Manifold& m) {

        ofstream os(filename.data());
        if (!os) {
            cerr << "stl_save: could not open file " << filename << endl;
            return false;
        }

        // Build the header
        os << "solid " << filename << endl;

        // Write the faces
        for (FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f) {
            os << "facet normal ";
            Vec3d n = normal(m, *f);
            os << n[0] << " " << n[1] << " " << n[2] << endl;
            os << "\touter loop" << endl;
            for (Walker w = m.walker(*f); !w.full_circle();
                 w        = w.circulate_face_ccw()) {
                Vec3d p = m.pos(w.vertex());
                os << "\t\tvertex " << p[0] << " " << p[1] << " " << p[2]
                   << endl;
            }
            os << "\tendloop" << endl;
            os << "endfacet" << endl;
        }

        // End the file
        os << "endsolid " << filename << endl;

        return true;
    }

    bool stl_save_binary(const string& filename, Manifold& m) {

        ofstream os(filename.data(), ios::binary);
        if (!os) {
            cerr << "stl_save: could not open file " << filename << endl;
            return false;
        }

        // Write the header
        char header[80];
        for (int i = 0; i < 80; ++i) header[i] = 0;
        os.write(header, 80);

        // Write the number of faces (4-byte little-endian unsigned integer)
        unsigned int n_faces = m.no_faces();
        os.write(reinterpret_cast<char*>(&n_faces), 4);

        // Write the faces
        for (FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f) {
            Vec3d n = normal(m, *f);
            normalize(n);

            Vec3f float_n(n[0], n[1], n[2]);
            os.write(reinterpret_cast<char*>(&float_n[0]), 12);
            for (Walker w = m.walker(*f); !w.full_circle();
                 w        = w.circulate_face_ccw()) {
                Vec3d p = m.pos(w.vertex());
                Vec3f float_p(p[0], p[1], p[2]);
                os.write(reinterpret_cast<char*>(&float_p[0]), 12);
            }
            unsigned short attr = 0;
            os.write(reinterpret_cast<char*>(&attr), 2);
        }
        return false;
    }

    bool stl_save(const string& filename, Manifold& m, bool is_binary) {
        if (is_binary)
            return stl_save_binary(filename, m);
        else
            return stl_save_ascii(filename, m);
    }
} // namespace HMesh
