/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <vector>
#include <GEL/HMesh/ply_load.h>

#include <GEL/Geometry/rply.h>
#include <GEL/Geometry/ply_load.h>
#include <GEL/Geometry/TriMesh.h>

#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/load.h>
#include <GEL/HMesh/cleanup.h>

using namespace CGLA;
using namespace std;

namespace HMesh
{
    using std::string;
    using Geometry::ply_load;
    vector<Vec3d> vertices;
    Manifold* mani;

    int vertex_cb(p_ply_argument argument) {
        static int idx=0;
        static Vec3d p;
        int eol;
        ply_get_argument_user_data(argument, NULL, &eol);
        if(idx<3)
            p[idx] = ply_get_argument_value(argument);
        ++idx;
        if (eol)
        {
            vertices.push_back(p);
            idx=0;
        }
        return 1;
    }

    int face_cb(p_ply_argument argument) {
        static vector<int> f;
        int length, value_index;
        ply_get_argument_property(argument, NULL, &length, &value_index);
        if(value_index==-1)
            f.resize(length);
        else
        {
            f[value_index] = ply_get_argument_value(argument);
            if(value_index==(length-1))
            {
                vector<Vec3d> pts(f.size());
                int i=0;
                for(int idx: f)
                    pts[i++] = vertices[idx];
                mani->add_face(pts);
            }
        }
        return 1;
    }

    bool ply_load(const string& fn, Manifold& m) {
        mani = &m;
        p_ply ply = ply_open(fn.c_str(), NULL);
        if (!ply) return false;
        if (!ply_read_header(ply)) return false;
        ply_set_read_cb(ply, "vertex", "x", vertex_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "y", vertex_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "z", vertex_cb, NULL, 1);
        ply_set_read_cb(ply, "face", "vertex_indices", face_cb, NULL, 0);
        ply_read(ply);
        ply_close(ply);
        stitch_mesh(m, 1e-30);
        m.cleanup();
        return true;
    }
}
