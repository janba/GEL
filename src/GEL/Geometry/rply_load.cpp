/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "rply.h"
#include "ply_load.h"
#include "../CGLA/Vec4f.h"

namespace Geometry
{
    using namespace std;
    using namespace CGLA;

    namespace 
    {
        TriMesh *mesh;

        int vertex_cb(p_ply_argument argument) {
            static int idx=0;
            static Vec3f p;
            int eol;
            ply_get_argument_user_data(argument, NULL, &eol);
            if(idx<3)
                p[idx] = ply_get_argument_value(argument);
            ++idx;
            if (eol) 
            {
                mesh->geometry.add_vertex(p);
                idx=0;
            }
            return 1;
        }

        int face_cb(p_ply_argument argument) {
            static Vec3i f;
            int length, value_index;
            ply_get_argument_property(argument, NULL, &length, &value_index);
            if(value_index >= 0)
            {
                if(value_index < 2)
                    f[value_index] = ply_get_argument_value(argument);	
                else
                {
                    f[2] = ply_get_argument_value(argument);
                    mesh->mat_idx.push_back(0);
                    mesh->geometry.add_face(f);
                    f[1] = f[2];
                }
            }
            return 1;
        }
    }

    bool ply_load(const std::string& fn, Geometry::TriMesh& _mesh)
    {
        mesh = &_mesh;

        _mesh.materials.resize(1);
        _mesh.materials[0].diffuse[0] = 172.0f/256.0f; 
        _mesh.materials[0].diffuse[1] = 48.0f/256.0f;
        _mesh.materials[0].diffuse[2] = 72.0f/256.0f;
        _mesh.materials[0].specular[0] = 0.6f; 
        _mesh.materials[0].specular[1] = 0.6f;
        _mesh.materials[0].specular[2] = 0.6f;
        _mesh.materials[0].shininess = 128.0f;

        p_ply ply = ply_open(fn.c_str(), NULL);
        if (!ply) return false;
        if (!ply_read_header(ply)) return false;
        ply_set_read_cb(ply, "vertex", "x", vertex_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "y", vertex_cb, NULL, 0);
        ply_set_read_cb(ply, "vertex", "z", vertex_cb, NULL, 1);
        ply_set_read_cb(ply, "face", "vertex_indices", face_cb, NULL, 0);
        ply_read(ply);
        ply_close(ply);
        return true;
    }
}