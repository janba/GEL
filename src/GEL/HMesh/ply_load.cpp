/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <vector>
#include <GEL/HMesh/ply_load.h>

#include <GEL/Geometry/ply_load.h>
#include <GEL/Geometry/TriMesh.h>

#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/load.h>

namespace HMesh
{
    using std::string;
    using Geometry::ply_load;
    using Geometry::TriMesh;
    using namespace CGLA;
    bool ply_load(const string& filename, Manifold& m)
    {
        TriMesh mesh;
        if(Geometry::ply_load(filename, mesh))
        {
            std::vector<Vec3d> vertices_as_doubles(mesh.geometry.no_vertices());
            for(auto idx=0;idx<mesh.geometry.no_vertices(); ++idx)
                vertices_as_doubles[idx] = Vec3d(mesh.geometry.vertex(idx));
            
            std::vector<int> faces(mesh.geometry.no_faces(),3);

            build(m, static_cast<size_t>(mesh.geometry.no_vertices()),
                    reinterpret_cast<const double*>(&vertices_as_doubles[0]),
                    static_cast<size_t>(mesh.geometry.no_faces()),
                    static_cast<const int*>(&faces[0]),
                    reinterpret_cast<const int*>(&mesh.geometry.face(0)));
            return true;
        }
        return false;
    }
}
