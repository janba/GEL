/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "ply_load.h"

#include "../Geometry/ply_load.h"
#include "../Geometry/TriMesh.h"

#include "Manifold.h"

namespace HMesh
{
    using std::string;
    using Geometry::ply_load;
    using Geometry::TriMesh;
    bool ply_load(const string& filename, Manifold& m)
    {
        TriMesh mesh;
        if(Geometry::ply_load(filename, mesh))
        {
            m.build(mesh);
            return true;
        }
        return false;
    }
}
