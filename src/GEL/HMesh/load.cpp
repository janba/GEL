/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */



#include <GEL/HMesh/load.h>

#include <GEL/HMesh/Manifold.h>

#include <GEL/HMesh/ply_load.h>
#include <GEL/HMesh/x3d_load.h>
#include <GEL/HMesh/obj_load.h>
#include <GEL/HMesh/off_load.h>
#include <GEL/HMesh/stl_load.h>
#include <GEL/HMesh/cleanup.h>

#include <GEL/Util/Serialization.h>

using namespace std;
using namespace CGLA;
using namespace Util;

namespace HMesh
{
    using std::string;
    bool load(const string& file_name, Manifold& mani)
    {
        if(file_name.length()<5){
            return false;
        }
        string ext = file_name.substr(file_name.length()-4,file_name.length());
        if(ext==".obj"){
            return obj_load(file_name, mani);
        }
        else if(ext==".x3d"){
            return x3d_load(file_name, mani);
        }
        else if(ext==".ply"){
            return ply_load(file_name, mani);
        }
        else if(ext==".off"){
            return off_load(file_name, mani);
        }
        else if(ext==".off"){
            return off_load(file_name, mani);
        }
        else if(ext==".bhm") {
            Serialization ser(file_name, std::ios_base::in);
            mani.deserialize(ser);
            return true;
        }
        else if(ext==".stl") {
            return stl_load(file_name, mani);
        }
        return false;
    }
}
