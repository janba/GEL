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
#include <GEL/HMesh/cleanup.h>

using namespace std;
using namespace CGLA;

namespace HMesh
{
    using std::string;
    bool load(const string& file_name, Manifold& mani)
    {
        if(file_name.length()<5){
            return false;
        }
        if(file_name.substr(file_name.length()-4,file_name.length())==".obj"){
            return obj_load(file_name, mani);
        }
        else if(file_name.substr(file_name.length()-4,file_name.length())==".x3d"){
            return x3d_load(file_name, mani);
        }
        else if(file_name.substr(file_name.length()-4,file_name.length())==".ply"){
            return ply_load(file_name, mani);
        }
        else if(file_name.substr(file_name.length()-4,file_name.length())==".off"){
            return off_load(file_name, mani);
        }
        return false;
    }
}
