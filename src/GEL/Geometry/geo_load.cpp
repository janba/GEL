/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/Geometry/load.h>
#include <GEL/Geometry/obj_load.h>
#include <GEL/Geometry/ply_load.h>

using namespace std;
using namespace CGLA;
//using namespace HMesh;

namespace Geometry
{
	
	void load(const string& fn, TriMesh &mesh)
	{
		if(fn.substr(fn.length()-4,fn.length())==".obj")
		{
			obj_load(fn, mesh);
		}
		else if(fn.substr(fn.length()-4,fn.length())==".ply")
		{
			ply_load(fn, mesh);
		}	
		else
		{
			cout << "Either the format was unrecognized or the file did not have the appropriate extension" << endl;
			exit(0);
		}
	}
	
}

