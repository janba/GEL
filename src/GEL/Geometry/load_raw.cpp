/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <fstream>
#include "../Geometry/load_raw.h"

using namespace CGLA;
using namespace std;

namespace Geometry
{

	template<class T>
	bool load_raw(const string& file, RGrid<T>& grid)
	{
			int sz = grid.get_size();
			ifstream f(file.c_str(),ios::binary);
			if(f)
			{
					f.read(reinterpret_cast<char*>(grid.get()),sz*sizeof(T));
					return true;
			}
			cerr << "Could not open volume :" << file << endl;
			return false;
	}
		
		template bool load_raw(const string&, RGrid<unsigned char>& grid);
		template bool load_raw(const string&, RGrid<unsigned short>& grid);
		template bool load_raw(const string&, RGrid<short>& grid);
		template bool load_raw(const string&, RGrid<float>& grid);

}
