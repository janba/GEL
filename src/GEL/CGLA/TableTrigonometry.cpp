/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "TableTrigonometry.h"

namespace CGLA {

	using namespace std;

	namespace TableTrigonometry
	{
		const CosTable& COS_TABLE()
		{
			static CosTable table;
			return table;
		}
	}

}
