/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <time.h>
#include <iostream>
#include "../CGLA/CGLA.h"
#include "HashKey.h"

using namespace CGLA;

namespace Util
{
	
	int randoms1[UCHAR_MAX];
	int randoms2[UCHAR_MAX];
	int randoms3[UCHAR_MAX];


	
	bool init_randoms()
	{
		gel_srand(1);
		int i;
		for(i=0;i<UCHAR_MAX-1;i++)
			randoms1[i] = gel_rand();
		for(i=0;i<UCHAR_MAX-1;i++)
		randoms2[i] = gel_rand();
		for(i=0;i<UCHAR_MAX-1;i++)
			{
				randoms3[i] = gel_rand();
			}
		return true;
	}

	void do_init_randoms()
	{
		init_randoms();
	}

}
