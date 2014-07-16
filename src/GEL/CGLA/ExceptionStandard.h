/** @file ExceptionStandard.h
 * @brief Exceptions are not much used in CGLA, but these classes define what we throw.
 */

/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * @brief Exceptions are not much used in CGLA, but these classes define what we throw.
 * ----------------------------------------------------------------------- */

#ifndef __CGLA_EXCEPTIONSTANDARD_H__
#define __CGLA_EXCEPTIONSTANDARD_H__

#include <string>
#include <iostream>

namespace CGLA
{

	class CGLAMotherException
	{
		std::string str;
	public:
		CGLAMotherException(const std::string s)
			{
				str = s;
			}
  
		void print(std::ostream& os) const 
		{
			os << str << std::endl; 
		}
	};

#define CGLA_DERIVEEXCEPTION(nameoe)															\
	class nameoe: public CGLAMotherException									\
	{																													\
	public:																										\
		nameoe(const std::string& s): CGLAMotherException(s) {}	\
	};																												\

}

#endif
