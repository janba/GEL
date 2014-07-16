/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Parse.h
 * @brief Parsing various entites from a string.
 */

#ifndef __UTIL_PARSE_H__
#define __UTIL_PARSE_H__

#include <string>
#include <sstream>
#include <vector>
#include "../CGLA/Vec4f.h"
#include "../CGLA/Vec3f.h"
#include "../CGLA/Vec2f.h"

namespace Util {
	std::string floatToString(float value);
	void parse(const char* str,bool& x);
	void parse(const char* str,std::string& x);
	void parse(const char* str,int& x);
	void parse(const char* str,CGLA::Vec2i&);
	void parse(const char* str,float& x);
	void parse(const char* str,CGLA::Vec2f&);
	void parse(const char* str,CGLA::Vec3f& vec);
	void parse(const char* str,CGLA::Vec4f&);
	void parse(const char* str,std::vector<float>& v);
	void parse(const char* str,std::vector<double>& v);
	void parse(const char* str,std::vector<CGLA::Vec2f>& v);
	void parse(const char* str,std::vector<CGLA::Vec3f>& v);
	void parse(const char* str,std::vector<int>& v);
	void parseSMPT(const char* str, float& x);
}

#endif // __PARSE_H__
