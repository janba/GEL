/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "Parse.h"
#include <cstring>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace CGLA;

namespace
{
	const char seps[]   = " ,\t\n\r";

	inline char* next_etoken(const char* btoken)
	{
		return const_cast<char*>(btoken)+strcspn(btoken, seps);
	}
	
	inline const char* next_btoken(char* etoken)
	{
		return const_cast<const char*>(etoken+strspn(etoken, seps));
	}
}

namespace Util {
	string floatToString(float value) {
		stringstream ret;
		string stringret;
		ret << value;
		getline(ret,stringret);
		return stringret;
	}

	void parse(const char* str,bool& x) {
		x = (strcmp(str,"true")==0) || (strcmp(str,"TRUE")==0);
	}

	void parse(const char* str,string& x) {
		x=string(str);
	}

	void parse(const char* str,int& x) {
		/* Establish string and get the first token: */
		x = strtol(str,0,10);
	} 

	void parse(const char* str,Vec2i& vec)
	{   /* Establish string and get the first token: */
		const char *btoken=next_btoken(const_cast<char*>(str));
		char *etoken=next_etoken(btoken);
		vec[0] = strtol(btoken,0,10);
		btoken= next_btoken(etoken);
		etoken= next_etoken(btoken);
		vec[1] = strtol(btoken,0,10);
	}


	void parse(const char* str,float& x) 
	{
		/* Establish string and get the first token: */
		x = strtod(str,0);
	}

	void parse(const char* str,Vec2f& vec) {
		const char *btoken=next_btoken(const_cast<char*>(str));
		char *etoken= next_etoken(btoken);
		vec[0] = strtod(btoken,0);
		btoken=next_btoken(etoken);
		etoken=next_etoken(btoken);
		vec[1] = strtod(btoken,0);
	}

	void parse(const char* str,Vec3f& vec) {
   /* Establish string and get the first token: */
		const char *btoken=next_btoken(const_cast<char*>(str));
		char *etoken= next_etoken(btoken);
		vec[0] = strtod(btoken,0);
		btoken=next_btoken(etoken);
		etoken=next_etoken(btoken);
		vec[1] = strtod(btoken,0);
		btoken=next_btoken(etoken);
		etoken=next_etoken(btoken);
		vec[2] = strtod(btoken,0);
	}


	void parse(const char* str,Vec4f& vec) {
   /* Establish string and get the first token: */
		const char *btoken=next_btoken(const_cast<char*>(str));
		char *etoken= next_etoken(btoken);
		vec[0] = strtod(btoken,0);
		btoken=next_btoken(etoken);
		etoken=next_etoken(btoken);
		vec[1] = strtod(btoken,0);
		btoken=next_btoken(etoken);
		etoken=next_etoken(btoken);
		vec[2] = strtod(btoken,0);
		btoken=next_btoken(etoken);
		etoken=next_etoken(btoken);
		vec[3] = strtod(btoken,0);
	}

	void parse(const char* str,vector<float>& v) {
		const char *btoken=next_btoken(const_cast<char*>(str));
		char *etoken= next_etoken(btoken);
		while(etoken>btoken)
		{
			v.push_back(strtod(btoken,0));
			btoken=next_btoken(etoken);
			etoken=next_etoken(btoken);
		}
	}

	void parse(const char* str,vector<double>& v) {
		const char *btoken=next_btoken(const_cast<char*>(str));
		char *etoken= next_etoken(btoken);
		while(etoken>btoken)
		{
			v.push_back(strtod(btoken,0));
			btoken=next_btoken(etoken);
			etoken=next_etoken(btoken);
		}
	}

  void parse(const char* str,vector<Vec2f>& v) {
		const char *btoken=next_btoken(const_cast<char*>(str));
		char *etoken= next_etoken(btoken);
		while(etoken>btoken)
		{
			Vec2f vec;
			vec[0] = strtod(btoken,0);
			btoken=next_btoken(etoken);
			etoken=next_etoken(btoken);
			vec[1] = strtod(btoken,0);
			btoken=next_btoken(etoken);
			etoken=next_etoken(btoken);
			v.push_back(vec);
		}
	}

	void parse(const char* str,vector<Vec3f>& v) {
   /* Establish string and get the first token: */
		const char *btoken=next_btoken(const_cast<char*>(str));
		char *etoken= next_etoken(btoken);
		while(etoken>btoken)
		{
			Vec3f vec;
			vec[0] = strtod(btoken,0);
			btoken=next_btoken(etoken);
			etoken=next_etoken(btoken);
			vec[1] = strtod(btoken,0);
			btoken=next_btoken(etoken);
			etoken=next_etoken(btoken);
			vec[2] = strtod(btoken,0);
			btoken=next_btoken(etoken);
			etoken=next_etoken(btoken);
			v.push_back(vec);
		}
	}

	void parse(const char* str,vector<int>& v) 
	{
		const char *btoken=next_btoken(const_cast<char*>(str));
		char *etoken= next_etoken(btoken);
		while(etoken>btoken)
			{
				v.push_back(strtol(btoken,0,10));
				btoken=next_btoken(etoken);
				etoken=next_etoken(btoken);
			}
	}
	
}
