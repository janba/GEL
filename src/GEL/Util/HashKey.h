/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file HashKey.h
 * @brief Kehy class for hash tables.
 */
#ifndef __UTIL_HASHKEY_H
#define __UTIL_HASHKEY_H

#include <stdlib.h>
#include <limits.h>
#include "../CGLA/Vec3uc.h"
#include "../CGLA/Vec3usi.h"

namespace Util
{

	
	extern int randoms1[UCHAR_MAX];
	extern int randoms2[UCHAR_MAX];
	extern int randoms3[UCHAR_MAX];

	bool init_randoms();
	void do_init_randoms();
	
	struct HashKey3uc
	{
		CGLA::Vec3uc key;

		HashKey3uc() {do_init_randoms();}
		HashKey3uc(CGLA::Vec3uc _key): key(_key) {do_init_randoms();}
		HashKey3uc(CGLA::Vec3i _key): key(_key) {do_init_randoms();}
		
		int hash(int use_size) const
		{
			return int((randoms1[key[0]] >> (key[1]&0x0f)) +
								 (randoms2[key[1]] >> (key[2]&0x0f)) +
								 (randoms3[key[2]] >> (key[0]&0x0f))) & (use_size-1);
		}
		
		bool operator==(const HashKey3uc& k2) const {return key==k2.key;}
		bool operator!=(const HashKey3uc& k2) const {return !(key==k2.key);}
	};

	struct HashKey3usi
	{
		CGLA::Vec3usi key;

		HashKey3usi() {do_init_randoms();}
		HashKey3usi(CGLA::Vec3usi _key): key(_key) {do_init_randoms();}
		HashKey3usi(CGLA::Vec3i _key): key(_key) {do_init_randoms();}

		int hash(int use_size) const
		{
			return int(
								 ((randoms1[key[0]&0xff00>>8] * randoms2[key[1]&0xff] >> (key[2]&0x0f))
									+ (randoms2[key[1]&0xff00>>8] * randoms1[key[2]&0xff] >> (key[0]&0x0f))
									+ (randoms3[key[2]&0xff00>>8] * randoms3[key[0]&0xff] >> (key[1]&0x0f)))
								 & (use_size-1));
		}
	
		bool operator==(const HashKey3usi& k2) const {return key==k2.key;}
		bool operator!=(const HashKey3usi& k2) const {return !(key==k2.key);}
	};

	struct HashKey1c
	{
		unsigned char key;

		HashKey1c() {do_init_randoms();}
		HashKey1c(unsigned char _key): key(_key) {do_init_randoms();}

		int hash(int use_size) const
		{
			return int(randoms1[key] & (use_size-1));
		}
	
		bool operator==(const HashKey1c& k2) const {return key==k2.key;}
		bool operator!=(const HashKey1c& k2) const {return !(key==k2.key);}
	};
}


#endif
