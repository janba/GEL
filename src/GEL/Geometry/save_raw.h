/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file save_raw.h
 * @brief Save a voxel grid.
 */

#ifndef __GEOMETRY_SAVE_RAW_H
#define __GEOMETRY_SAVE_RAW_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include "GridAlgorithm.h"

namespace Geometry 
{
	template<class G>
		class VolSaver
		{
		public:
			typedef typename G::DataType DataType;
			
		private:
			std::ofstream of;
			const float min_val, max_val;
			const float diff;
			float old;
		public:
			
 			VolSaver(const std::string& name, float _min_val, float _max_val): 
				of(name.c_str(), std::ios::binary),
				min_val(_min_val), 
				max_val(_max_val),
				diff(max_val - min_val) {}
			
			void operator()(const CGLA::Vec3i& pi, const float& vox_val)
				{
					float scaled = (vox_val-min_val) / diff;
					float clamped = 255.0f *(std::min(1.0f,std::max(0.0f,scaled)));
					unsigned char x = static_cast<unsigned char>(clamped);
					of.write((char*) &x, 1);
				}	
		};

	template<class G>
		class VolSaverAscii
		{
		public:
			typedef typename G::DataType DataType;
			
		private:
			std::ofstream of;
			const float min_val, max_val;
			float old;
		public:
			
 			VolSaverAscii(const std::string& name, float _min_val, float _max_val): 
				of(name.c_str()),
				min_val(_min_val), 
				max_val(_max_val) {} 
			
			void operator()(const CGLA::Vec3i& pi, const float& vox_val)
				{
					if(vox_val > min_val && vox_val < max_val)
						of << pi[0] << " " << pi[1] << " " << pi[2]
							 << " " << vox_val << std::endl;
				}	
		};
	

	template<class G>
		class VolSaverFloat
		{
		public:
			typedef typename G::DataType DataType;
			
		private:
			std::ofstream of;
		public:
			
			VolSaverFloat(const std::string& name): 
				of(name.c_str(), std::ios::binary) {}
			
			void operator()(const CGLA::Vec3i& pi, const float& vox_val)
				{
					of.write((char*) &vox_val, sizeof(float));
				}	
		};
	
	template<class G>
		void save_raw_float(const std::string& name, G& grid)
		{
			VolSaverFloat<G> vs(name);
			for_each_voxel_ordered_const(grid, vs);	
		}
	
	template<class G>
		void save_raw_byte(const std::string& name, G& grid,
											 const typename G::DataType& min_val,
											 const typename G::DataType& max_val)
		{
			VolSaver<G> vs(name, min_val, max_val);
			for_each_voxel_ordered_const(grid, vs);	
		}

	template<class G>
		void save_raw_ascii(const std::string& name, G& grid,
												const typename G::DataType& min_val,
												const typename G::DataType& max_val)
		{
			VolSaverAscii<G> vs(name, min_val, max_val);
			for_each_voxel_ordered_const(grid, vs);	
		}

	
}

#endif
