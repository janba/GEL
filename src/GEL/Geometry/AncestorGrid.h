/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file AncestorGrid.h
 * @brief Abstract class for voxel grids.
 */

#ifndef __GEOMETRY_ANCESTORGRID_H
#define __GEOMETRY_ANCESTORGRID_H

#include "../CGLA/Vec3i.h"

namespace Geometry 
{

	/** \brief Class template is used as abstract ancestor of 
			voxel grids. 

			Strictly speaking, this class is not 
			abstract, since it does not have any virtual functions.
			However, operator[]() and store() simply call 
			functions in derived classes. To do so, you must pass
			the derived class as a template argument to this class
			when you define the derived class. This is called the
			Barton and Nackman trick. See Todd Veldhuizen, 
			"Techniques for Scientific C++" 1.3.3.
	*/
 	template<typename T, class ChildT>
	class AncestorGrid
	{
	public:
		typedef T DataType;

	private:
		/// xyz dimensions of grid.
		CGLA::Vec3i dims;
  
	public:
  
		/// Construct a grid of specified xyz dimensions.
		AncestorGrid(int _x_dim, int _y_dim, int _z_dim): 
			dims(_x_dim,_y_dim,_z_dim) {}
  
		/// Construct a grid of specified xyz dimensions.
		AncestorGrid(const CGLA::Vec3i& _dims): dims(_dims) {}
		
		/** Check if voxel is within grid bounds.
				This function is passed a Vec3i, p, and returns true
				if p is within the voxel grid. */
 		bool in_domain(const CGLA::Vec3i& p) const
		{
			for(int i=0; i<3; i++)
				if (p[i]<0 || p[i] >= dims[i])
					return false;
			return true;
		}
  
		/** Get dimensions of grid.

				This function returns a Vec3i with the dimensions
				of the grid. */
 		const CGLA::Vec3i& get_dims() const {return dims;}

		/// Get the corner having smallest coordinates.
 		const CGLA::Vec3i get_lo_corner() const {return CGLA::Vec3i(0);}

		/// Get the corner having greatest coordinates. 
 		const CGLA::Vec3i& get_hi_corner() const {return dims;}

		/** Access (read only) a voxel in a grid. 

				This is the operator[] which is passed a Vec3i 
				and returns a const reference to a voxel.
				This function is "statically virtual", i.e.
				it simply calls the store function of a derived 
				class.
				See below why there is no non-const operator[] 
		*/
		const T& operator[](const CGLA::Vec3i& p) const 
		{
			return static_cast<const ChildT&>(*this)[p];
		}


		/** Store a voxel in grid. 

				This function returns nothing but is passed a 
				Vec3i p and T value t and stores t at p in the 
				grid. This function is "statically virtual", i.e.
				it simply calls the store function of a derived 
				class. 
				
				Yes, it would be simpler to provide a non-const
				operator[], however, a non-const operator[] will
				often be called even when no writing takes place.
				(Scott Meyers, "More Effective C++, p. 218)
				If a grid implementation allocates memory when 
				a voxel is accessed for writing, then it is a problem
				that we cannot be sure a non-const operator[] is 
				called only if we are writing. We might then allocate
				memory even if we just want to read. 
		*/
 		void store(const CGLA::Vec3i& p, const T& t)
		{
			return static_cast<ChildT&>(*this).store(p,t);
		}
	};
}

#endif
