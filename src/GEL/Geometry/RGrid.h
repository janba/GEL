/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 @file RGrid.h
 Regular voxel grid - just an array of voxels. This is a template since we
 might want different types of voxel grids.
 */

#ifndef __GEOMETRY_RGRID_H
#define __GEOMETRY_RGRID_H

#include <vector>
#include "../CGLA/Vec3i.h"
#include "AncestorGrid.h"

namespace Geometry 
{

	/** \brief Regular voxel grid. 

			This class template can be instantiated and used directly.
			This is just a regular voxel grid stored as a linear array
			with functions to access its contents. */			
	template<class T>
	class RGrid: public AncestorGrid<T,RGrid<T> >
	{
	public:
		typedef T DataType;

	private:
		/// x size of grid.
		int x_dim;
		
		/// x size times y size of grid. Stored for efficiency reasons.
		int xy_dim;

		/// Vector containing the actual data.
		std::vector<T> data;

		/// Convert xyz index into an index in a linear array.
		int grid_idx(const CGLA::Vec3i& idx) const
		{
			return  idx[2] * xy_dim + idx[1] * x_dim + idx[0];
		}

		/// The default grid value, used to clear grid. 
		DataType default_val;

	public:
	
		/** Construct a regular voxel grid. This function
				is passed a Vec3i _dims and an optional 
				initialization value, val. It creates a grid
				of specified dimensions, and initializes the 
				value of all voxels to val. */
		RGrid(CGLA::Vec3i _dims, const T& val = T()):
			AncestorGrid<T,RGrid<T> >(_dims), 
			x_dim(_dims[0]), xy_dim(_dims[0]*_dims[1]),
			data(_dims[0]*_dims[1]*_dims[2],val),
			default_val(val)
		{}	

		/** Construct a grid of dimensions 0,0,0 */
		RGrid(): AncestorGrid<T,RGrid<T> >(CGLA::Vec3i(0)), 
			x_dim(0), xy_dim(0),
			data(0,0), default_val(0)
		{}	

		/** Store a voxel in a regular grid. */
		void store(const CGLA::Vec3i& p, const T& t) 
		{
			assert(this->in_domain(p));
			data[grid_idx(p)] = t;
		}

		/** Read/write access to voxel grid. This is 
				a non-const overloaded operator[]. In a regular 
				grid, we have reserved room for all voxels in 
				advance. Hence, it is possible to create a non-const
				operator[]. See AncestorGrid::operator[]. */
		T& operator[](const CGLA::Vec3i& p) 
		{
			assert(this->in_domain(p));
			return data[grid_idx(p)];
		}

		/// Read only access to voxel grid through const operator[]
		const T& operator[](const CGLA::Vec3i& p) const 
		{
			assert(this->in_domain(p));
			return data[grid_idx(p)];
		}

		/// Const function to get a pointer to the first voxel in grid.
		const T* get() const {return &data[0];}

		/// Non-const function to get a pointer to the first voxel in grid.
		T* get() {return &data[0];}

		/// Get x dimensions of grid.
		int get_x_dim() const { return x_dim;}

		/// Get x size times y size of grid.
		int get_xy_dim() const { return xy_dim;}

		/// Get length of linear array actually containing voxels.
		int get_size() const { return data.size();}

		void clear()
		{
			int N = data.size();
			for(int i=0;i<N;++i)
				data[i] = default_val;
		}


	};

}

#endif
