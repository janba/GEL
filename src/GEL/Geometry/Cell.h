/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Cell.h
 * @brief Cell class for a voxel grid in a 2 level grid.
 */

#ifndef __GEOMETRY_CELL_H
#define __GEOMETRY_CELL_H

#include <vector>
#include "../CGLA/Vec3i.h"
#include "../CGLA/BitMask.h"
#include "RGrid.h"

namespace Geometry 
{
	/** \brief Class template for a cell in a hierarchical grid.

			The template arguments are the type and the cell dimension.
			A Cell is much like an RGrid - except that a Cell may 
			be coalesced into a single value and then split again
			when a new value is inserted. Also cells are constrained to
			be of size NxNxN where N=2^m for some m. This makes some things
			faster.

			The reason for making CELL_DIM a template argument rather 
			than a constructor argument is that we generally have many 
			cells	and do not wish that each cell contains its dimension. 
	*/
	template<class T, int CELL_DIM, class ChildT>
	class Cell
	{
	public:
		typedef T DataType;

	private:
		bool touched;

		/// Linear array containing data of this Cell.
		std::vector<T> data;
		
		/** Return a bit mask used to mask indices. 
				Some functions in this class are passed a Vec3i 
				whose lower bits index into the cell. The mask is
				used to select these lower bits before the 3D
				index is converted to	a linear index.
				
				If you are wondering about why the static variable 
				inside the function is not simply a static member of 
				the class, check out item 47 og "Effective C++" by 
				Scott Meyers.
		*/
		static const CGLA::BitMask& get_bit_mask()
		{
			static const CGLA::BitMask bit_mask(CELL_DIM);
			return bit_mask;
		}

	public:

		/** Get dimensions of Cell. Returns only one value since 
				x, y, and z dimensions are the same. */
		static int get_dim() 
		{
			return CELL_DIM;
		}

	private:
		
		/** Get index. Computes an index into the linear array representation
				of the voxels in the cell from a 3D grid index. The top bits (which
				have been used to find this Cell) are masked out. */
 		int get_idx(const CGLA::Vec3i& idx) const
		{
			CGLA::Vec3i bot_idx = get_bit_mask().mask(idx);
 			return (bot_idx[2]*get_dim()+bot_idx[1])*get_dim()+bot_idx[0];
		}

	protected:

		/** Store a value in the Cell. If the Cell is 
				coalesced, it is first resized to contain CELL_DIMS 
				cubed voxels, and then the voxel is inserted. */
 		void store_priv(const CGLA::Vec3i& p, const T& new_val) 
		{
			if(is_coalesced()) split();
			touched = true;
			data[get_idx(p)] = new_val;
		}

		
	public:
	
		/** Create empty grid cell. A Cell contains initially a single value.
				Reading (using operator[]) any voxel will yield this value. */
		Cell(const T& val): data(1,val), touched(false)
		{
			assert(CELL_DIM==get_dim());
		}

		/** Clear cell. Removes Cell contents replacing it with a single
				specified value. */
 		void coalesce(const T& val) 
		{
			data=std::vector<T>(1,val);
			touched = true;
		}

		/// Split cell - causes memory to be reserved for CELL_DIM^3 voxels.
		void split()
		{
			T val = data[0];
			data.resize(CGLA::qbe(CELL_DIM), val);
			touched = true;
		}

		/// Check if the cell is coalesced
		bool is_coalesced() const 
		{
			// We rely on size() being constant time below.
			// Probably this function should be changed.... it would be bad
			// if it was slow ... very bad....
			return data.size()==1;
		}

		
		/** Read access to voxel grid. This function is passed
				a Vec3i and returns the corresponding voxel. If the
				cell is coalesced, the single value stored is returned. */
		const T& operator[](const CGLA::Vec3i& p) const 
		{
			if(is_coalesced()) return data[0];
			return *(&data[get_idx(p)]);
		}

 		void store(const CGLA::Vec3i& p, const T& new_val) 
		{
			return static_cast<ChildT&>(*this).store(p,new_val);
		}
		
	
		/** Const get a pointer to the first element in the linear
				array representation of the grid. */
 		const T* get() const {return &data[0];}

		/** Non-const get a pointer to the first element in the linear
				array representation of the grid. */
		T* get() 
		{
			touched = true;
			return &data[0];
		}

		void untouch() {touched=false;}
		void touch() {touched=true;}
		bool is_touched() const {return touched;}

	};

	
	template<class T, int CELL_DIM>
	class DefaultCell: public Cell<T,CELL_DIM,DefaultCell<T, CELL_DIM> >
	{
	public:
		DefaultCell(const T& val): Cell<T,CELL_DIM,DefaultCell>(val){}
 		void store(const CGLA::Vec3i& p, const T& new_val) 
		{
			store_priv(p, new_val);
		}
		
	};


}
#endif
