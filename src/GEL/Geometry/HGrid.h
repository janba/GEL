/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file HGrid.h
 * @brief Hierarchical voxel grid - space saving.
 */

#ifndef __GEOMETRY_HGRID_H
#define __GEOMETRY_HGRID_H
// Author: J. Andreas Bærentzen,
// Created: Wed Jan 24 18:29:0

#include <vector>
#include "AncestorGrid.h"
#include "Cell.h"

namespace Geometry 
{
	/** \brief Hierarchical voxel grid. 

			In many cases we wish to save on the storage requirements of volumes.
			A hierarchical voxel grid is a volume representation where the volume 
			is divided into box shaped regions, and each region is represented only 
			if it contains voxels. This class template is for such a grid.
	*/
	template<class T, class CellT=DefaultCell<T,8> >
	class HGrid: public AncestorGrid<T,HGrid<T,CellT> > 
	{
	public:
		typedef T DataType;		
		typedef CellT CellType;

	private:

		/// Dimensions of top level grid.
		const CGLA::Vec3i top_dims;

		/// Top level grid. I.e. vector of sub grids.
		std::vector<CellT> top_grid;

		/// The default grid value, used to clear grid. 
		DataType default_val;
		
		/// Size of the top grid (number of cells x*y*z)
		int top_grid_size;

	public:

		const CGLA::Vec3i& get_top_dims() const {return top_dims;}
		
		int get_bottom_dim() const {return CellT::get_dim();}

		
	private:
 
		/// Get index into top level grid from int vector position.
		int get_top_index(const CGLA::Vec3i& idx) const
		{
			const CGLA::Vec3i top_idx = idx/get_bottom_dim();
			return (top_idx[2]*top_dims[1]+top_idx[1])*top_dims[0]+top_idx[0];
		}

	public:

		/// Construct grid of specified dimensions
		HGrid(const CGLA::Vec3i& dims, const T& val = T()):
			AncestorGrid<T,HGrid<T,CellT> >(dims),
			top_dims(dims/CellT::get_dim()+
							 CGLA::Vec3i(dims[0]%CellT::get_dim()?1:0,
													 dims[1]%CellT::get_dim()?1:0,
													 dims[2]%CellT::get_dim()?1:0)),
			top_grid(top_dims[0]*top_dims[1]*top_dims[2],CellT(val)),
			default_val(val), top_grid_size(top_grid.size())
		{}

		/** Store a voxel vox at position p in grid. The Cell will
				automatically subdivide if it is not already subdivided. */
 		void store(const CGLA::Vec3i& p, const T& vox)
		{
			assert(this->in_domain(p));
			top_grid[get_top_index(p)].store(p, vox);
		}


		/** Read only access to a voxel in the grid. */
 		const T& operator[](const CGLA::Vec3i& p) const
		{
			assert(this->in_domain(p));
			return top_grid[get_top_index(p)][p];
		}


//  		bool get(const CGLA::Vec3i& p, T* voxel)
// 		{
// 			assert(in_domain(p));
// 			CellT* cell = top_grid[get_top_index(p)];
// 			if(cell->is_coalesced()) return false;
// 			voxel = cell[p];
// 		}

		CellT& get_cell(const CGLA::Vec3i& p)
		{
			return top_grid[(p[2]*top_dims[1]+p[1])*top_dims[0]+p[0]];
		}

		const CellT& get_cell(const CGLA::Vec3i& p) const
		{
			return top_grid[(p[2]*top_dims[1]+p[1])*top_dims[0]+p[0]];
		}

		CellT& get_cell(int i) 
		{
			return top_grid[i];
		}

		const CellT& get_cell(int i) const
		{
			return top_grid[i];
		}

		void clear()
		{
			int N = top_grid.size();
			for(int i=0;i<N;++i)
				top_grid[i].coalesce(default_val);
		}
		
	};

}
#endif
