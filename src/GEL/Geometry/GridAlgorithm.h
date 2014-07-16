/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#ifndef __GEOMETRY_GRIDALGORITHM_H
#define __GEOMETRY_GRIDALGORITHM_H

/**
 * @file GridAlgorithm.h
 * A number of algorithms for traversing a voxel grid. The algorithms are 
 * generally transformative and allow us to invoke a function on each.
 */

/*
Functions:

  void for_each_voxel(Grid& g, F& f);
  void for_each_voxel(Grid& g,	const Vec3i& p0,	const Vec3i& p7, F& f);
  void for_each_voxel_ordered(Grid& g, const Vec3i& p0, const Vec3i& p7, F& f);
  void for_each_voxel_ordered(Grid& g,	F&);
  void for_each_voxel_const(const Grid& g, 
														const Vec3i& p0, const Vec3i& p7, F& f);
  void for_each_voxel_const(const Grid& g,	F& f);
  void for_each_voxel_ordered_const(Grid& g, 
																		const Vec3i& p0, const Vec3i& p7, 
																		F& f);
  void for_each_voxel_ordered_const(Grid& g, F& f);

Purpose:  
----------

Visit all voxels (or voxels in a region og interest) in grid. A function 
is invoked on each voxel. 

The "ordered" functions traverse the grid in a systematic way: A
slice at a time and for each slice one row at a time. The other
functions make no guarantees about how the volume (roi) is traversed. 

Types: 
----------

Grid - a grid type, either HGrid<T,CellT> or RGrid<T>

F - functor type a funtion or class with the function call operator 
overloaded. The functor should look either like this

void fun(const Vec3i&, T& x)

or this

void fun(const Vec3i&, const T& x)

in the case of the const functions.

Arguments:
----------

g  - The voxel grid, we wish to traverse.
p0 - (xmin, ymin, zmin) coordinates of the window we wish to
     traverse.
p7 - (xmax, ymax, zmax) coordinates of the window we wish to
     traverse.
f  - functor.

*/

#include <iostream>
#include "RGrid.h"
#include "HGrid.h"

namespace Geometry 
{
	template<class T, class F>
	void _for_each_voxel(T* data, 
											 int x_dim, int xy_dim,
											 const CGLA::Vec3i& p0, 
											 const CGLA::Vec3i& p7,
											 F& functor,
											 const CGLA::Vec3i& offset = CGLA::Vec3i(0))
	{
		const int Amin = p0[2]*xy_dim;
		const int Amax = p7[2]*xy_dim;
		int Bmin = Amin +p0[1]*x_dim;
		int Bmax = Amin +p7[1]*x_dim; 
		CGLA::Vec3i p0o = p0+offset;
		CGLA::Vec3i p(p0o);
		for(int A=Amin; A<Amax; A+=xy_dim, ++p[2])
			{
				p[1] = p0o[1];
				for(int B = Bmin; B<Bmax; B+=x_dim, ++p[1])
					{
						p[0] = p0o[0];
						int Cmin = B+p0[0];
						int Cmax = B+p7[0];
						for(int C=Cmin; C<Cmax; ++C, ++p[0])
							functor(p, data[C]);
					}
				Bmin += xy_dim;
				Bmax += xy_dim;
			}
	}

	template<class T, class F>
	void _for_each_voxel(T* data, 
											 const CGLA::Vec3i& dims,
											 F& functor,
											 const CGLA::Vec3i& offset = CGLA::Vec3i(0))
	{
		int l=0;
		CGLA::Vec3i p(offset);
		CGLA::Vec3i p7(offset);
		p7 += dims;

		for(;                   p[2]<p7[2]; ++p[2])
			for(p[1]=offset[1];   p[1]<p7[1]; ++p[1])
				for(p[0]=offset[0]; p[0]<p7[0]; ++p[0])
					functor(p, data[l++]);
	}


	/** Loop over all voxels in a sub-region (slice) of an
			RGrid and invoke a functor on each voxel. 
			The grid is the first argument, the slice is
			specified by the two subsequent args, and the functor 
			is the last argument. */
	template<class T, class F>
	void for_each_voxel(RGrid<T>& grid,
											const CGLA::Vec3i& p0, 
											const CGLA::Vec3i& p7,
											F& functor)
	{
		_for_each_voxel(grid.get(), grid.get_x_dim(), grid.get_xy_dim(),
										CGLA::v_max(p0, CGLA::Vec3i(0)),
										CGLA::v_min(p7, grid.get_dims()), 
										functor);
	}

	/** Loop over all voxels in an entire	RGrid. 
			Grid is the first argument, and a functor is 
			the second.	For each voxel, an operation 
			specified by the functor is performed. */
	template<class T, class F>
	void for_each_voxel(RGrid<T>& grid,	F& functor)
	{
		_for_each_voxel(grid.get(), grid.get_dims(), functor);
	}


	/** For each voxel (ordered). The idea of ordered traversal is 
			that we traverse the volume in a systematic fashion as opposed
			to traversing simply according to the memory layout of the volume
			data structure. This is important e.g. if we want to save the 
			volume in raw format. 
			For an RGrid, there is no difference though.
	*/ 
	template<class T, class F>
	void for_each_voxel_ordered(RGrid<T>& grid,
															const CGLA::Vec3i& p0, 
															const CGLA::Vec3i& p7,
															F& functor)
	{
		_for_each_voxel(grid.get(), grid.get_x_dim(), grid.get_xy_dim(),
										CGLA::v_max(p0, CGLA::Vec3i(0)),
										CGLA::v_min(p7, grid.get_dims()), 
										functor);
	}

	template<class T, class F>
	void for_each_voxel_ordered(RGrid<T>& grid,	F& functor)
	{
		_for_each_voxel(grid.get(), grid.get_dims(), functor);
	}


	template<class T, class CellT, class F>
	void for_each_cell(HGrid<T,CellT>& grid,
										 const CGLA::Vec3i& p0, 
										 const CGLA::Vec3i& p7,
										 F& functor)
	{
		CGLA::Vec3i p0t = p0/grid.get_bottom_dim();
		CGLA::Vec3i p7t = CGLA::v_min(p7/grid.get_bottom_dim()+
																	CGLA::Vec3i(1),
																	grid.get_top_dims());
		for(CGLA::Vec3i pt(p0t); pt[2]<p7t[2]; ++pt[2])
			for(pt[1]=p0t[1];      pt[1]<p7t[1]; ++pt[1])
				for(pt[0]=p0t[0];    pt[0]<p7t[0]; ++pt[0])
					functor(pt*CellT::get_dim(), grid.get_cell(pt));
	}

	template<class T, class CellT, class F>
	void for_each_cell(HGrid<T,CellT>& grid,
										 F& functor)
	{
		CGLA::Vec3i p0t;
		CGLA::Vec3i p7t =	grid.get_dims();
		const int inc = CellT::get_dim();
		int l=0;
		for(CGLA::Vec3i pt(p0t); pt[2]<p7t[2]; pt[2]+=inc)
			for(pt[1]=p0t[1];      pt[1]<p7t[1]; pt[1]+=inc)
				for(pt[0]=p0t[0];    pt[0]<p7t[0]; pt[0]+=inc)
					functor(pt, grid.get_cell(l++));
	}

	
	template<class CellT, class F>
	class _HGridCellFunctor
	{
		const CGLA::Vec3i p0;
		const CGLA::Vec3i p7;
		F& functor;

	public:
		_HGridCellFunctor(const CGLA::Vec3i _p0,
											const CGLA::Vec3i _p7,
											F& _functor): p0(_p0), p7(_p7), functor(_functor) {}
			
		void operator()(const CGLA::Vec3i& offset, 
										CellT& cell)
		{
			CGLA::Vec3i p0c = CGLA::v_max(p0-offset, CGLA::Vec3i(0));
			CGLA::Vec3i p7c = CGLA::v_min(p7-offset, CGLA::Vec3i(CellT::get_dim()));

			if(cell.is_coalesced())
				cell.split();

			_for_each_voxel(cell.get(), 
											CellT::get_dim(), 
											CGLA::sqr(CellT::get_dim()),
											p0c, p7c,	functor, offset);
		}
	};


	template<class T, class CellT, class F>
	void for_each_voxel(HGrid<T,CellT>& grid,
											const CGLA::Vec3i& _p0, 
											const CGLA::Vec3i& _p7,
											F& functor)
	{
		CGLA::Vec3i p0 = CGLA::v_max(_p0, CGLA::Vec3i(0));
		CGLA::Vec3i p7 = CGLA::v_min(_p7, grid.get_dims());
		_HGridCellFunctor<CellT,F> cell_functor(p0, p7, functor);
		for_each_cell(grid, p0, p7, cell_functor);
	}

	template<class T, class CellT, class F>
	void for_each_voxel(HGrid<T,CellT>& grid,	F& functor)
	{
		_HGridCellFunctor<CellT,F> cell_functor(CGLA::Vec3i(0), 
																						grid.get_dims(), functor);
		for_each_cell(grid, cell_functor);
	}

	template<class T, class CellT, class F>
	void for_each_voxel_ordered(HGrid<T,CellT>& grid,
															const CGLA::Vec3i& _p0, 
															const CGLA::Vec3i& _p7,
															F& functor)
	{
		CGLA::Vec3i p0 = CGLA::v_max(_p0, CGLA::Vec3i(0));
		CGLA::Vec3i p7 = CGLA::v_min(_p7, grid.get_dims());
		for(int k=p0[2];k<p7[2];++k)
			for(int j=p0[1];j<p7[1];++j)
				for(int i=p0[0];i<p7[0];++i)
					{
						CGLA::Vec3i p(i,j,k);
						float val = grid[p];
						float nval = val;
						functor(p, nval);
						if(nval != val)
							grid.store(p, nval);
					}
	}

	template<class T, class CellT, class F>
	void for_each_voxel_ordered(HGrid<T,CellT>& grid,	F& functor)
	{
		for_each_voxel_ordered(grid, CGLA::Vec3i(0), grid.get_dims(), functor);
	}


	template<typename T>
	class _AssignFun
	{
		T val;
	public:
		_AssignFun(const T& _val): val(_val) {}
		void operator()(const CGLA::Vec3i& pi, T& vox_val)
		{
			vox_val = val;
		}
	};
	
 	template<class G>
	void clear_region(G& grid, const typename G::DataType& value)
	{
		_AssignFun<typename G::DataType> afun(value);
		for_each_voxel(grid, afun);
	}

	template<class G>
	void clear_region(G& grid, 
										const CGLA::Vec3i& p0,
										const CGLA::Vec3i& p7,
										const typename G::DataType& value)
	{
		_AssignFun<typename G::DataType>afun(value) ;
		for_each_voxel(grid, p0, p7, afun);
	}


	//----------------------------------------------------------------------
	// const versions.

	/** Loop over all voxels in a sub-region (slice) of an
			RGrid and invoke a functor on each voxel. 
			The grid is the first argument, the slice is
			specified by the two subsequent args, and the functor 
			is the last argument. */
	template<class T, class F>
	void for_each_voxel_const(const RGrid<T>& grid,
														const CGLA::Vec3i& p0, 
														const CGLA::Vec3i& p7,
														F& functor)
	{
		_for_each_voxel(grid.get(), grid.get_x_dim(), grid.get_xy_dim(),
										CGLA::v_max(p0, CGLA::Vec3i(0)),
										CGLA::v_min(p7, grid.get_dims()), 
										functor);
	}

	/** Loop over all voxels in an entire	RGrid. 
			Grid is the first argument, and a functor is 
			the second.	For each voxel, an operation 
			specified by the functor is performed. */
	template<class T, class F>
	void for_each_voxel_const(const RGrid<T>& grid,	F& functor)
	{
		_for_each_voxel(grid.get(), grid.get_dims(), functor);
	}


	/** For each voxel (ordered). The idea of ordered traversal is 
			that we traverse the volume in a systematic fashion as opposed
			to traversing simply according to the memory layout of the volume
			data structure. This is important e.g. if we want to save the 
			volume in raw format. 
			For an RGrid, there is no difference though.
	*/ 
	template<class T, class F>
	void for_each_voxel_ordered_const(const RGrid<T>& grid,
																		const CGLA::Vec3i& p0, 
																		const CGLA::Vec3i& p7,
																		F& functor)
	{
		_for_each_voxel(grid.get(), grid.get_x_dim(), grid.get_xy_dim(),
										CGLA::v_max(p0, CGLA::Vec3i(0)),
										CGLA::v_min(p7, grid.get_dims()), 
										functor);
	}

	template<class T, class F>
	void for_each_voxel_ordered_const(const RGrid<T>& grid,	F& functor)
	{
		_for_each_voxel(grid.get(), grid.get_dims(), functor);
	}


	template<class T, class CellT, class F>
	void for_each_cell_const(const HGrid<T,CellT>& grid,
													 const CGLA::Vec3i& p0, 
													 const CGLA::Vec3i& p7,
													 F& functor)
	{
		CGLA::Vec3i p0t = p0/grid.get_bottom_dim();
		CGLA::Vec3i p7t = CGLA::v_min(p7/grid.get_bottom_dim()+
																	CGLA::Vec3i(1),
																	grid.get_top_dims());
		for(CGLA::Vec3i pt(p0t); pt[2]<p7t[2]; ++pt[2])
			for(pt[1]=p0t[1];      pt[1]<p7t[1]; ++pt[1])
				for(pt[0]=p0t[0];    pt[0]<p7t[0]; ++pt[0])
					functor(pt*CellT::get_dim(), grid.get_cell(pt));
	}

	template<class T, class CellT, class F>
	void for_each_cell_const(const HGrid<T,CellT>& grid,
													 F& functor)
	{
		CGLA::Vec3i p0t;
		CGLA::Vec3i p7t =	grid.get_dims();
		const int inc = CellT::get_dim();
		int l=0;
		for(CGLA::Vec3i pt(p0t); pt[2]<p7t[2]; pt[2]+=inc)
			for(pt[1]=p0t[1];      pt[1]<p7t[1]; pt[1]+=inc)
				for(pt[0]=p0t[0];    pt[0]<p7t[0]; pt[0]+=inc)
					functor(pt, grid.get_cell(l++));
	}

	
	template<class CellT, class F>
	class _HGridCellFunctorConst
	{
		const CGLA::Vec3i p0;
		const CGLA::Vec3i p7;
		F& functor;

	public:
		_HGridCellFunctorConst(const CGLA::Vec3i _p0,
													 const CGLA::Vec3i _p7,
													 F& _functor): p0(_p0), p7(_p7), functor(_functor) {}
			
		void operator()(const CGLA::Vec3i& offset, 
										const CellT& cell)
		{
			CGLA::Vec3i p0c = CGLA::v_max(p0-offset, CGLA::Vec3i(0));
			CGLA::Vec3i p7c = CGLA::v_min(p7-offset, CGLA::Vec3i(CellT::get_dim()));

			if(cell.is_coalesced())
				{
					typename CellT::DataType val = *cell.get();
					for(CGLA::Vec3i p(p0c); p[2]<p7c[2]; ++p[2])
						for(p[1]=p0c[1];      p[1]<p7c[1]; ++p[1])
							for(p[0]=p0c[0];    p[0]<p7c[0]; ++p[0])
								functor(p+offset, val);
				}
			_for_each_voxel(cell.get(), 
											CellT::get_dim(), 
											CGLA::sqr(CellT::get_dim()),
											p0c, p7c,	functor, offset);
		}
	};


	template<class T, class CellT, class F>
	void for_each_voxel_const(const HGrid<T,CellT>& grid,
														const CGLA::Vec3i& _p0, 
														const CGLA::Vec3i& _p7,
														F& functor)
	{
		CGLA::Vec3i p0 = CGLA::v_max(_p0, CGLA::Vec3i(0));
		CGLA::Vec3i p7 = CGLA::v_min(_p7, grid.get_dims());
		_HGridCellFunctorConst<CellT,F> cell_functor(p0, p7, functor);
		for_each_cell_const(grid, p0, p7, cell_functor);
	}

	template<class T, class CellT, class F>
	void for_each_voxel_const(const HGrid<T,CellT>& grid,	F& functor)
	{
		_HGridCellFunctorConst<CellT,F> cell_functor(CGLA::Vec3i(0), 
																								 grid.get_dims(), functor);
		for_each_cell_const(grid, cell_functor);
	}

	template<class T, class CellT, class F>
	void for_each_voxel_ordered_const(const HGrid<T,CellT>& grid,
																		const CGLA::Vec3i& _p0, 
																		const CGLA::Vec3i& _p7,
																		F& functor)
	{
		CGLA::Vec3i p0 = CGLA::v_max(_p0, CGLA::Vec3i(0));
		CGLA::Vec3i p7 = CGLA::v_min(_p7, grid.get_dims());
		for(int k=p0[2];k<p7[2];++k)
			for(int j=p0[1];j<p7[1];++j)
				for(int i=p0[0];i<p7[0];++i)
					{
						CGLA::Vec3i p(i,j,k);
						functor(p, grid[p]);
					}
	}

	template<class T, class CellT, class F>
	void for_each_voxel_ordered_const(const HGrid<T,CellT>& grid,	F& functor)
	{
		for_each_voxel_ordered_const(grid, 
																 CGLA::Vec3i(0), 
																 grid.get_dims(), 
																 functor);
	}


}

#endif
