/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file TrilinFilter.h
 * @brief Trilinear filter to apply to voxel grids.
 */

#ifndef __GEOMETRY_VOXELGRID_TRILINFILTER_H
#define __GEOMETRY_VOXELGRID_TRILINFILTER_H

#include "../CGLA/Vec3f.h"
#include "RGrid.h"
#include "Neighbours.h"

namespace Geometry
{
    /// Trilinear filter that can be applied to voxel grids.
	template<class GridT>
	class TrilinFilter
	{
	public:
		typedef typename GridT::DataType DataType;
	private:
		const GridT* grid;
	public:
		TrilinFilter(const GridT* _grid):	grid(_grid) {}

		bool in_domain(const CGLA::Vec3f&) const;

		float operator()(const CGLA::Vec3f& p) const;

		CGLA::Vec3f grad(const CGLA::Vec3f& v) const;
	};
}


#endif
