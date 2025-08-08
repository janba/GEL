/**
 * @file load_raw.h
 * @brief contains function template for loading raw voxel data.
 */

/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#ifndef GEOMETRY_VOXELGRID_LOAD_RAW_H
#define GEOMETRY_VOXELGRID_LOAD_RAW_H

#include <string>
#include <GEL/CGLA/Mat.h>
#include <GEL/Geometry/RGrid.h>

namespace Geometry
{
    /// Function template for loading raw voxel data. Template arg is the grid type.
	template<class T>
		bool load_raw(const std::string& file, RGrid<T>& grid);
	
}

#endif
