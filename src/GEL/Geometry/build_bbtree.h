/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file build_bbtree.h
 * @brief Contains functions for building a bounding box tree.
 */

#ifndef __GEOMETRY_BUILD_BBTREE_H
#define __GEOMETRY_BUILD_BBTREE_H

#include "BoundingTree.h"

namespace HMesh
{
    class Manifold;
}

namespace Geometry
{
typedef BoundingTree<OBox> OBBTree;
typedef BoundingTree<AABox> AABBTree;

void build_OBBTree(HMesh::Manifold& m, OBBTree& tree);
void build_AABBTree(HMesh::Manifold& m, AABBTree& tree);

}
#endif
