/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file BSPTree.h
 * @brief Binary space partitioning tree.
 */

#ifndef __GEOMETRY_BSPTREE_H__
#define __GEOMETRY_BSPTREE_H__

#include <vector>

#include "../CGLA/Mat4x4f.h"
#include "../Geometry/TriMesh.h"

#include "BBox.h"

namespace Geometry 
{
  struct BSPNode;

  // Initialization structure
  struct BSPNode 
  {
    unsigned char axis_leaf; // 00 = axis 0, 01 = axis 1, 10 = axis 2, 11 = leaf
    double plane;
    BSPNode *left, *right;
    size_t id;
    int count;
    unsigned int ref;
  };

  // Faster structure
  struct BSPLeaf {
    // bits 0..30 : offset to first son
    // bit 31 (sign) flag whether node is a leaf
    unsigned int flagAndOffset;
    unsigned int count;
  };

  struct BSPInner {
    unsigned int flagAndOffset;
    // bits 0..1 : splitting dimension
    // bits 2..30: offset bits
    // bit 31 (sign) flag whether node is a leaf
    float splitCoordinate;
  };

  union FastBSPNode
  {
    BSPLeaf leaf;
    BSPInner inner;
  };

    /** BSPTree class. Mostly used to accelerate ray casting. 
     */
  class BSPTree 
  {
    bool b_is_build;

	  static int node_calls;
    static int tri_calls;
    BSPNode* root;
    BBox bbox;
    std::vector<const Geometry::TriMesh*> trimesh;
    std::vector<CGLA::Mat4x4f> transforms;
    
    std::vector<ISectTri> isecttris;
    std::vector<TriAccel> triaccel;
    std::vector<ISectTri*> all_objects;
    std::vector<TriAccel*> all_triaccel;
    std::vector<FastBSPNode> fast_tree;
    
    unsigned int max_objects;
    unsigned int max_level;
	
  public:
    BSPTree();
    ~BSPTree();

    void init(std::vector<const Geometry::TriMesh*>& trimesh, 
              int max_objects, int max_level);
    void init(std::vector<const Geometry::TriMesh*>& _trimesh, 
              std::vector<CGLA::Mat4x4f>& _transforms, 
              int _max_objects, int _max_level);
    void init(const Geometry::TriMesh* mesh, CGLA::Mat4x4f transform, 
              std::vector<int> &trilist, 
              int _max_objects, int _max_level);

    void build();
    bool is_build();
    void clear();

    bool intersect(Ray &ray) const;
	
  private:
    void delete_node(BSPNode *node);
    void subdivide_node(BSPNode &node, BBox &bbox, 
                        unsigned int level, 
                        std::vector<ISectTri*>& objects, 
                        std::vector<TriAccel*>& tri_objects);
    void init();

    bool intersect_node(Ray &ray, const BSPNode &node, 
                        double t_min, double t_max) const;
    void print(BSPNode *node, int depth);
    int size(BSPNode *node);
    int size();

    void make_fast_tree(BSPNode *node);
    void push_fast_bsp_node(BSPNode *node, int id);
    void intersect_fast_node(Ray &ray, 
                             const FastBSPNode *node, 
                             double t_min, double t_max) const;
    bool intersect(Ray &ray, const ISectTri &isecttri, double t_max) const;
  };
}
#endif

