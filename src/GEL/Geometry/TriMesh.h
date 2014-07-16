/* ----------------------------------------------------------------------- *
* This file is part of GEL, http://www.imm.dtu.dk/GEL
* Copyright (C) the authors and DTU Informatics
* For license and list of authors, see ../../doc/intro.pdf
* ----------------------------------------------------------------------- */

/**
* @file TriMesh.h
* @brief Triangle mesh class based on indexed face set data structure.
*/

#ifndef __GEOMETRY_TRIMESH_H__
#define __GEOMETRY_TRIMESH_H__

#include "../CGLA/Mat4x4f.h"
#include "IndexedFaceSet.h"
#include "Material.h"

namespace Geometry
{
  /** \brief A Triangle Mesh struct. 

      This struct contains three indexed face sets, one for geometry,
      textures, and normals. It also contains a vector of materials and a
      vector of texture maps.

      A number of functions are defined allowing for rendering of the triangle
      mesh. */
  class TriMesh 
  {
  public:
		
    // ------- DATA -------

    /// Name of model
    std::string name;
		
    /// Indexed face set for the actual geometry
    IndexedFaceSet geometry;

    /// Indexed face set for the normals
    IndexedFaceSet normals;

    /// Indexed face set for the texture coordinates.
    IndexedFaceSet texcoords;

    /// Material indices for all faces
    std::vector<int> mat_idx;

    /// Texture indices for all faces
    std::vector<int> tex_idx;

    /// Vector of materials
    std::vector<Material> materials;

    /// Vector of triangle face areas
    std::vector<float> face_areas;

    /// Tabulated cumulative distribution function for sampling triangles according to area
    std::vector<float> face_area_cdf;

    /// Total surface area of the triangle mesh
    float surface_area;

    // -------- FUNCTIONS -----------

    /// Get an axis aligned bounding box for the model.
    bool get_bbox(CGLA::Vec3f& p0, CGLA::Vec3f& p7) const;

    /// Get a bounding sphere for the model.
    bool get_bsphere(CGLA::Vec3f& c, float& r) const;

    /// Returns true if at least one normal has been defined.
    bool has_normals() const { return normals.no_faces() > 0; }

    /// Find a material from its name
    int find_material(const std::string&) const;

    /// Compute normals for the mesh. Does not check if there are normals.
    void compute_normals();

    /// Compute areas for all faces and total surface area.
    void compute_areas();

    /// Apply a transformation matrix to the mesh
    void transform(const CGLA::Mat4x4f& m);

    /// Apply a transformation matrix to the texture coordinates
    void tex_transform(const CGLA::Mat4x4f& m);

    /// Apply a transformation matrix to the texture coordinates associated with a particular material
    void tex_transform(const CGLA::Mat4x4f& m, const std::string& material);
  };
}
#endif
