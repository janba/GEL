/* ----------------------------------------------------------------------- *
* This file is part of GEL, http://www.imm.dtu.dk/GEL
* Copyright (C) the authors and DTU Informatics
* For license and list of authors, see ../../doc/intro.pdf
* ----------------------------------------------------------------------- */

#include "../CGLA/Vec3f.h"
#include <stdio.h>
#include <iostream>
#include <set>

#include "TriMesh.h"

using namespace std;
using namespace CGLA;

namespace Geometry 
{
  int TriMesh::find_material(const string& name) const
  {
    for(size_t i = 0; i < materials.size(); ++i)
      if(materials[i].name == name)
        return i;
    return 0;
  }
	
  void TriMesh::compute_normals()
  {
    // By default the normal faces are the same as the geometry faces
    // and there are just as many normals as vertices, so we simply
    // copy.
    normals = geometry;

    const int NV = normals.no_vertices();
    // The normals are initialized to zero.
    int i;
    for(i=0;i<NV; ++i)
      normals.vertex_rw(i) = Vec3f(0);

    // For each face
    int NF = geometry.no_faces();
    for(i=0;i<NF; ++i)
    {
      // Compute the normal
      const Vec3i& f  = geometry.face(i);
      const Vec3f& p0 = geometry.vertex(f[0]);
      const Vec3f& a  = geometry.vertex(f[1]) - p0;
      const Vec3f& b  = geometry.vertex(f[2]) - p0;
      Vec3f face_normal = cross(a,b);
      float l = sqr_length(face_normal);
      if(l > 0.0f)
        face_normal /= sqrt(l);

      // Add the angle weighted normal to each vertex
      for(int j=0;j<3; ++j)
      {
        const Vec3f& p0 = geometry.vertex(f[j]);
        Vec3f a = geometry.vertex(f[(j+1)%3]) - p0;
        float l_a = sqr_length(a);
        if(l_a > 0.0f)
          a /= sqrt(l_a);
        Vec3f b = geometry.vertex(f[(j+2)%3]) - p0;
        float l_b = sqr_length(b);
        if(l_b > 0.0f)
          b /= sqrt(l_b);
        float d = max(-1.0f, min(1.0f, dot(a,b)));
        normals.vertex_rw(f[j]) += face_normal * acos(d);
      }
    }

    // Normalize all normals
    for(i=0;i<NV; ++i)
    {
      float l_vert_rw = sqr_length(normals.vertex_rw(i));
      if(l_vert_rw > 0.0f)
        normals.vertex_rw(i) /= sqrt(l_vert_rw);
    }
  }

  void TriMesh::compute_areas()
  {
    int no_of_faces = geometry.no_faces();
    surface_area = 0.0f;
    face_areas.resize(no_of_faces);
    face_area_cdf.resize(no_of_faces);
    for(int i = 0; i < no_of_faces; ++i)
    {
      const Vec3i& f  = geometry.face(i);
      const Vec3f& p0 = geometry.vertex(f[0]);
      const Vec3f& a  = geometry.vertex(f[1]) - p0;
      const Vec3f& b  = geometry.vertex(f[2]) - p0;
      face_areas[i] = 0.5f*cross(a, b).length();
      face_area_cdf[i] = surface_area + face_areas[i];
      surface_area += face_areas[i];
    }
    if(surface_area > 0.0f)
      for(int i = 0; i < no_of_faces; ++i)
        face_area_cdf[i] /= surface_area;
  }

  bool TriMesh::get_bbox(Vec3f& p0, Vec3f& p7) const
  {
    if(geometry.no_vertices() == 0)
      return false;

    int i;
    p0 = geometry.vertex(0);
    p7 = geometry.vertex(0);
    for(i = 1; i < geometry.no_vertices(); ++i) 
    {
      p0 = v_min(geometry.vertex(i), p0);
      p7 = v_max(geometry.vertex(i), p7);
    }
    return true;
  }

  bool TriMesh::get_bsphere(Vec3f& c, float& r) const
  {
    Vec3f p0,p7;
    if(!get_bbox(p0, p7))
      return false;

    Vec3f rad = (p7 - p0)/2.0;
    c = p0 + rad;
    r = rad.length();
    return true;
  }

  void TriMesh::transform(const Mat4x4f& m)
  {
    for(int i = 0; i < geometry.no_vertices(); ++i)
      geometry.vertex_rw(i) = m.mul_3D_point(geometry.vertex(i));
    for(int i = 0; i < normals.no_vertices(); ++i)
      normals.vertex_rw(i) = normalize(m.mul_3D_vector(normals.vertex(i)));
  }

  void TriMesh::tex_transform(const Mat4x4f& m)
  {
    for(int i = 0; i < texcoords.no_vertices(); ++i)
      texcoords.vertex_rw(i) = m.mul_3D_point(texcoords.vertex(i));
  }

  void TriMesh::tex_transform(const Mat4x4f& m, const string& material)
  {
    set<int> v_touched;
    int m_idx = find_material(material);
    for(int i = 0; i < texcoords.no_faces(); ++i)
      if(mat_idx[i] == m_idx)
      {
        Vec3i face = texcoords.face(i);
        for(int j = 0; j < 3; ++j)
        {
          int v_idx = face[j];
          if(v_touched.count(v_idx) == 0)
          {
            v_touched.insert(v_idx);
            texcoords.vertex_rw(v_idx) = m.mul_3D_point(texcoords.vertex(v_idx));
          }
        }
      }
  }
}
