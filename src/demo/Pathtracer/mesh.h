#ifndef __PATHTRACER_MESH__H__
#define __PATHTRACER_MESH__H__

#include "luminaire.h"
#include "material.h"

namespace Geometry { class TriMesh; class IndexedFaceSet; }

class mesh : public luminaire
{
public:
    //ctor
    mesh(void);
    mesh(const std::string&);
    ~mesh(void);

    //geometry
    const Geometry::TriMesh& get_trimesh(void) const;
    void set_trimesh(const Geometry::TriMesh&);
    const CGLA::Vec3f& triangle_normal(size_t i) const;

    //material
    const material* get_material(void) const;
    void set_material(const material*);

    //exitance (W/m2) (assuming diffuse emission)
    CGLA::Vec3f exitance(void) const;
    void set_exitance(const CGLA::Vec3f&);

    //sampling of mesh luminaire
    bool sample(const ray&, const hit_info&, CGLA::Vec3f& L, 
        CGLA::Vec3f& w) const;

protected:
    Geometry::TriMesh* msh_;
    std::vector<CGLA::Vec3f> triangle_normals_;
    Geometry::IndexedFaceSet* tangents_;
    const material* mat_;
    CGLA::Vec3f exitance_;
};

#endif

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
