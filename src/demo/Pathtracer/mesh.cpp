#include "mesh.h"
#include "scene.h"

#include <GEL/Geometry/obj_load.h>

using namespace CGLA;
using namespace Geometry;

extern scene* current;

mesh::mesh(void)
{
    msh_ = new TriMesh;
    tangents_ = new IndexedFaceSet;
    mat_ = 0;
    exitance_ = Vec3f(0.f);
}

mesh::mesh(const std::string& fn)
{
    msh_ = new TriMesh;
    tangents_ = new IndexedFaceSet;
    mat_ = 0;
    exitance_ = Vec3f(0.f);

    TriMesh msh;
    obj_load(fn, msh);
    set_trimesh(msh);
}

mesh::~mesh(void)
{
    delete tangents_;
    delete msh_;
}

const CGLA::Vec3f& mesh::triangle_normal(size_t i) const
{
    return triangle_normals_.at(i);
}

const TriMesh& mesh::get_trimesh(void) const
{
    return *msh_;
}

void mesh::set_trimesh(const TriMesh& tm)
{
    *msh_ = tm;

    msh_->compute_normals();
    //
    //int tris = msh_->geometry.no_faces();

    //std::vector<Vec3f> tan1;
    //tan1.resize(tris);
    //triangle_normals_.resize(tris);

    //for (int t=0; t<tris; t++)
    //{
    //  const Vec3i& f  = msh_->geometry.face(t);
    //  const Vec3f p0 = msh_->geometry.vertex(f[0]);
    //  const Vec3f a  = msh_->geometry.vertex(f[1]) - p0;
    //  const Vec3f b  = msh_->geometry.vertex(f[2]) - p0;
    //  Vec3f tri_normal = cross(a,b);
    //  float l = sqr_length(tri_normal);
    //  if(l > 0.0f)
    //      tri_normal /= sqrt(l);
    //  triangle_normals_[t] = trs_.mul_3D_vector(tri_normal);

    //}
#if 0
    if (msh_->texcoords.no_faces() == 0)
        return;

    for (long a = 0; a < triangleCount; a++)
    {
        long i1 = triangle->index[0];
        long i2 = triangle->index[1];
        long i3 = triangle->index[2];

        const Point3D& v1 = vertex[i1];
        const Point3D& v2 = vertex[i2];
        const Point3D& v3 = vertex[i3];

        const Point2D& w1 = texcoord[i1];
        const Point2D& w2 = texcoord[i2];
        const Point2D& w3 = texcoord[i3];

        float x1 = v2.x - v1.x;
        float x2 = v3.x - v1.x;
        float y1 = v2.y - v1.y;
        float y2 = v3.y - v1.y;
        float z1 = v2.z - v1.z;
        float z2 = v3.z - v1.z;

        float s1 = w2.x - w1.x;
        float s2 = w3.x - w1.x;
        float t1 = w2.y - w1.y;
        float t2 = w3.y - w1.y;

        float r = 1.0F / (s1 * t2 - s2 * t1);
        Vector3D sdir((t2 * x1 - t1 * x2) * r, (t2 * y1 - t1 * y2) * r,
            (t2 * z1 - t1 * z2) * r);
        Vector3D tdir((s1 * x2 - s2 * x1) * r, (s1 * y2 - s2 * y1) * r,
            (s1 * z2 - s2 * z1) * r);

        tan1[i1] += sdir;
        tan1[i2] += sdir;
        tan1[i3] += sdir;

        tan2[i1] += tdir;
        tan2[i2] += tdir;
        tan2[i3] += tdir;

        triangle++;
    }
#endif
}

const material* mesh::get_material(void) const
{
    return mat_;
}

void mesh::set_material(const material* m)
{
    mat_ = m;
}

bool mesh::sample(const ray& r,
                  const hit_info& hi,
                  Vec3f& Li,
                  Vec3f& w) const
{
    //skip non-emitters
    if (intensity(exitance_) == 0.f)
        return false;

    //pick a random face
    int face = gel_rand() % msh_->geometry.no_faces();
    float pface = 1.f / float(msh_->geometry.no_faces());

    //pick a random position in triangle
    float r1 = std::sqrt(mt_random());
    float r2 = mt_random();
    Vec3f barycentric(1.f - r1, (1.f - r2) * r1, r1 * r2);

    //find position
    Vec3i fidx = msh_->geometry.face(face);
    Vec3f vertices[3];
    Vec3f light_pos(0.f);
    for (int i=0; i<3; ++i)
    {
        vertices[i] = msh_->geometry.vertex(fidx[i]);
        light_pos += barycentric[i] * msh_->geometry.vertex(fidx[i]);
    }

    light_pos = trs_.mul_3D_point(light_pos);

    float parea = 2.f / length(cross(vertices[2]-vertices[1], 
        vertices[0]-vertices[1]));

    //find normal
    Vec3i nidx = msh_->normals.face(face);
    Vec3f light_normal(0.f);
    for (int i=0; i<3; ++i)
        light_normal += barycentric[i] * msh_->normals.vertex(nidx[i]);
    light_normal = trs_.mul_3D_vector(light_normal);
    light_normal.normalize();

    //find w
    w = normalize(light_pos - hi.position);
    //float costheta1 = dot(hi.shading_normal, w);
    //if (costheta1 < 0.f)
    //  return false;

    float costheta2 = dot(light_normal, -w);
    if (costheta2 < 0.f)
        return false;

    float d2 = sqr_length(light_pos - hi.position);

    //trace shadow feeler..
    ray s;
    s.origin = hi.position + epsilon * w;
    s.direction = w;
    s.distance = std::sqrt(d2) - 2.f * epsilon;
    s.depth = r.depth + 1;
    if (current->intersect(s))
        return false;

    //compute emitted
    Li = exitance_ / float(M_PI) * costheta2 / (d2 * pface * parea);

    return true;
}

Vec3f mesh::exitance(void) const
{
    return exitance_;
}

void mesh::set_exitance(const CGLA::Vec3f& E)
{
    exitance_ = E;
}

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
