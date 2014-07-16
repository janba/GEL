#include <GEL/Geometry/BSPTree.h>

#include "scene.h"
#include "mesh.h"

using namespace Geometry;
using namespace CGLA;

scene::scene(void) : cam_(0)
{
    accel_ = new BSPTree;
}

scene::~scene(void)
{
    delete accel_;
}

void scene::insert(const luminaire* obj)
{
    objs_.push_back(obj);
}

size_t scene::luminaires(void) const
{
    return objs_.size();
}

const luminaire* scene::get_luminaire(size_t i) const
{
    return objs_.at(i);
}

bool scene::intersect(const ray& r)
{
    Ray s;
    s.direction = r.direction;
    s.origin = r.origin;
    s.dist = r.distance;
    s.trace_depth = r.depth;

    return accel_->intersect(s);
}

bool scene::intersect(const ray& r, hit_info& hi)
{
    Ray s;
    s.direction = r.direction;
    s.origin = r.origin;
    s.trace_depth = r.depth;
    s.dist = r.distance;

    bool hit = accel_->intersect(s);

    if (hit)
    {
        //find the mesh which was hit
        const mesh& msh = dynamic_cast<const mesh&>(*objs_[s.id]);

        //set hit position and normals
        hi.distance = s.dist;
        hi.position = s.hit_pos;
        hi.shading_normal = s.hit_normal;

        const IndexedFaceSet& geometry = msh.get_trimesh().geometry;
        const Vec3i& f  = geometry.face(s.hit_face_id);
        const Vec3f p0 = geometry.vertex(f[0]);
        const Vec3f a  = geometry.vertex(f[1]) - p0;
        const Vec3f b  = geometry.vertex(f[2]) - p0;
        Vec3f face_normal = cross(a,b);
        float l = sqr_length(face_normal);
        if(l > 0.0f)
            face_normal /= sqrt(l);

        hi.geometric_normal = msh.transform().mul_3D_vector(face_normal);
        //hi.geometric_normal = msh.triangle_normal(s.hit_face_id);

        hi.inside = dot(hi.geometric_normal, -s.direction) < 0.f;
        if (hi.inside)
        {
            hi.geometric_normal = -hi.geometric_normal;
            hi.shading_normal = -hi.shading_normal;
        }

        bool has_texcoords = msh.get_trimesh().texcoords.no_faces() > 0;

        if (has_texcoords)
        {
            //compute texture coords..
            const IndexedFaceSet& texcoords = msh.get_trimesh().texcoords;
            const Vec3i& t = texcoords.face(s.hit_face_id);
            Vec3f uv0 = texcoords.vertex(t[0]);
            Vec3f uv1 = texcoords.vertex(t[1]);
            Vec3f uv2 = texcoords.vertex(t[2]);
            hi.texcoords(0) = (1.f-s.u-s.v)*uv0(0) + s.u*uv1(0) + s.v*uv2(0);
            hi.texcoords(1) = (1.f-s.u-s.v)*uv0(1) + s.u*uv1(1) + s.v*uv2(1);

            //compute tangent vectors..


        }
        else
        {
            //hmm.. no uvs. not good!
            hi.texcoords = Vec2f(0.f);
            //orthogonal(hi.shading_normal, hi.tangent, hi.bitangent);
        }

        //sample material to get bsdf
        const material* mat = msh.get_material();
        if (mat)
            mat->sample(r, hi);

        //convert from radiant exitance to radiance
        if (!hi.inside)
            hi.emitted = msh.exitance() / float(M_PI);
    }

    return hit;
}

const camera* scene::get_camera(void) const
{
    return cam_;
}

void scene::set_camera(const camera* c)
{
    cam_ = c;
}

void scene::initialize(int mo, int md)
{
    assert(cam_);
    assert(!objs_.empty());

    std::cout << "Building BSP tree.." << std::endl;

    std::vector<const TriMesh*> msh;
    std::vector<Mat4x4f> trs;

    std::vector<const luminaire*>::iterator it;
    for (it=objs_.begin(); it!=objs_.end(); ++it)
    {
        const luminaire* obj = *it;

        const mesh* mxx = dynamic_cast<const mesh*>(obj);
        if (!mxx)
            continue;

        msh.push_back(&mxx->get_trimesh());
        trs.push_back(mxx->transform());
    }

    accel_->init(msh, trs, mo, md);
    accel_->build();

    std::cout << "Done!" << std::endl;
}

//02566 framework, Anders Wang Kristensen, awk@imm.dtu.dk, 2007
