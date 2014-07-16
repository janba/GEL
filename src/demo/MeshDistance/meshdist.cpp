#include <iostream>

#include <GEL/Util/Timer.h>
#include <GEL/Util/ArgExtracter.h>

#include <GEL/CGLA/Mat4x4d.h>

#include <GEL/Geometry/RGrid.h>
#include <GEL/Geometry/HGrid.h>
#include <GEL/Geometry/save_raw.h>
#include <GEL/Geometry/GridAlgorithm.h>
#include <GEL/Geometry/build_bbtree.h>
#include <GEL/Geometry/AABox.h>

#include <GEL/HMesh/triangulate.h>
#include <GEL/HMesh/Manifold.h>

#include <GEL/HMesh/load.h>
#include <GEL/HMesh/x3d_load.h>
#include <GEL/HMesh/x3d_save.h>


using namespace std;
using namespace HMesh;
using namespace Geometry;
using namespace CGLA;
using namespace Util;

namespace
{
    
	Vec3i vol_dim(64);
	
    const Mat4x4d fit_bounding_volume(const Vec3d& p0,
                                      const Vec3d& p7,
                                      float buf_reg)
	{
		Vec3d sz = p7 - p0;
		Vec3i dims = vol_dim;
		Vec3d scal_vec = (Vec3d(dims)-Vec3d(2*buf_reg+2))/sz;
		float scal = min(scal_vec[0], min(scal_vec[1],scal_vec[2]));
		
		Mat4x4d m = translation_Mat4x4d(Vec3d(0)+Vec3d(buf_reg+1));
		m *= scaling_Mat4x4d(Vec3d(scal));
		m *= translation_Mat4x4d(-p0);
		return m;
	}
    
	bool do_ray_tests = false;
	bool do_obb = true;
	bool do_aabb = false;
	bool flip_normals = false;
    
}



template<class BBTree>
class DistCompCache
{
	BBTree *T;
	float old_d;
	Vec3i old_p;
public:
    
	DistCompCache(BBTree* _T): T(_T), old_p(-99999) {}
    
	void operator()(const CGLA::Vec3i& pi, float& vox_val)
	{
		Vec3f p(pi);
 		if(sqr_length(pi-old_p)==1)
            vox_val = T->compute_signed_distance(p,CGLA::sqr(1.001+fabs(old_d)));
 		else
			vox_val = T->compute_signed_distance(p);
		if(flip_normals) vox_val = -vox_val;
		old_p = pi;
		old_d = vox_val;
	}
};

template<class BBTree>
class DistComp
{
	BBTree *T;
public:
    
	DistComp(BBTree* _T): T(_T) {}
    
	void operator()(const CGLA::Vec3i& pi, float& vox_val)
	{
		Vec3d p(pi);
		vox_val =  T->compute_signed_distance(p);
		if(flip_normals) vox_val = -vox_val;
	}
};


template<class BBTree>
class RayCast
{
	BBTree *T;
    
public:
	
	RayCast(BBTree* _T): T(_T) {}
    
	void operator()(const CGLA::Vec3i& pi, float& vox_val)
	{
		int n = T->intersect_cnt(Vec3f(pi), Vec3f(1,0,0));
		if(n%2==0)
			vox_val=1000;
		else
			vox_val=-1000;
	}
};


typedef RGrid<float> RGridf;


int main(int argc, char** argv)
{
	// LOAD OBJ
    Manifold m;
    string fn = "../../data/bunny-little.x3d";
    if(argc>1)
	{
		ArgExtracter ae(argc, argv);
		
		do_aabb = ae.extract("-A");
		do_obb = ae.extract("-O");
		ae.extract("-x", vol_dim[0]);
		ae.extract("-y", vol_dim[1]);
		ae.extract("-z", vol_dim[2]);
		do_ray_tests = ae.extract("-R");
		flip_normals = ae.extract("-f");
		fn = ae.get_last_arg();
    }
    if(!load(fn, m))
    {
        cout << "Failed to load mesh" << endl;
        exit(1);
    }
    string file_prefix = fn.substr(0, fn.length()-4) + "-";
	cout << "Volume dimensions " << vol_dim << endl;
	if(!valid(m))
	{
		cout << "Not a valid manifold" << endl;
		exit(0);
	}
	triangulate_by_edge_face_split(m);
	
	Vec3d p0, p7;
	bbox(m, p0, p7);

  cout << p0 << endl;
  cout << p7 << endl;
	
	Mat4x4d T = fit_bounding_volume(p0,p7,10);
    
  cout << "Transformation " << T << endl;
	
	for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v)
		m.pos(*v) = T.mul_3D_point(m.pos(*v));
	
	
 	RGridf grid(vol_dim,FLT_MAX);
	Util::Timer tim;
	
	
	float T_build_obb=0, T_build_aabb=0, T_dist_obb=0,
    T_dist_aabb=0, T_ray_obb=0, T_ray_aabb=0;
	
	if(do_obb)
	{
		cout << "Building OBB Tree" << endl;
		tim.start();
		OBBTree obb_tree;
		build_OBBTree(m, obb_tree);
		T_build_obb = tim.get_secs();
		
		cout << "Computing distances from OBB Tree" << endl;
		tim.start();
		DistCompCache<OBBTree> dist(&obb_tree);
		for_each_voxel(grid, dist);
		T_dist_obb = tim.get_secs();
		
		cout << "Saving distance field" << endl;
		save_raw_float(file_prefix+"obb_dist.raw", grid);
		
		if(do_ray_tests)
		{
			cout << "Ray tests on OBB Tree" << endl;
			tim.start();
			RayCast<OBBTree> ray(&obb_tree);
			for_each_voxel(grid, ray);
			T_ray_obb = tim.get_secs();
			
			cout << "Saving ray volume" << endl;
			save_raw_float(file_prefix+"obb_ray.raw", grid);
		}
	}
	
	if(do_aabb)
	{
		cout << "Building AABB Tree" << endl;
		tim.start();
		AABBTree aabb_tree;
		build_AABBTree(m, aabb_tree);
		T_build_aabb = tim.get_secs();
		
		cout << "Computing distances from AABB Tree" << endl;
		tim.start();
		DistCompCache<AABBTree> dist(&aabb_tree);
		for_each_voxel(grid, dist);
		T_dist_aabb = tim.get_secs();
		
		cout << "Saving distance field" << endl;
		save_raw_float(file_prefix+"aabb_dist.raw", grid);
		
		if(do_ray_tests)
		{
			cout << "Ray tests on AABB tree" << endl;
			tim.start();
			RayCast<AABBTree> ray(&aabb_tree);
			for_each_voxel(grid, ray);
			T_ray_aabb = tim.get_secs();
			
			cout << "Saving ray volume" << endl;
			save_raw_float(file_prefix+"aabb_ray.raw", grid);
		}
	}
	cout.width(10);
	cout << "Poly";
	cout.width(11);
	cout <<"build_obb";
	cout.width(12);
	cout << "build_aabb";
	cout.width(10);
	cout << "dist_obb" ;
	cout.width(10);
	cout << "dist_aabb";
	cout.width(10);
	cout << "ray_obb" ;
	cout.width(10);
	cout << "ray_aabb";
	cout << endl;
	
	cout.precision(4);
	cout.width(10);
	cout << m.no_faces() << " ";
	cout.width(10);
	cout << T_build_obb;
	cout.width(12);
	cout << T_build_aabb;
	cout.width(10);
	cout << T_dist_obb;
	cout.width(10);
	cout << T_dist_aabb;
	cout.width(10);
	cout << T_ray_obb;
	cout.width(10);
	cout << T_ray_aabb;
	cout << endl;
  system("PAUSE");
}