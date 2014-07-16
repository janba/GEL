/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include "../CGLA/Vec3d.h"
#include "../Geometry/TriMesh.h"

#include "BSPTree.h"

using namespace std;
using namespace CGLA;

namespace Geometry
{
  int BSPTree::node_calls;
  int BSPTree::tri_calls;

  void create_tri_accel(Vec3f A, Vec3f B, Vec3f C, TriAccel &tri_accel) 
  {
    Vec3f N = normalize(cross(B-A, C-A));
    // Find projection dir
    int k,u,v;
    if (fabs(N[0]) > fabs(N[1]))
      if (fabs(N[0]) > fabs(N[2])) 
    k = 0; /* X */ 
      else 
    k=2;   /* Z */
    else
      if (fabs(N[1]) > fabs(N[2])) 
    k = 1; /* Y */ 
      else 
    k=2;   /* Z */
    
    u = (k+1)% 3; 
    v = (k+2)% 3;

    Vec3f b = C - A;
    Vec3f c = B - A;
    
    double div = (b[u]*c[v]-b[v]*c[u]);

    tri_accel.n_u = N[u]/N[k];
    tri_accel.n_v = N[v]/N[k];
    tri_accel.n_d = dot(A, N/N[k]);
    tri_accel.k = k;
    
    tri_accel.b_nu = -b[v]/div;
    tri_accel.b_nv = b[u]/div;
    tri_accel.b_d = (b[v]*A[u]-b[u]*A[v])/div;
    
    tri_accel.c_nu = c[v]/div;
    tri_accel.c_nv = -c[u]/div;
    tri_accel.c_d = (c[u]*A[v]-c[v]*A[u])/div;
  }
  
  // Most of this file is a direct copy from Henrik's notes
  BSPTree::BSPTree() 
  {
    root=0;
    b_is_build = false;
  }

  BSPTree::~BSPTree() 
  {
    clear();        
  }

  void BSPTree::clear() 
  {
    if (root!=0) 
    {
      delete_node(root);
      root=0;
      b_is_build=false;
    }
    isecttris.clear();
    all_objects.clear();
    all_triaccel.clear();
  }

  void BSPTree::delete_node(BSPNode *node) 
  {
    if (node->left!=0)
      delete_node(node->left);
    if (node->right!=0)
      delete_node(node->right);
    delete node;
  }

  void BSPTree::subdivide_node(BSPNode &node, BBox &bbox, 
                               unsigned int level, 
                               vector<ISectTri*>& objects, 
                               vector<TriAccel*>& tri_objects) 
  {
    const int TESTS = 2;
    
    if (objects.size()<=max_objects || level==max_level) 
    {
      node.axis_leaf = 4; // Means that this is a leaf
      node.plane = 0;
      node.left = 0;
      node.right = 0;
      node.id = all_objects.size();
      node.count = objects.size();
        
      for(unsigned int i = 0; i < objects.size(); ++i) 
      {
        all_objects.push_back(objects[i]);
        all_triaccel.push_back(tri_objects[i]);
      }  
    } 
    else 
    {
      bool right_zero=false;
      bool left_zero=false;
      unsigned int i;
      BSPNode* left_node  = new BSPNode();
      BSPNode* right_node = new BSPNode();
      vector<ISectTri*> left_objects;
      vector<ISectTri*> right_objects;
      vector<TriAccel*> tri_left_objects;
      vector<TriAccel*> tri_right_objects;
      
      node.left = left_node;
      node.right = right_node;

      int new_axis=-1;
      double min_cost=CGLA::BIG;
      int new_pos = -1;      

      for(i=0;i<3;i++) 
      {
        for(int k=1;k<TESTS;k++) 
        {
          BBox left_bbox = bbox;
          BBox right_bbox = bbox;
                    
          double center = (bbox.max_corner[i]- bbox.min_corner[i])*(double)k/(double)TESTS + bbox.min_corner[i];
          node.plane = center;
          
          left_bbox.max_corner[i] = center; 
          right_bbox.min_corner[i] = center; 

          // Try putting the triangles in the left and right boxes
          int left_count = 0;
          int right_count = 0;
          for(unsigned int j=0;j<objects.size();j++) 
          {
            ISectTri* tri = objects[j];
            left_count += left_bbox.intersect_triangle(*tri);
            right_count += right_bbox.intersect_triangle(*tri);
          }

          //double len = bbox.max_corner[i] - bbox.min_corner[i];
          double cost = left_count*left_bbox.area() + right_count*right_bbox.area(); // - len*len;
          if(cost < min_cost) 
          {
            min_cost = cost;
            new_axis = i;
            new_pos = k;
            right_zero = (right_count==0);
            left_zero = (left_count==0);
          }
        }
      }
      node.axis_leaf = new_axis;
      left_node->axis_leaf = static_cast<unsigned char>(-1); 
      right_node->axis_leaf = static_cast<unsigned char>(-1); 

      // Now chose the right splitting plane
      BBox left_bbox = bbox;
      BBox right_bbox = bbox;

      double size = bbox.max_corner[node.axis_leaf]- bbox.min_corner[node.axis_leaf];
      double center = size*(double)new_pos/(double)TESTS + bbox.min_corner[node.axis_leaf];
      double diff = f_eps < size/8.0 ? size/8.0 : f_eps;
      
      if (left_zero) 
      {
        // Find min position of all triangle vertices and place the center there
        center = bbox.max_corner[node.axis_leaf];
        for(unsigned int j=0;j<objects.size();j++) 
        {
          ISectTri* tri = objects[j];
          if (tri->point0[node.axis_leaf]<center)
            center=tri->point0[node.axis_leaf];
          if (tri->point1[node.axis_leaf]<center)
            center=tri->point1[node.axis_leaf];
          if (tri->point2[node.axis_leaf]<center)
            center=tri->point2[node.axis_leaf];
        }
        center -= diff;
      }
      if (right_zero) 
      {
        // Find max position of all triangle vertices and place the center there
        center = bbox.min_corner[node.axis_leaf];
        for(unsigned int j=0;j<objects.size();j++) 
        {
          ISectTri* tri = objects[j];
          if (tri->point0[node.axis_leaf]>center)
            center=tri->point0[node.axis_leaf];
          if (tri->point1[node.axis_leaf]>center)
            center=tri->point1[node.axis_leaf];
          if (tri->point2[node.axis_leaf]>center)
            center=tri->point2[node.axis_leaf];
        }
        center += diff;
      }

      node.plane = center;
      left_bbox.max_corner[node.axis_leaf] = center; 
      right_bbox.min_corner[node.axis_leaf] = center;  
            
      // Now put the triangles in the right and left node
      for(i=0;i<objects.size();i++) 
      {
        ISectTri* tri = objects[i];
        TriAccel *tri_accel = tri_objects[i];
        if (left_bbox.intersect_triangle(*tri)) 
        {
          left_objects.push_back(tri);
          tri_left_objects.push_back(tri_accel);
        }
        if (right_bbox.intersect_triangle(*tri)) 
        {
          right_objects.push_back(tri);
          tri_right_objects.push_back(tri_accel);
        }
      }
    //if (left_zero||right_zero)
    //  cout << left_objects.size() << "," << right_objects.size() << "," << level << endl;

      objects.clear();
      subdivide_node(*left_node , left_bbox , level+1, left_objects, tri_left_objects);
      subdivide_node(*right_node, right_bbox, level+1, right_objects, tri_right_objects);
    }
  }

  void BSPTree::init() 
  {
    root = new BSPNode();
    bbox.compute_bbox(isecttris);
    bbox.min_corner-=Vec3f(1.0);
    bbox.max_corner+=Vec3f(1.0);
  }

  void BSPTree::init(vector<const TriMesh*>& _trimesh, 
                     vector<Mat4x4f>& _transforms, 
                     int _max_objects, int _max_level) 
  {
    trimesh = _trimesh;
    transforms = _transforms;
    for(unsigned int i=0;i<trimesh.size();i++) 
    {
      const TriMesh *mesh = trimesh[i];
      // Loop through all triangles and add them to intersection structure
      for(int j=0;j<mesh->geometry.no_faces();j++) 
      {
        Vec3i face = mesh->geometry.face(j);
        ISectTri new_tri;
        new_tri.point0 = transforms[i].mul_3D_point(mesh->geometry.vertex(face[0]));
        new_tri.point1 = transforms[i].mul_3D_point(mesh->geometry.vertex(face[1]));
        new_tri.point2 = transforms[i].mul_3D_point(mesh->geometry.vertex(face[2]));
        new_tri.edge0 = new_tri.point1 - new_tri.point0;
        new_tri.edge1 = new_tri.point2 - new_tri.point0;
        new_tri.mesh_id = i;
        new_tri.tri_id = j;
        isecttris.push_back(new_tri);
        TriAccel ta;
        create_tri_accel(new_tri.point0, new_tri.point1, new_tri.point2, ta);
        ta.mesh_id = i;
        ta.tri_id = j;
        triaccel.push_back(ta);
      }
    }

    max_objects = _max_objects;
    max_level = _max_level;
    init();
  }

  void BSPTree::init(const TriMesh* mesh, Mat4x4f transform, 
             vector<int> &trilist, 
             int _max_objects, int _max_level) 
  {
    trimesh.push_back(mesh);
    transforms.push_back(transform);
    // Loop through all triangles and add them to intersection structure
    for(unsigned int j=0;j<trilist.size();j++) 
    {
      Vec3i face = mesh->geometry.face(trilist[j]);
      ISectTri new_tri;
      new_tri.point0 = transform.mul_3D_point(mesh->geometry.vertex(face[0]));
      new_tri.point1 = transform.mul_3D_point(mesh->geometry.vertex(face[1]));
      new_tri.point2 = transform.mul_3D_point(mesh->geometry.vertex(face[2]));
      new_tri.edge0 = new_tri.point1 - new_tri.point0;
      new_tri.edge1 = new_tri.point2 - new_tri.point0;
      new_tri.mesh_id = 0;
      new_tri.tri_id = trilist[j];
      isecttris.push_back(new_tri);
      TriAccel ta;
      create_tri_accel(new_tri.point0, new_tri.point1, new_tri.point2, ta);
      ta.mesh_id = 0;
      ta.tri_id = trilist[j];
      triaccel.push_back(ta);
    }

    max_objects = _max_objects;
    max_level = _max_level;
    init();
  }

  void BSPTree::build() 
  {
    if (!b_is_build) 
    {
      vector<ISectTri*> objects;
      vector<TriAccel*> tri_objects;
      for(unsigned int i=0;i<isecttris.size();i++) 
      {
        ISectTri* tri = &isecttris[i];
        TriAccel* tri_accel = &triaccel[i];
        objects.push_back(tri);
        tri_objects.push_back(tri_accel);
      }
      subdivide_node(*root, bbox, 0, objects, tri_objects);
      b_is_build = true;
    }
    make_fast_tree(root);
  }

  bool BSPTree::is_build() 
  {
    return b_is_build;
  }

  void BSPTree::print(BSPNode *node, int depth) 
  {
    if (node==0)
      return;
    for(int i=0;i<depth;i++)
      cout << " ";
//  cout << "axis:" << node->axis_leaf << ", count:" << node->objects.size() << ", plane:" << node->plane << ", " << endl;
    print(node->left, depth+1);
    print(node->right, depth+1);
  }

  int BSPTree::size(BSPNode *node) 
  {
    if (node==0)
      return 0;
    int s = sizeof(BSPNode);
    s+= node->count * sizeof(ISectTri);
    s+=size(node->left);
    s+=size(node->right);
    return s;
  }

  int BSPTree::size() 
  {
    return size(root);
  }

/*__declspec(align(16))*/ static const unsigned int modulo[] = {0,1,2,0,1};

  inline bool intersect2(Ray &ray, const TriAccel &acc, double t_max) 
  {
//inline bool Intersect(TriAccel &acc,Ray &ray)
#define ku modulo[acc.k+1]
#define kv modulo[acc.k+2]
    // don’t prefetch here, assume data has already been prefetched
    // start high-latency division as early as possible
    const double nd = 1.0/((double)ray.direction[acc.k] + (double)acc.n_u * (double)ray.direction[ku] + (double)acc.n_v * (double)ray.direction[kv]);
    const double f = ((double)acc.n_d - (double)ray.origin[acc.k]   - (double)acc.n_u * (double)ray.origin[ku] - (double)acc.n_v * (double)ray.origin[kv]) * nd;
    // check for valid distance.
    if (!(t_max > f && f > 0.001)||ray.dist<f) return false;
    // compute hitpoint positions on uv plane
    const double hu = (ray.origin[ku] + f * ray.direction[ku]);
    const double hv = (ray.origin[kv] + f * ray.direction[kv]);
    // check first barycentric coordinate
    const double lambda = (hu * (double)acc.b_nu + hv * (double)acc.b_nv + (double)acc.b_d);
    if (lambda < 0.0) return false;
    // check second barycentric coordinate
    const double mue = (hu * (double)acc.c_nu + hv * (double)acc.c_nv + (double)acc.c_d);
    if (mue < 0.0) return false;
    // check third barycentric coordinate
    if (lambda+mue > 1.0) return false;
    // have a valid hitpoint here. store it.
    ray.dist = f;
    ray.u = lambda;
    ray.v = mue;
    ray.hit_object = (TriMesh*)acc.mesh_id;
    ray.hit_face_id = acc.tri_id;
    ray.has_hit=true;
    return true;
  }

  bool BSPTree::intersect_node(Ray &ray, const BSPNode &node, double t_min, double t_max) const 
  {
    node_calls++;    
    if (node.axis_leaf==4) 
    {
      bool found = false; 
      for(int i=0; i < node.count; ++i) 
      {
        const ISectTri* tri = all_objects[node.id+i];
        if (intersect(ray, *tri, t_max))  
          found=true;
        //const TriAccel* tri2 = all_triaccel[node.id+i];
        //if (intersect2(ray, *tri2, t_max))  
        //  found=true;
      }
      if (found)
        return true;
      else 
        return false;
    } 
    else 
    {
      BSPNode *near_node;
      BSPNode *far_node;
      if (ray.direction[node.axis_leaf]>=0) 
      {
        near_node = node.left;
        far_node = node.right;
      } 
      else 
      {
        near_node = node.right;
        far_node = node.left;
      }

      // In order to avoid instability
      double t;
      if (fabs(ray.direction[node.axis_leaf])<d_eps)
        t = (node.plane - ray.origin[node.axis_leaf])/d_eps;// intersect node plane;
      else
        t = (node.plane - ray.origin[node.axis_leaf])/ray.direction[node.axis_leaf];// intersect node plane;
      
      if (t>t_max) 
        return intersect_node(ray, *near_node, t_min, t_max);      
      else if (t<t_min) 
        return intersect_node(ray, *far_node, t_min, t_max);
      else 
      {
        if (intersect_node(ray, *near_node, t_min, t))
          return true;
        else 
          return intersect_node(ray, *far_node, t, t_max);
      }
    }
  }

  bool BSPTree::intersect(Ray &ray) const 
  {
    double t_min, t_max;
    bbox.intersect_min_max(ray, t_min, t_max);
    if (t_min>t_max)
      return false;

    if (!intersect_node(ray, *root, t_min, t_max))
      return false;
    //intersect_fast_node(ray, &fast_tree[0], t_min, t_max);
    //if (!ray.has_hit)
    //  return false;
    else 
    {
      // Calculate the normal at the intersection
      ray.id = reinterpret_cast<size_t>(ray.hit_object);
      ray.hit_object = trimesh[ray.id];
      
      const Vec3i& face = ray.hit_object->normals.face(ray.hit_face_id);
      const Vec3f& normal0 = ray.hit_object->normals.vertex(face[0]);
      const Vec3f& normal1 = ray.hit_object->normals.vertex(face[1]);
      const Vec3f& normal2 = ray.hit_object->normals.vertex(face[2]);
      ray.hit_normal = transforms[ray.id].mul_3D_vector(
        normalize(normal0*(1 - ray.u - ray.v) + normal1*ray.u + normal2*ray.v));
      ray.hit_pos = ray.origin + ray.direction*ray.dist;
/*
      const Vec3i& face = ray.hit_object->normals.face(ray.hit_face_id);
      const Vec3f& normal0 = ray.hit_object->normals.vertex(face[0]);
      const Vec3f& normal1 = ray.hit_object->normals.vertex(face[1]);
      const Vec3f& normal2 = ray.hit_object->normals.vertex(face[2]);
      ray.hit_normal = normalize(normal0*(1 - ray.u - ray.v) + normal1*ray.u + normal2*ray.v);
      ray.hit_pos = ray.origin + ray.direction*ray.dist;
*/
      return true;
    }
  }

  const int MAX_DEPTH=25;

  void BSPTree::make_fast_tree(BSPNode *node) 
  {
    fast_tree.resize(1000000); // 
    FastBSPNode f_node;
    fast_tree.push_back(f_node);
    push_fast_bsp_node(node,0);
  }

  void BSPTree::push_fast_bsp_node(BSPNode *node, int id) 
  {
    if (node->axis_leaf==4)  // It is a leaf
    {
      //assert(false);
      //TODO: cant compile on 64 bit gcc

      //fast_tree[id].leaf.flagAndOffset = (unsigned int)1<<31 | (unsigned int)(&all_triaccel[node->id]);
      fast_tree[id].leaf.count = node->count;
    } 
    else // It is an inner node
    { 
      FastBSPNode fnode;
      int p_l = fast_tree.size();
      fast_tree.push_back(fnode); // left
      fast_tree.push_back(fnode); // right
      push_fast_bsp_node(node->left, p_l);
      push_fast_bsp_node(node->right, p_l+1);

      //assert(false);
      //TODO: gcc64 bit failure

      //fast_tree[id].inner.flagAndOffset = (unsigned int) &fast_tree[p_l] | node->axis_leaf;
      fast_tree[id].inner.splitCoordinate = node->plane;
      node->ref = fast_tree[id].inner.flagAndOffset;
    }
  }

#define ABSP_ISLEAF(n)       (n->inner.flagAndOffset & (unsigned int)1<<31)
#define ABSP_DIMENSION(n)    (n->inner.flagAndOffset & 0x3)
#define ABSP_OFFSET(n)       (n->inner.flagAndOffset & (0x7FFFFFFC))
#define ABSP_NEARNODE(n)     (FastBSPNode*)(ray.direction[dimension]>=0?ABSP_OFFSET(node):ABSP_OFFSET(node)+sizeof(*node))
#define ABSP_FARNODE(n)      (FastBSPNode*)(ray.direction[dimension]>=0?ABSP_OFFSET(node)+sizeof(*node):ABSP_OFFSET(node))
  
  struct Stack 
  {
    FastBSPNode *node;
    double t_min;
    double t_max;
  };

  inline void IntersectAlltrianglesInLeaf(const BSPLeaf* leaf, Ray &ray, double t_max) {
    TriAccel** tri_acc_ptr = reinterpret_cast<TriAccel**>(leaf->flagAndOffset & (0x7FFFFFFF));
    for(unsigned int i = 0; i < leaf->count; ++i)
      intersect2(ray, *(*tri_acc_ptr + i), t_max);
  }

  void BSPTree::intersect_fast_node(Ray &ray, const FastBSPNode *node, double t_min, double t_max) const 
  {
    Stack stack[MAX_DEPTH];
    int stack_id=0;
    double t;
    // Precalculate one over dir
    double one_over_dir[3];
    for(int i=0;i<3;i++) 
    {
      if (ray.direction[i]!=0)
        one_over_dir[i]=1.0/ray.direction[i];
      else
        one_over_dir[i]=1.0/d_eps;
    }

    int dimension;
    while(1) 
    {
      while(!ABSP_ISLEAF(node)) 
      {
        dimension = ABSP_DIMENSION(node);
        t = (node->inner.splitCoordinate - ray.origin[dimension])*one_over_dir[dimension];
        if (t>=t_max) 
          node = ABSP_NEARNODE(node);
        else if (t<=t_min)
          node = ABSP_FARNODE(node);
        else 
        {
          // Stack push
          stack[stack_id].node = ABSP_FARNODE(node);
          stack[stack_id].t_min = t;
          stack[stack_id++].t_max = t_max;
          // Set current node to near side
          node = ABSP_NEARNODE(node);
          t_max = t;
        }
      }
      
      IntersectAlltrianglesInLeaf(&node->leaf, ray, t_max);
      if (ray.dist<t_max)
        return;
      if (stack_id==0)
        return;
      // Stack pop
      
      node = stack[--stack_id].node;
      t_min = stack[stack_id].t_min;
      t_max = stack[stack_id].t_max;
    }
  }

  bool BSPTree::intersect(Ray &ray, const ISectTri &isecttri, double t_max) const 
  {
    tri_calls++;

    // This is the Möller-Trumbore method
    Vec3d direction(ray.direction);
    Vec3d edge0(isecttri.edge0);
    Vec3d edge1(isecttri.edge1);

    // Ray-triangle intersection
    Vec3d p = cross(direction, edge1);
    double a = dot(edge0, p);
    if(a > -d_eps && a < d_eps)
      return false;

    // Just delay these 
    Vec3d origin(ray.origin);
    Vec3d point0(isecttri.point0);    
    double f = 1.0/a;
    Vec3d s = origin - point0;
    double u = f*dot(s, p);
    if(u < 0.0 || u > 1.0)
      return false;

    Vec3d q = cross(s, edge0);
    double v = f*dot(direction, q);  
    if(v < 0.0 || u + v > 1.0)
      return false;

    double t = f*dot(edge1, q);
    if(t < f_eps || t*t < 1.0e-9)
      return false;
    if(t > t_max)
      return false;
    if(t > ray.dist)
      return false;
  
    ray.dist = t;
    ray.u = u;
    ray.v = v;
    ray.hit_object = (TriMesh*)isecttri.mesh_id;
    ray.hit_face_id = isecttri.tri_id;
    ray.has_hit=true;
    return true; 
  }
}
