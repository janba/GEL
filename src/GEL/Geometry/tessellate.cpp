/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <map>
#include <algorithm>
#include <queue>

#include "../Util/Grid2D.h"
#include "../CGLA/Vec3f.h"
#include "../CGLA/Vec2f.h"

#include "tessellate.h"

using namespace std;
using namespace CGLA;
using namespace Util;

namespace Geometry
{
	
	float MAX_ERR = 5;
	float MAX_DIST = 5;
	int  ADAPTIVE = false;
	
#define mymax(x,y) ((x)>(y)?(x):(y))
	
	// Test surfaces ---
	
	
	// ----------- CONSTANTS AND GLOBALS ------------
	
	namespace
	{
		float maximum_edge_error=0;
		
		IndexedFaceSet* face_set_ptr = 0;
		ParSurf* the_surf = 0;
		
		class TreeNode;
		typedef TreeNode* TreeNodePtr;
		
		vector<TreeNode*> root_nodes;
		
		class Point 
		{
			friend class Edge;
			friend float subdiv(int level, const Point&, const Point&, TreeNodePtr&);
			friend float compute_mid_pt_err(const Point&, const Point&, Point&);
			
			
			Vec2f ppos;
			Vec3f pos;
			mutable int vertex_idx;
			
public:
				
				Point(const Vec2f& _ppos, const Vec3f& _pos): 
				ppos(_ppos), pos(_pos), vertex_idx(-1)
			{
			}
			
			Point() {}
			
			void create_vertex() const
			{
				if(vertex_idx == -1) 
					vertex_idx = face_set_ptr->add_vertex(pos);
			}
		};
		
		
		
		inline const Point create_point(float u, float v)
		{
			Vec2f ppos(u,v);
			Vec3f pos = (*the_surf)(u,v);
			return Point(ppos, pos);
		}
		
		float subdiv(int level, const Point& lp, const Point& rp, TreeNodePtr& node);
		
		class TreeNode
		{
			friend class Edge;
			friend float subdiv(int level, const Point&, const Point&, TreeNodePtr&);
			
			const Point p;
			const float error;
			
			const TreeNodePtr left;
			const TreeNodePtr right;
			
public:
				
				TreeNode(const Point& _p, float _error,
						 TreeNodePtr _left, TreeNodePtr _right):
				p(_p), error(_error), left(_left), right(_right)  {}
			
			~TreeNode() 
			{
				if(left)
					delete left; 
				
				if(right)
					delete right;
			}
			
		};
		
		class Edge
		{
			friend Edge create_edge(const Point&, const Point&, bool);
			
			Edge(const Point& _left_point, 
				 const Point& _right_point, 
				 TreeNodePtr _root,
				 bool _cw): 
				left_point(_left_point), right_point(_right_point), root(_root), cw(_cw)
			{
			}
			
			
			Point left_point;
			Point right_point;
			TreeNodePtr root;
			bool cw;
			
public:
				
				Edge() {}
			
			float get_error() const
			{
				if(root)
					return root->error;
				return 0.0f;
			}
			
			float get_length() const
			{
				return (left_point.pos-right_point.pos).length();
			}
			
			bool is_simple() const 
			{
				if(root->left==0)
				{
					assert(root->left==0);
					assert(root->right==0);
					return true;
				}
				return false;
			}
			
			const Edge get_sub_edge(int edge_no) const
			{
				assert(edge_no == 0 || edge_no == 1);
				assert(!is_simple());
				
				if(!cw) edge_no = 1-edge_no;
				
				if(edge_no ==0)
					return Edge(left_point, root->p, root->left, cw);
				else
					return Edge(root->p, right_point, root->right, cw);
			}
			
			const Edge get_opp_edge() const
			{
				return Edge(left_point, right_point, root, !cw);
			}
			
			
			const Point& get_split_point() const
			{
				assert(!is_simple());
				return root->p;
			}
			
			const Point& get_point(int i) const
			{
				assert(i==0 || i==1);
				if(!cw) i=1-i;
				if(i == 0) 
					return left_point;
				return right_point;
			}
			
			int get_point_idx(int i) const
			{
				assert(i==0 || i==1);
				if(!cw) i=1-i;
				if(i == 0) 
					return left_point.vertex_idx;
				return right_point.vertex_idx;
			}
			
		};
		
		
		float compute_mid_pt_err(const Point&, const Point&, Point&);
		
		Edge create_edge(const Point& left_point, const Point& right_point, 
						 bool cw=true)
		{
			TreeNodePtr root=0;
			if(ADAPTIVE)
				subdiv(0, left_point,right_point,root);
			else
			{
				Point mid_point;
				float err = compute_mid_pt_err(left_point, right_point, mid_point);
				root = new TreeNode(mid_point,err,0,0);
			}
			if(root)
				root_nodes.push_back(root);
			
			return Edge(left_point, right_point, root, cw);
		}
		
		class Triangle
		{
			Edge a,b,c;
public:
			Triangle(const Edge& _a, const Edge& _b, const Edge& _c):
			a(_a),  b(_b), c(_c) {}
			
			const Edge& operator[](int i) const
		{
			switch(i)
			{
				case 0: return a;
				case 1: return b;
				case 2: return c;
			}
			return c;
		}
			
		};
		
		float compute_mid_pt_err(const Point& lp, const Point& rp, 
								 Point& mid_point)
		{
			Vec3f m = (lp.pos - rp.pos)/2.0f + rp.pos;
			Vec2f uvm = (lp.ppos -rp.ppos)/2.0f + rp.ppos;
			
			mid_point = create_point(uvm[0], uvm[1]);
			Vec3f diff = mid_point.pos - m;
			return diff.length();
		}
		
		float subdiv(int level, const Point& lp, const Point& rp, TreeNodePtr& node)
		{
			Point mid_point;
			float err = compute_mid_pt_err(lp,rp,mid_point);
			
//			Vec3f dist = lp.pos-rp.pos;
			if (level>10) 
			{
				node = new TreeNode(mid_point,err,0,0);
				return err;
			}
			
			TreeNodePtr lnode;
			float lerr = subdiv(level+1,lp, mid_point, lnode);
			
			TreeNodePtr rnode;
			float rerr = subdiv(level+1,mid_point, rp, rnode);
			
			float max_err = mymax(err, mymax(lerr,rerr));
			
			if(max_err < MAX_ERR)
			{
				if(lnode) delete lnode;
				if(rnode) delete rnode;
				
				node = new TreeNode(mid_point, max_err, 0,0);
			}
			else
			{
				mid_point.create_vertex();
				node = new TreeNode(mid_point, max_err, lnode,rnode);
			}
			
			return max_err;
		}
		
		
		
		void split(const Triangle orig_tri)
		{
			queue<Triangle> tri_queue;
			tri_queue.push(orig_tri);
			
			
			while(!tri_queue.empty())
			{
				const Triangle tri = tri_queue.front();
				tri_queue.pop();
				
				int no_split_edges = 0;
				
				for(int i=0;i<3;++i)
					if(!tri[i].is_simple())
						++no_split_edges;
				
				if(no_split_edges==0 /* || level > 2000*/) 
				{
					int a = tri[0].get_point_idx(0);
					int b = tri[1].get_point_idx(0);
					int c = tri[2].get_point_idx(0);
					face_set_ptr->add_face(Vec3i(a,b,c));
					maximum_edge_error = mymax(mymax(tri[0].get_error(),
													 tri[1].get_error()),
											   tri[2].get_error());
				}
				else if(no_split_edges==3)
				{
					const Point& p0 = tri[0].get_point(0);
					const Point& p1 = tri[1].get_point(0);
					const Point& p2 = tri[2].get_point(0);
					
					const Point sp[3] = {tri[0].get_split_point(),
						tri[1].get_split_point(),
						tri[2].get_split_point()};
					
					Edge sp_edg[3] = {create_edge(sp[0],p2),
						create_edge(sp[1],p0),
						create_edge(sp[2],p1)};
					
					int i_min=0;
					float l_min = sp_edg[0].get_length();
					int i;
					for(i=1;i<3;++i)
					{
						float len = sp_edg[i].get_length();
						if(len<l_min) 
						{
							l_min = len;
							i_min = i;
						}
					}
					
					i = i_min;
					int j = (i+1)%3;
					int k = (i+2)%3;
					
					Edge spi_spj = create_edge(sp[i],sp[j]);
					Edge spi_spk = create_edge(sp[i],sp[k]);
					
					Triangle tri0(tri[i].get_sub_edge(1),
								  tri[j].get_sub_edge(0),
								  spi_spj.get_opp_edge());
					Triangle tri1(spi_spj,
								  tri[j].get_sub_edge(1),
								  sp_edg[i].get_opp_edge());
					Triangle tri2(sp_edg[i],
								  tri[k].get_sub_edge(0),
								  spi_spk.get_opp_edge());
					Triangle tri3(spi_spk,
								  tri[k].get_sub_edge(1),
								  tri[i].get_sub_edge(0));
					
					tri_queue.push(tri0);
					tri_queue.push(tri1);
					tri_queue.push(tri2);
					tri_queue.push(tri3);
				}
				else if(no_split_edges==2)
				{
					int i;
					if(tri[2].is_simple()) i = 0;
					else if(tri[1].is_simple()) i=2;
					else i = 1;
					
					int j = (i+1)%3;
					int k = (i+2)%3;
					
					const Point& spi = tri[i].get_split_point();
					const Point& spj = tri[j].get_split_point();
					
					const Point& pi = tri[i].get_point(0);
					const Point& pk = tri[k].get_point(0);
					
					Edge spi_spj = create_edge(spi, spj);
					Edge pi_spj = create_edge(pi, spj);
					Edge pk_spi = create_edge(pk, spi);
					
					Triangle tri0(tri[i].get_sub_edge(1),
								  tri[j].get_sub_edge(0),
								  spi_spj.get_opp_edge());
					
					
					if(pi_spj.get_length() < pk_spi.get_length())
					{
						Triangle tri1(spi_spj,
									  pi_spj.get_opp_edge(),
									  tri[i].get_sub_edge(0));
						Triangle tri2(pi_spj,
									  tri[j].get_sub_edge(1),
									  tri[k]);
						tri_queue.push(tri1);					
						tri_queue.push(tri2);	
					}
					else
					{
						Triangle tri1(tri[i].get_sub_edge(0),
									  pk_spi.get_opp_edge(),
									  tri[k]);
						Triangle tri2(spi_spj,
									  tri[j].get_sub_edge(1),
									  pk_spi);
						tri_queue.push(tri1);					
						tri_queue.push(tri2);
					}
					tri_queue.push(tri0);
				}
				else if(no_split_edges==1)
				{
					int i;
					if(!tri[0].is_simple()) i=0;
					else if(!tri[1].is_simple()) i=1;
					else i= 2;
					
					int j = (i+1)%3;
					int k = (i+2)%3;
					
					const Point& pk = tri[k].get_point(0);
					const Point& spi = tri[i].get_split_point();
					
					Edge spi_pk = create_edge(spi, pk);
					
					Triangle tri0(tri[i].get_sub_edge(1),
								  tri[j],
								  spi_pk.get_opp_edge());
					
					Triangle tri1(tri[i].get_sub_edge(0),
								  spi_pk,
								  tri[k]);
					
					tri_queue.push(tri0);
					tri_queue.push(tri1);
					
				}
			}
		}
	}
	
	
	void tessellate(IndexedFaceSet& face_set, ParSurf& s, Grid2D<Vec3f>& inigrid)
	{	
		face_set_ptr = &face_set;
		int i;
		int n=inigrid.get_xdim()-1;
		int m=inigrid.get_ydim()-1;
		the_surf = &s;
		Grid2D<Point> points(n+1,m+1);
		for(i=0;i<=n;++i)
			for(int j=0;j<=m;++j)
			{
				Vec3f p = inigrid(i,j);
				points(i,j) = create_point(p[0], p[1]);
				points(i,j).create_vertex();
			}
				
				
				Grid2D<Edge> h_edges(n,m+1);
		Grid2D<Edge> v_edges(n+1,m);
		for(i=0;i<=n;++i)
			for(int j=0;j<=m;++j)
			{
				if(i<n)
					h_edges(i,j) = create_edge(points(i,j), points(i+1,j));
				if(j<m)
					v_edges(i,j) = create_edge(points(i,j), points(i,j+1));
			}
				
				for(i=0;i<n;++i)
					for(int j=0;j<m;++j)
					{
						Edge diag = create_edge(points(i,j), points(i+1,j+1));
						Triangle tri0(h_edges(i,j), v_edges(i+1,j), diag.get_opp_edge());
						Triangle tri1(diag, h_edges(i,j+1).get_opp_edge(), 
									  v_edges(i,j).get_opp_edge());
						split(tri0);
						split(tri1);
					}
						for(i=0;i<static_cast<int>(root_nodes.size());++i)
						{
							if(root_nodes[i])
							{
								delete root_nodes[i];
							}
						}
						root_nodes.clear();
		maximum_edge_error=0;
		the_surf = 0;
	}
	
	
	void tessellate(IndexedFaceSet& face_set, ParSurf& s,
					float u_min, float u_max, float v_min, float v_max, 
					int n, int m)
	{	
		face_set_ptr = &face_set;
		int i;
		the_surf = &s;
		Grid2D<Point> points(n+1,m+1);
		float step_x = (u_max-u_min)/float(n);
		float step_y = (v_max-v_min)/float(m);
		for(i=0;i<=n;++i)
			for(int j=0;j<=m;++j)
			{
				points(i,j) = create_point(u_min+i*step_x,v_min+j*step_y);
				points(i,j).create_vertex();
			}
				
				
				Grid2D<Edge> h_edges(n,m+1);
		Grid2D<Edge> v_edges(n+1,m);
		for(i=0;i<=n;++i)
			for(int j=0;j<=m;++j)
			{
				if(i<n)
					h_edges(i,j) = create_edge(points(i,j), points(i+1,j));
				if(j<m)
					v_edges(i,j) = create_edge(points(i,j), points(i,j+1));
			}
				
				for(i=0;i<n;++i)
					for(int j=0;j<m;++j)
					{
						Edge diag = create_edge(points(i,j), points(i+1,j+1));
						Triangle tri0(h_edges(i,j), v_edges(i+1,j), diag.get_opp_edge());
						Triangle tri1(diag, h_edges(i,j+1).get_opp_edge(), 
									  v_edges(i,j).get_opp_edge());
						split(tri0);
						split(tri1);
					}
						for(i=0;i<static_cast<int>(root_nodes.size());++i)
						{
							if(root_nodes[i])
							{
								delete root_nodes[i];
							}
						}
						root_nodes.clear();
		maximum_edge_error=0;
		the_surf = 0;
	}
	
	
	/** The Edge2VertexMap is a simple database class.
		The database stores integer values referenced by keys formed by
		a pair of indices. The intended use is to let the key be a pair
		of vertices sharing an edge and the value be a new vertex inserted on
		that edge. */
	class EdgeMap
	{
		typedef pair<int,int>      VertPair;
		typedef map<VertPair, Edge> E2VMap;
		typedef E2VMap::iterator   E2VIter;
		
		E2VMap edge_to_vertex;
		
public:
			
			void set_edge(int ind0, int ind1, const Edge& new_edge)
		{
				edge_to_vertex[VertPair(ind0,ind1)] = new_edge;
		}
		
		bool find_edge(int ind0, int ind1, Edge& edg)
		{
			E2VIter iter = edge_to_vertex.find(VertPair(ind0,ind1));
			
			if(iter == edge_to_vertex.end())
				return false;
			
			edg = iter->second;
			return true;
		}
		
	};
	
	void tessellate(IndexedFaceSet& face_set, ParSurf& s, 
					std::vector<Vec2f> uv_points,
					std::vector<Vec3i> triangles)
	{		
		face_set_ptr = &face_set;
		size_t i;
		the_surf = &s;
		size_t N = uv_points.size();
		vector<Point> points(N);
		for(i=0;i<N;++i)
		{
			points[i] = create_point(uv_points[i][0], uv_points[i][1]);
			points[i].create_vertex();
		}
		EdgeMap edges;
		vector<Triangle> tris;
		for(size_t j=0;j<triangles.size();++j)
		{
			Edge e[3];
			int v[3];
			
			v[0] = triangles[j][0];
			v[1] = triangles[j][1];
			v[2] = triangles[j][2];
			
			for(int i=0;i<3;++i)
			{
				int k = (i+1)%3;
				if(!edges.find_edge(v[i],v[k],e[i]))
				{
					e[i] = create_edge(points[v[i]], points[v[k]]);
					edges.set_edge(v[k],v[i],e[i].get_opp_edge());
				}
			}
			split(Triangle(e[0],e[1],e[2]));
			
		}
		for(i=0;i<static_cast<int>(root_nodes.size());++i)
		{
			if(root_nodes[i])
			{
				delete root_nodes[i];
			}
		}
		root_nodes.clear();
		maximum_edge_error=0;
		the_surf = 0;
	}
	
}
