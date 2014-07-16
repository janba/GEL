/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <cfloat>
#include <queue>
#include <algorithm>
#include <vector>
#include "AABox.h"
#include "OBox.h"
#include "BoundingTree.h"

using namespace std;
using namespace CGLA;


namespace
{
	template <class _Tp>
	inline const _Tp& my_min(const _Tp& __a, const _Tp& __b)
	{
		return __b < __a ? __b : __a;
	}
	template <class _Tp>
	inline const _Tp& my_max(const _Tp& __a, const _Tp& __b) 
	{
		return  __a < __b ? __b : __a;
	}

}

namespace Geometry
{

template<class BoxType>
bool BoundingTree<BoxType>::intersect(const CGLA::Vec3f& p , const CGLA::Vec3f& dir,
																float& tmin) const 
{
	return root->intersect(p,dir,tmin);
}

template<class BoxType>
void BoundingTree<BoxType>::intersect(Ray& r) const 
{
		root->intersect(r);
}

template<class BoxType>
int BoundingTree<BoxType>::intersect_cnt(const CGLA::Vec3f& p, 
																	 const CGLA::Vec3f& dir) const
{
	return root->intersect_cnt(p,dir);
}

template<class BoxType>
void BoundingTree<BoxType>::build(std::vector<Triangle>& triangles)
{
	delete root;
	root = Node::build(triangles);
}


template<class Node>
class HE
{
	const Node* node;
	float sq_dist_min;
	float sq_dist_max;
	float sgn;

public:

	HE() {}

	HE(const Vec3f& p, const Node*_node): 
		node(_node), sq_dist_min(FLT_MAX), sq_dist_max(FLT_MAX), sgn(0)
																				
	{
		node->sq_distance(p,sq_dist_min,sq_dist_max, sgn);
	}
	
	float get_sq_dist_min() const {return sq_dist_min;}
	float get_sq_dist_max() const {return sq_dist_max;}

	float get_dist() const {return sgn * sqrt(sq_dist_min);}

	const Node* get_node() const 
	{
		return node;
	}

	bool operator<(const HE<Node>& r) const
	{
		return r.sq_dist_min< sq_dist_min;
	}

};



template<class BoxType>
float BoundingTree<BoxType>::compute_signed_distance(const CGLA::Vec3f& p,
																										 float minmax) const
{
	int N=100;
	vector<HE<Node> > Q(N);
	Q[0] = HE<Node>(p,root);
	
	HE<Node> *Q_beg = &Q[0];
	int Q_end = 1;
	int pushes = 1;
	while(const IntNode* n = dynamic_cast<const IntNode*>(Q[0].get_node()))
		{
			float q0_max= Q[0].get_sq_dist_max();
			//float q0_min= Q[0].get_sq_dist_min();
			pop_heap(Q_beg, Q_beg + Q_end);
			--Q_end;
			

			HE<Node> hl(p,n->get_left());
			if(hl.get_sq_dist_min() < (minmax + DIST_THRESH))
				{
					if(hl.get_sq_dist_max() < minmax)
							minmax = hl.get_sq_dist_max();
					
					Q[Q_end++] = hl;
					push_heap(Q_beg, Q_beg + Q_end);
					if(Q_end == N) 
						{
							Q.resize(N=2*N);
							Q_beg = &Q[0];
						}
					++pushes;
				}

			HE<Node> hr(p,n->get_right());
			if(hr.get_sq_dist_min() < (minmax + DIST_THRESH))
				{
					if(hr.get_sq_dist_max() < minmax)
							minmax = hr.get_sq_dist_max();

					Q[Q_end++] = hr;
					push_heap(Q_beg, Q_beg + Q_end);
					if(Q_end == N)
						{
							Q.resize(N=2*N);
							Q_beg = &Q[0];
						}
					++pushes;
				}

//  			if((hr.get_sq_dist_min() > (q0_max + DIST_THRESH)) &&
// 				 (hl.get_sq_dist_min() > (q0_max + DIST_THRESH)) )
// 				{
// 					cout.precision(4);
// 					cout << q0_min << " " << q0_max << endl;
// 					cout << hl.get_sq_dist_min() << endl;
// 					cout << hr.get_sq_dist_min() << endl;
// 					cout << typeid(*n).name() << endl;
// 					if(const LeafNode* ll =dynamic_cast<const LeafNode*>(hl.get_node()))
// 						{
// 							ll->get_tri().print();
// 							cout << sqr_length(p-ll->get_tri().get_v0()) << endl;
// 							cout << sqr_length(p-ll->get_tri().get_v1()) << endl;
// 							cout << sqr_length(p-ll->get_tri().get_v2()) << endl;
// 							float d=FLT_MAX, s;
// 							ll->get_tri().signed_distance(p,d,s);
// 							cout << "Dist " << d << endl;
// 						}
// 					if(const LeafNode* lr =dynamic_cast<const LeafNode*>(hr.get_node()))
// 						{
// 							lr->get_tri().print();
// 							cout << sqr_length(p-lr->get_tri().get_v0()) << endl;
// 							cout << sqr_length(p-lr->get_tri().get_v1()) << endl;
// 							cout << sqr_length(p-lr->get_tri().get_v2()) << endl;
// 							float d=FLT_MAX, s;
// 							lr->get_tri().signed_distance(p,d,s);
// 							cout << "Dist " << d << endl;
// 						}
// 					cout << "P=" << p<< endl;
// 				}

 			assert((hr.get_sq_dist_min() < (q0_max + DIST_THRESH)) ||
 						 (hl.get_sq_dist_min() < (q0_max + DIST_THRESH)) );
			assert(Q_end > 0);
			assert(Q_end <=N);
		}
	return Q[0].get_dist();
}

template class BoundingTree<AABox>;
template class BoundingTree<OBox>;

}
