/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file KDTree.h
 * @brief KD Tree implementation based on a binary heap.
 */

#ifndef __GEOMETRY_KDTREE_H
#define __GEOMETRY_KDTREE_H

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "../CGLA/CGLA.h"
#include "../CGLA/ArithVec.h"

#if (_MSC_VER >= 1200)
#pragma warning (push)
#pragma warning (disable: 4018)
#endif

namespace Geometry
{
	/** \brief A classic K-D tree. 

			A K-D tree is a good data structure for storing points in space
			and for nearest neighbour queries. It is basically a generalized 
			binary tree in K dimensions. */
	template<class KeyT, class ValT>
	class KDTree
	{
		typedef typename KeyT::ScalarType ScalarType;
		typedef KeyT KeyType;
		typedef std::vector<KeyT> KeyVectorType;
		typedef std::vector<ValT> ValVectorType;
	
		/// KDNode struct represents node in KD tree
		struct KDNode
		{
			KeyT key;
			ValT val;
			short dsc;

			KDNode(): dsc(0) {}

			KDNode(const KeyT& _key, const ValT& _val):
				key(_key), val(_val), dsc(-1) {}

			ScalarType dist(const KeyType& p) const 
			{
				KeyType dist_vec = p;
				dist_vec  -= key;
				return dot(dist_vec, dist_vec);
			}
		};

		typedef std::vector<KDNode> NodeVecType;
		bool is_built;
		NodeVecType init_nodes;
		NodeVecType nodes;

		/** Comp is a class used for comparing two keys. Comp is constructed
				with the discriminator - i.e. the coordinate of the key that is used
				for comparing keys - Comp objects are passed to the sort algorithm.*/
		class Comp
		{
			const int dsc;
		public:
			Comp(int _dsc): dsc(_dsc) {}
			bool operator()(const KeyType& k0, const KeyType& k1) const
			{
				int dim=KeyType::get_dim();
				for(int i=0;i<dim;i++)
					{
						int j=(dsc+i)%dim;
						if(k0[j]<k1[j])
							return true;
						if(k0[j]>k1[j])
							return false;
					}
				return false;
			}

			bool operator()(const KDNode& k0, const KDNode& k1) const
			{
				return (*this)(k0.key,k1.key);
			}
		};


		/** Passed a vector of keys, this function will construct an optimal tree.
				It is called recursively */
		void optimize(int, int, int);

		/** Finde nearest neighbour. */
		int closest_point_priv(int, const KeyType&, ScalarType&) const;

							
		void in_sphere_priv(int n, 
												const KeyType& p, 
												const ScalarType& dist,
												std::vector<KeyT>& keys,
												std::vector<ValT>& vals) const;
	
		/** Finds the optimal discriminator. There are more ways, but this 
				function traverses the vector and finds out what dimension has
				the greatest difference between min and max element. That dimension
				is used for discriminator */
		int opt_disc(int,int) const;
	
	public:

		/** Build tree from vector of keys passed as argument. */
		KDTree(): is_built(false), init_nodes(1) {}

		/** Insert a key value pair into the tree. Note that the tree needs to 
				be built - by calling the build function - before you can search. */
		void insert(const KeyT& key, const ValT& val)
		{
				if(is_built)
				{
						assert(init_nodes.size()==1);
						init_nodes.swap(nodes);
						is_built=false;
				}
				init_nodes.push_back(KDNode(key,val));
		}

		/** Build the tree. After this function has been called, it is no longer 
				legal to insert elements, but you can perform searches. */
		void build()
		{
			assert(!is_built);
			nodes.resize(init_nodes.size());
			if(init_nodes.size() > 1)	
				optimize(1,1,init_nodes.size());
			NodeVecType v(1);
			init_nodes.swap(v);
			is_built = true;
		}

		/** Find the key value pair closest to the key given as first 
				argument. The second argument is the maximum search distance. Upon
				return this value is changed to the distance to the found point.
				The final two arguments contain the closest key and its 
				associated value upon return. */
		bool closest_point(const KeyT& p, ScalarType& dist, KeyT&k, ValT&v) const
		{
			assert(is_built);
			if(nodes.size()>1)
			{
					ScalarType max_sq_dist = CGLA::sqr(dist);
					if(int n = closest_point_priv(1, p, max_sq_dist))
					{
							k = nodes[n].key;
							v = nodes[n].val;
							dist = std::sqrt(max_sq_dist);
							return true;
					}
			}
			return false;
		}

		/** Find all the elements within a given radius (second argument) of
				the key (first argument). The key value pairs inside the sphere are
				returned in a pair of vectors passed as the two last arguments.
				Note that we don't resize the two last arguments to zero - so either
				they should be empty vectors or you should desire appending the newly
				found elements onto these vectors.				
		*/
		int in_sphere(const KeyType& p, 
									ScalarType dist,
									std::vector<KeyT>& keys,
									std::vector<ValT>& vals) const
		{
			assert(is_built);
			if(nodes.size()>1)
			{
					ScalarType max_sq_dist = CGLA::sqr(dist);
					in_sphere_priv(1,p,max_sq_dist,keys,vals);
					return keys.size();
			}
			return 0;
		}
		

	};

	template<class KeyT, class ValT>
	int KDTree<KeyT,ValT>::opt_disc(int kvec_beg,  
																	int kvec_end) const 
	{
		KeyType vmin = init_nodes[kvec_beg].key;
		KeyType vmax = init_nodes[kvec_beg].key;
		for(int i=kvec_beg;i<kvec_end;i++)
			{
				vmin = CGLA::v_min(vmin,init_nodes[i].key);
				vmax = CGLA::v_max(vmax,init_nodes[i].key);
			}
		int od=0;
		KeyType ave_v = vmax-vmin;
		for(int i=1;i<KeyType::get_dim();i++)
			if(ave_v[i]>ave_v[od]) od = i;
		return od;
	} 

	template<class KeyT, class ValT>
	void KDTree<KeyT,ValT>::optimize(int cur,
																	 int kvec_beg,  
																	 int kvec_end)
	{
		// Assert that we are not inserting beyond capacity.
		assert(cur < nodes.size());

		// If there is just a single element, we simply insert.
		if(kvec_beg+1==kvec_end) 
			{
				nodes[cur] = init_nodes[kvec_beg];
				nodes[cur].dsc = -1;
				return;
			}
	
		// Find the axis that best separates the data.
		int disc = opt_disc(kvec_beg, kvec_end);

		// Compute the median element. See my document on how to do this
		// www.imm.dtu.dk/~jab/publications.html
		int N = kvec_end-kvec_beg;
		int M = 1<< (CGLA::two_to_what_power(N));
		int R = N-(M-1);
		int left_size  = (M-2)/2;
		int right_size = (M-2)/2;
		if(R < M/2)
			{
				left_size += R;
			}
		else
			{
				left_size += M/2;
				right_size += R-M/2;
			}

		int median = kvec_beg + left_size;

		// Sort elements but use nth_element (which is cheaper) than
		// a sorting algorithm. All elements to the left of the median
		// will be smaller than or equal the median. All elements to the right
		// will be greater than or equal to the median.
		const Comp comp(disc);
		std::nth_element(&init_nodes[kvec_beg], 
										 &init_nodes[median], 
										 &init_nodes[kvec_end], comp);

		// Insert the node in the final data structure.
		nodes[cur] = init_nodes[median];
		nodes[cur].dsc = disc;

		// Recursively build left and right tree.
		if(left_size>0)	
			optimize(2*cur, kvec_beg, median);
		
		if(right_size>0) 
			optimize(2*cur+1, median+1, kvec_end);
	}

	template<class KeyT, class ValT>
	int KDTree<KeyT,ValT>::closest_point_priv(int n, const KeyType& p, 
																						ScalarType& dist) const
	{
		int ret_node = 0;
		ScalarType this_dist = nodes[n].dist(p);

		if(this_dist<dist)
			{
				dist = this_dist;
				ret_node = n;
			}
		if(nodes[n].dsc != -1)
			{
				int dsc         = nodes[n].dsc;
				ScalarType dsc_dist  = CGLA::sqr(nodes[n].key[dsc]-p[dsc]);
				bool left_son   = Comp(dsc)(p,nodes[n].key);

				if(left_son||dsc_dist<dist)
					{
						int left_child = 2*n;
						if(left_child < nodes.size())
							if(int nl=closest_point_priv(left_child, p, dist))
								ret_node = nl;
					}
				if(!left_son||dsc_dist<dist)
					{
						int right_child = 2*n+1;
						if(right_child < nodes.size())
							if(int nr=closest_point_priv(right_child, p, dist))
								ret_node = nr;
					}
			}
		return ret_node;
	}

	template<class KeyT, class ValT>
	void KDTree<KeyT,ValT>::in_sphere_priv(int n, 
																				 const KeyType& p, 
																				 const ScalarType& dist,
																				 std::vector<KeyT>& keys,
																				 std::vector<ValT>& vals) const
	{
		ScalarType this_dist = nodes[n].dist(p);
		assert(n<nodes.size());
		if(this_dist<dist)
			{
				keys.push_back(nodes[n].key);
				vals.push_back(nodes[n].val);
			}
		if(nodes[n].dsc != -1)
			{
				const int dsc         = nodes[n].dsc;
				const ScalarType dsc_dist  = CGLA::sqr(nodes[n].key[dsc]-p[dsc]);

				bool left_son = Comp(dsc)(p,nodes[n].key);

				if(left_son||dsc_dist<dist)
					{
						int left_child = 2*n;
						if(left_child < nodes.size())
							in_sphere_priv(left_child, p, dist, keys, vals);
					}
				if(!left_son||dsc_dist<dist)
					{
						int right_child = 2*n+1;
						if(right_child < nodes.size())
							in_sphere_priv(right_child, p, dist, keys, vals);
					}
			}
	}
}
namespace GEO = Geometry;

#if (_MSC_VER >= 1200)
#pragma warning (pop)
#endif


#endif
