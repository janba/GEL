/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file KDTree.h
 * @brief KD Tree. Internal nodes form a contiguous preorder tree; leaves are short point buckets.
 */

#ifndef GEOMETRY_KDTREE_H
#define GEOMETRY_KDTREE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>
#include <GEL/CGLA/CGLA-util.h>
#include <GEL/CGLA/ArithVec.h>

#if (_MSC_VER >= 1200)
#pragma warning (push)
#pragma warning (disable: 4018)
#endif

namespace Geometry
{
namespace detail
{

/** Bounded max-heap of at most `sz` items, stored in the caller's vector.
 The largest item (the one that compares greatest) sits at the front. The
 vector is a heap when it is full, and after `finish` in every case. */
template<class T>
class NQueue {
    std::vector<T>& q;
    const size_t max_sz;

    void sift_down()
    {
        size_t i = 0;
        const size_t n = q.size();
        for (;;) {
            size_t largest = i;
            const size_t left = 2 * i + 1;
            const size_t right = left + 1;
            if (left < n && q[largest] < q[left])
                largest = left;
            if (right < n && q[largest] < q[right])
                largest = right;
            if (largest == i)
                break;
            std::swap(q[i], q[largest]);
            i = largest;
        }
    }

public:
    /// Keep at most `sz` items in `storage`. Previous contents are dropped.
    NQueue(std::vector<T>& storage, size_t sz): q(storage), max_sz(sz)
    {
        q.clear();
        if (max_sz > 0 && q.capacity() < max_sz)
            q.reserve(max_sz);
    }

    /// Push an item. Once full, a worse item is ignored and a better one
    /// replaces the current worst.
    void push(const T& item) {
        if (max_sz == 0)
            return;
        if (q.size() < max_sz) {
            q.push_back(item);
            if (q.size() == max_sz)
                std::make_heap(q.begin(), q.end());
            return;
        }
        if (!(item < q.front()))
            return;
        q.front() = item;
        sift_down();
    }

    size_t size() const {
        return q.size();
    }

    bool at_capacity() const {
        return max_sz > 0 && q.size() == max_sz;
    }

    /// Get the biggest NQueue element. The queue must be full.
    const T& top() const {
        return q.front();
    }

    /// Make `storage` a heap when fewer than `sz` items were inserted.
    void finish() {
        if (q.size() < max_sz)
            std::make_heap(q.begin(), q.end());
    }
};

template<class Key>
inline typename Key::ScalarType sqr_distance(const Key& a, const Key& b)
{
    using Scalar = typename Key::ScalarType;
    const Scalar* aa = a.get();
    const Scalar* bb = b.get();
    Scalar sum = 0;
    constexpr unsigned dim = Key::get_dim();
    for (unsigned i = 0; i < dim; ++i) {
        const Scalar d = aa[i] - bb[i];
        sum += d * d;
    }
    return sum;
}

} // namespace detail

    template<class KeyT, class ValT>
    struct KDTreeRecord {
        using ScalarType = typename KeyT::ScalarType;
        ScalarType d;
        KeyT k;
        ValT v;
        KDTreeRecord(ScalarType _d, KeyT _k, ValT _v): d(_d), k(_k), v(_v) {}
        KDTreeRecord() {}
    };
    
    template<class KeyT, class ValT>
    bool operator < (const KDTreeRecord<KeyT, ValT>& a, const KDTreeRecord<KeyT, ValT>& b) {
        return a.d < b.d;
    }

    template<typename KeyType>
    inline unsigned int comp(int dsc, const KeyType& k0, const KeyType& k1) {
        constexpr unsigned dim=KeyType::get_dim();
        for(unsigned i=0;i<dim;i++)
        {
            unsigned j=(dsc+i)%dim;
            if(k0[j]<k1[j])
                return 0;
            if(k0[j]>k1[j])
                return 1;
        }
        return 1;
    }
    
    
    /** \brief A classic K-D tree.
     
     A K-D tree is a good data structure for storing points in space
     and for nearest neighbour queries. It is basically a generalized
     binary tree in K dimensions. */
    template<class KeyT, class ValT>
    class KDTree
    {
    public:
        typedef typename KeyT::ScalarType ScalarType;
        typedef KeyT KeyType;
        
        /// KDNode struct represents node in KD tree
        struct KDNode
        {
            KeyT key;
            ValT val;
            short dsc = -1;
            
            KDNode() {}
            
            KDNode(const KeyT& _key, const ValT& _val):
            key(_key), val(_val) {}
            
            ScalarType dist(const KeyType& p) const
            {
                return detail::sqr_distance(key, p);
            }
        };
        
        typedef std::vector<KDNode> NodeVecType;
        typedef typename NodeVecType::const_iterator NodeVecConstIterType;
        
    private:
        
        bool is_built;
        NodeVecType init_nodes;
        NodeVecType nodes;
        /// Split once a bucket would hold more points than this.
        /// 16 was the fastest leaf size measured for 3D nearest-neighbor queries.
        static constexpr unsigned leaf_size = 16;
        /// Internal node: index of the right subtree, or 0 when it has none.
        /// Leaf: number of bucket points starting at this node.
        /// Node 0 is unused, so 0 is never a child index.
        std::vector<unsigned> link;
        
        
        /** Passed a vector of keys, this function will construct an optimal tree.
         It is called recursively */
        void optimize(unsigned, unsigned, unsigned);
        
        /** Finde nearest neighbour. */
        unsigned closest_point_priv(unsigned, const KeyType&, ScalarType&) const;
        
        
        template<class Sink>
        void in_sphere_priv(unsigned n,
                            const KeyType& p,
                            const ScalarType& dist,
                            Sink& sink) const;

        /// `make` builds the queue record from a squared distance and a tree node.
        /// The id-only search passes a function that leaves the position key in the tree.
        template<class Record, class Make>
        void m_closest_priv(unsigned n,
                            const KeyType& p,
                            ScalarType& max_dist,
                            detail::NQueue<Record>& nq,
                            Make& make) const;

        template<class Record, class Make>
        void m_closest_into(unsigned m,
                            const KeyType& p,
                            ScalarType dist,
                            std::vector<Record>& out,
                            Make make) const;
        
        /** Finds the optimal discriminator. There are more ways, but this
         function traverses the vector and finds out what dimension has
         the greatest difference between min and max element. That dimension
         is used for discriminator */
        short opt_disc(unsigned,unsigned) const;
        
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
            link.assign(init_nodes.size(), 0);
            if(init_nodes.size() > 1)
                optimize(1,1,static_cast<unsigned>(init_nodes.size()));
            NodeVecType v(1);
            init_nodes.swap(v);
            is_built = true;
        }
        
        NodeVecConstIterType begin() const { return nodes.begin();}
        NodeVecConstIterType end() const {return nodes.end();}
        
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
                if(unsigned n = closest_point_priv(1, p, max_sq_dist))
                {
                    k = nodes[n].key;
                    v = nodes[n].val;
                    dist = std::sqrt(max_sq_dist);
                    return true;
                }
            }
            return false;
        }
        
        /** Replace `keys` and `vals` with the points inside the sphere.
         An empty tree leaves both vectors untouched and returns 0. */
        unsigned in_sphere(const KeyType& p,
                      ScalarType dist,
                      std::vector<KeyT>& keys,
                      std::vector<ValT>& vals) const
        {
            assert(is_built);
            if(nodes.size()>1)
            {
                keys.clear();
                vals.clear();
                ScalarType max_sq_dist = CGLA::sqr(dist);
                auto sink = [&](const KDNode& node) {
                    keys.push_back(node.key);
                    vals.push_back(node.val);
                };
                in_sphere_priv(1, p, max_sq_dist, sink);
                return static_cast<unsigned>(keys.size());
            }
            return 0;
        }

        /// Find all the elements within a given radius (second argument) of
        /// the key (first argument). The key value pairs inside the sphere are
        /// returned in a vector passed as the last argument.
        /// Note that we don't resize the vector to zero - so either
        /// it should be an empty vector or you should desire appending the newly
        /// found elements to this vector.
        unsigned in_sphere(const KeyType& p,
              ScalarType dist,
              std::vector<std::pair<KeyT, ValT>>& kvs) const
        {
            assert(is_built);
            if(nodes.size()>1)
            {
                ScalarType max_sq_dist = CGLA::sqr(dist);
                auto sink = [&](const KDNode& node) {
                    kvs.emplace_back(node.key, node.val);
                };
                in_sphere_priv(1, p, max_sq_dist, sink);
                return static_cast<unsigned>(kvs.size());
            }
            return 0;
        }

        /// Append values within the sphere. Position keys stay in the tree.
        /// The returned count is the number of values appended.
        unsigned in_sphere(const KeyType& p,
                           ScalarType dist,
                           std::vector<ValT>& vals) const
        {
            assert(is_built);
            if(nodes.size()>1)
            {
                const size_t before = vals.size();
                ScalarType max_sq_dist = CGLA::sqr(dist);
                auto sink = [&](const KDNode& node) {
                    vals.push_back(node.val);
                };
                in_sphere_priv(1, p, max_sq_dist, sink);
                return static_cast<unsigned>(vals.size() - before);
            }
            return 0;
        }
        
        /** Find the m elements closest to p and within a distance dist. This function returns a vector
         of KDTreeRecords in no particular order. This function is likely to be faster than simply
         finding all elements within a given radius using in_sphere and then sorting because once m elements have been
         found, the search radius can be narrowed. */
        std::vector<KDTreeRecord<KeyT, ValT>> m_closest(unsigned m, const KeyType& p, ScalarType dist) const {
            std::vector<KDTreeRecord<KeyT,ValT>> nv;
            m_closest(m, p, dist, nv);
            return nv;
        }

        /// Same search as `m_closest`, written into `out` so the caller can reuse its storage.
        void m_closest(unsigned m, const KeyType& p, ScalarType dist,
                       std::vector<KDTreeRecord<KeyT, ValT>>& out) const {
            auto make = [](ScalarType d, const KDNode& node) {
                return KDTreeRecord<KeyT, ValT>(d, node.key, node.val);
            };
            m_closest_into(m, p, dist, out, make);
        }

        /// One nearest-neighbor hit without the position key. `d` is squared distance.
        struct IdRecord {
            ScalarType d = 0;
            ValT v = {};
            IdRecord() = default;
            IdRecord(ScalarType dist, ValT val) : d(dist), v(val) {}
            bool operator<(const IdRecord& other) const { return d < other.d; }
        };

        /// The `m` closest values within `dist`, ordered as a heap by squared distance.
        /// The position keys stay in the tree.
        std::vector<IdRecord> m_closest_ids(unsigned m, const KeyType& p, ScalarType dist) const {
            std::vector<IdRecord> nv;
            m_closest_ids(m, p, dist, nv);
            return nv;
        }

        /// Same search as `m_closest_ids`, written into `out`.
        void m_closest_ids(unsigned m, const KeyType& p, ScalarType dist,
                           std::vector<IdRecord>& out) const {
            auto make = [](ScalarType d, const KDNode& node) {
                return IdRecord(d, node.val);
            };
            m_closest_into(m, p, dist, out, make);
        }
    };
    
    template<class KeyT, class ValT>
    short KDTree<KeyT,ValT>::opt_disc(unsigned kvec_beg,
                                    unsigned kvec_end) const
    {
        constexpr unsigned dim = KeyType::get_dim();
        const ScalarType* first = init_nodes[kvec_beg].key.get();
        std::array<ScalarType, dim> vmin;
        std::array<ScalarType, dim> vmax;
        for (unsigned a = 0; a < dim; ++a) {
            vmin[a] = first[a];
            vmax[a] = first[a];
        }
        for (unsigned i = kvec_beg + 1; i < kvec_end; ++i) {
            const ScalarType* coord = init_nodes[i].key.get();
            for (unsigned a = 0; a < dim; ++a) {
                const ScalarType c = coord[a];
                if (c < vmin[a])
                    vmin[a] = c;
                else if (c > vmax[a])
                    vmax[a] = c;
            }
        }
        short axis = 0;
        ScalarType best = vmax[0] - vmin[0];
        for (unsigned a = 1; a < dim; ++a) {
            const ScalarType span = vmax[a] - vmin[a];
            if (span > best) {
                best = span;
                axis = static_cast<short>(a);
            }
        }
        return axis;
    }
    
    template<class KeyT, class ValT>
    void KDTree<KeyT,ValT>::optimize(unsigned cur,
                                     unsigned kvec_beg,
                                     unsigned kvec_end)
{
        const unsigned count = kvec_end - kvec_beg;
        assert(cur + count <= nodes.size());

        // A short contiguous run is cheaper to scan than to search.
        if (count <= leaf_size)
        {
            for (unsigned i = 0; i < count; ++i) {
                nodes[cur + i] = init_nodes[kvec_beg + i];
                nodes[cur + i].dsc = -1;
            }
            link[cur] = count;
            return;
        }

        // Widest axis. The split only needs that one coordinate.
        const short disc = opt_disc(kvec_beg, kvec_end);

        // Left subtree is the half of the range before the median, stored
        // immediately after this node. Every internal node therefore has a
        // left child at cur+1. The right subtree, if any, follows the left.
        const unsigned left_size = count / 2;
        const unsigned median = kvec_beg + left_size;
        const unsigned right_size = count - left_size - 1;

        auto compare_func = [disc](const KDNode& a, const KDNode& b) {
            return a.key[disc] < b.key[disc];
        };
        KDNode* data = init_nodes.data();
        std::nth_element(data + kvec_beg, data + median, data + kvec_end, compare_func);

        nodes[cur] = init_nodes[median];
        nodes[cur].dsc = disc;

        optimize(cur + 1, kvec_beg, median);
        if (right_size > 0) {
            const unsigned right = cur + 1 + left_size;
            link[cur] = right;
            optimize(right, median + 1, kvec_end);
        }
    }
    
    template<class KeyT, class ValT>
    unsigned KDTree<KeyT,ValT>::closest_point_priv(unsigned n, const KeyType& p,
                                              ScalarType& dist) const
    {
        const short dsc = nodes[n].dsc;
        if (dsc < 0) {
            unsigned ret_node = 0;
            const unsigned count = link[n];
            const KDNode* leaf = nodes.data() + n;
            for (unsigned i = 0; i < count; ++i) {
                const ScalarType this_dist = detail::sqr_distance(leaf[i].key, p);
                if (this_dist < dist) {
                    dist = this_dist;
                    ret_node = n + i;
                }
            }
            return ret_node;
        }

        unsigned ret_node = 0;
        const ScalarType this_dist = detail::sqr_distance(nodes[n].key, p);
        if(this_dist<dist)
        {
            dist = this_dist;
            ret_node = n;
        }
        const ScalarType split = nodes[n].key[dsc];
        const ScalarType plane = CGLA::sqr(split - p[dsc]);
        const unsigned right = link[n];
        if (p[dsc] < split) {
            if (unsigned found = closest_point_priv(n + 1, p, dist))
                ret_node = found;
            if (right != 0 && plane < dist)
                if (unsigned found = closest_point_priv(right, p, dist))
                    ret_node = found;
        } else {
            if (right != 0)
                if (unsigned found = closest_point_priv(right, p, dist))
                    ret_node = found;
            if (plane < dist)
                if (unsigned found = closest_point_priv(n + 1, p, dist))
                    ret_node = found;
        }
        return ret_node;
    }

    template<class KeyT, class ValT>
    template<class Sink>
    void KDTree<KeyT,ValT>::in_sphere_priv(unsigned n,
                                           const KeyType& p,
                                           const ScalarType& dist,
                                           Sink& sink) const
    {
        const short dsc = nodes[n].dsc;
        if (dsc < 0) {
            const unsigned count = link[n];
            const KDNode* leaf = nodes.data() + n;
            for (unsigned i = 0; i < count; ++i) {
                if (detail::sqr_distance(leaf[i].key, p) < dist)
                    sink(leaf[i]);
            }
            return;
        }
        const ScalarType this_dist = detail::sqr_distance(nodes[n].key, p);
        if(this_dist<dist)
            sink(nodes[n]);
        const ScalarType split = nodes[n].key[dsc];
        const ScalarType plane = CGLA::sqr(split - p[dsc]);
        const unsigned right = link[n];
        if (p[dsc] < split) {
            in_sphere_priv(n + 1, p, dist, sink);
            if (right != 0 && plane < dist)
                in_sphere_priv(right, p, dist, sink);
        } else {
            if (right != 0)
                in_sphere_priv(right, p, dist, sink);
            if (plane < dist)
                in_sphere_priv(n + 1, p, dist, sink);
        }
    }

    template<class KeyT, class ValT>
    template<class Record, class Make>
    void KDTree<KeyT,ValT>::m_closest_priv(unsigned n,
                                           const KeyType& p,
                                           ScalarType& max_dist,
                                           detail::NQueue<Record>& nq,
                                           Make& make) const
    {
        const short dsc = nodes[n].dsc;
        if (dsc < 0) {
            const unsigned count = link[n];
            const KDNode* leaf = nodes.data() + n;
            for (unsigned i = 0; i < count; ++i) {
                const ScalarType dist = detail::sqr_distance(leaf[i].key, p);
                if (dist < max_dist) {
                    nq.push(make(dist, leaf[i]));
                    if (nq.at_capacity())
                        max_dist = nq.top().d;
                }
            }
            return;
        }
        const ScalarType dist = detail::sqr_distance(nodes[n].key, p);
        if(dist<max_dist) {
            nq.push(make(dist, nodes[n]));
            if(nq.at_capacity())
                max_dist = nq.top().d;
        }
        const ScalarType split = nodes[n].key[dsc];
        const ScalarType plane = CGLA::sqr(split - p[dsc]);
        const unsigned right = link[n];
        if (p[dsc] < split) {
            m_closest_priv(n + 1, p, max_dist, nq, make);
            if (right != 0 && plane < max_dist)
                m_closest_priv(right, p, max_dist, nq, make);
        } else {
            if (right != 0)
                m_closest_priv(right, p, max_dist, nq, make);
            if (plane < max_dist)
                m_closest_priv(n + 1, p, max_dist, nq, make);
        }
    }

    template<class KeyT, class ValT>
    template<class Record, class Make>
    void KDTree<KeyT,ValT>::m_closest_into(unsigned m,
                                           const KeyType& p,
                                           ScalarType dist,
                                           std::vector<Record>& out,
                                           Make make) const
    {
        assert(is_built);
        detail::NQueue<Record> nq(out, m);
        if (nodes.size() > 1 && m > 0) {
            ScalarType max_sq_dist = CGLA::sqr(dist);
            m_closest_priv(1, p, max_sq_dist, nq, make);
        }
        nq.finish();
    }
}
namespace GEO = Geometry;

#if (_MSC_VER >= 1200)
#pragma warning (pop)
#endif


#endif
