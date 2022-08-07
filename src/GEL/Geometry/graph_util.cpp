//
//  graph_abstraction.cpp
//  MeshEditE
//
//  Created by Jakob Andreas Bærentzen on 30/04/2018.
//  Copyright © 2018 J. Andreas Bærentzen. All rights reserved.
//

#include <future>
#include <thread>
#include <unordered_set>
#include <queue>
#include <vector>
#include <iostream>
#include <random>
#include <GEL/Util/AttribVec.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/build_bbtree.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/Geometry/GridAlgorithm.h>
#include <GEL/Geometry/bounding_sphere.h>
#include <GEL/HMesh/HMesh.h>
#include <GEL/Geometry/graph_util.h>

using namespace std;
using namespace CGLA;
using namespace Util;
using namespace HMesh;

namespace Geometry {

    using NodeID = AMGraph::NodeID;
    using NodeSet = AMGraph::NodeSet;
    using NodeSetUnordered = unordered_set<NodeID>;
    using NodeQueue = queue<NodeID>;
    using NodeSetVec = vector<pair<double,NodeSet>>;



    struct SkeletonPQElem {
        double pri;
        Geometry::AMGraph::NodeID n0, n1;
        SkeletonPQElem(double _pri, Geometry::AMGraph::NodeID _n0, Geometry::AMGraph::NodeID _n1);
        
        bool operator < (const SkeletonPQElem& pq1) const {
            return pri < pq1.pri;
        }
    };
    SkeletonPQElem::SkeletonPQElem(double _pri, NodeID _n0, NodeID _n1): pri(_pri), n0(_n0), n1(_n1) {}


    int test_intersection (const NodeSet& set1, const NodeSet& set2)
    {
        auto first1 = begin(set1);
        auto first2 = begin(set2);
        int matches=0;
        while (first1!=end(set1) && first2!=end(set2))
        {
            if (*first1<*first2) ++first1;
            else if (*first2<*first1) ++first2;
            else {
                ++matches;
                ++first1;
                ++first2;
            }
        }
        return matches;
    }

    NodeSetUnordered neighbors(const AMGraph3D& g, const NodeSetUnordered& s) {
        NodeSetUnordered _nbors;
        for(auto n: s)
            for(auto m: g.neighbors(n))
                _nbors.insert(m);
        
        NodeSetUnordered nbors;
        for(auto n: _nbors)
            if (s.count(n) == 0)
                nbors.insert(n);
        
        return nbors;
    }



    NodeSet order(NodeSetUnordered& s) {
        NodeSet _s;
        for(auto n : s)
            _s.insert(n);
        return _s;
    }


    Vec3d node_set_barycenter(const AMGraph3D& g, const NodeSet& ns) {
        Vec3d b(0);
        for(auto n: ns)
            b += g.pos[n];
        return b/ns.size();
    }

    NodeSet nodeset_dilation(const AMGraph3D& g, const NodeSet& s,
                             const NodeSet& excl) {
        
        NodeSet out = s;
        for(auto n : s) {
            for(auto m: g.neighbors(n)) {
                if(excl.count(m)==0 && out.count(m)==0)
                    out.insert(m);
            }
        }
        return out;
    }

    std::vector<NodeSetUnordered> connected_components(const AMGraph& g,
                                                       const NodeSetUnordered& s) {
        NodeSetUnordered s_visited;
        vector<NodeSetUnordered> component_vec;
        for(auto nf0 : s) {
            if(s_visited.count(nf0)==0)
            {
                NodeQueue Q;
                Q.push(nf0);
                s_visited.insert(nf0);
                NodeSetUnordered component;
                while(!Q.empty())
                {
                    NodeID nf = Q.front();
                    Q.pop();
                    
                    bool in_s = s.count(nf)>0;
                    if(in_s)
                        component.insert(nf);
                    for(auto nnf: g.neighbors(nf)) {
                        if (s_visited.count(nnf)==0) {
                            if( s.count(nnf)>0 ) {
                                Q.push(nnf);
                                s_visited.insert(nnf);
                            }
                        }
                    }
                    
                }
                component_vec.push_back(component);
            }
        }
        return component_vec;
    }

    void saturate_graph(AMGraph3D& g, int hops, double dist_frac, double rad) {
        AMGraph3D g2 = g;
        using NodeMap = std::map<NodeID, pair<int, double>>;
        for(NodeID n0: g.node_ids()) {
            queue<NodeID> Q;
            Q.push(n0);
            NodeMap node_map;
            node_map[n0] = make_pair(0, 0.0);
            
            while(!Q.empty()) {
                auto n = Q.front();
                auto [h_n, d_n] = node_map[n];
                Q.pop();
                for(auto m: g.neighbors(n)) {
                    double d_m = d_n + sqrt(g.sqr_dist(n, m));
                    if(node_map.count(m) == 0 || d_m < node_map[m].second) {
                        double d_n0_m = sqrt(g.sqr_dist(n0, m));
                        if (d_n0_m<dist_frac*d_m && d_n0_m < rad)
                            g2.connect_nodes(n0, m);
                        if(h_n+1<hops)
                            Q.push(m);
                        node_map[m] = make_pair(h_n+1, d_m);
                    }
                }
            }
        }
        g = g2;
    }

    Vec3d geometric_median(const vector<Vec3d>& pts) {
        Vec3d y(0.0);
        for(const auto& p: pts)
            y += p;
        y /= pts.size();
        for(int _ = 0;_<100;++_) {
            Vec3d yp(0);
            double wsum = 0;
            for(const auto& p: pts) {
                double w = 1.0/length(p-y);
                yp += p*w;
                wsum += w;
            }
            yp /= wsum;
            double err = sqr_length(y-yp);
            y = yp;
            if(err<1e-10)
                break;
        }
        return y;
    }


    void smooth_graph(AMGraph3D& g, const int iter, const float alpha) {
        auto lsmooth = [](AMGraph3D& g, float _alpha){
            AttribVec<AMGraph::NodeID, Vec3d> new_pos(g.no_nodes(), Vec3d(0));
            for(auto n: g.node_ids()) {
                double wsum = 0;
                auto N = g.neighbors(n);
                for(auto nn: N) {
                    double w = 1.0;
                    new_pos[n] += w*g.pos[nn];
                    wsum += w;
                }
                double alpha = N.size()==1 ? 0 : _alpha;
                new_pos[n] = (alpha) * new_pos[n] / wsum + (1.0-alpha) * g.pos[n];
            }
            return new_pos;
        };
        
        for(int i = 0;i<iter;++i) {
            auto npos = lsmooth(g, alpha);
            g.pos = npos;
        }
    }

    int graph_edge_contract(AMGraph3D& g, double dist_thresh) {
        using NodeID = AMGraph::NodeID;
        auto priority = [&](NodeID a, NodeID b) { return -g.sqr_dist(a,b); };
        priority_queue<SkeletonPQElem> Q;
        
        int cntr, total_work = 0;
        do {
            Util::AttribVec<AMGraph::NodeID, int> touched(g.no_nodes(),0);
            cntr = 0;
            for(auto n0 : g.node_ids()) {
                for(auto n1: g.neighbors(n0)) {
                    double pri = priority(n0,n1);
                    if(pri>-sqr(dist_thresh))
                        Q.push(SkeletonPQElem(pri, n0, n1));
                }
            }
            
            while(!Q.empty()) { 
                auto skel_rec = Q.top();
                Q.pop();
                if(touched[skel_rec.n0]==0 && touched[skel_rec.n1]==0) {
                    auto e = g.find_edge(skel_rec.n0, skel_rec.n1);
                    if( e != AMGraph::InvalidEdgeID) {
                        g.merge_nodes(skel_rec.n0,skel_rec.n1, true);
                        touched[skel_rec.n0] = 1;
                        touched[skel_rec.n1] = 1;
                        ++cntr;
                    }
                }
            }
            total_work += cntr;
        } while(cntr);
        
        g.cleanup();
        return total_work;
    }

    namespace  {
    const Vec3f& get_color(int i)
    {
        static Vec3f ctable[100000];
        static bool was_here;
        gel_srand(0);
        if(!was_here)
        {
            was_here = true;
            ctable[0] = Vec3f(0);
            for(int j=1;j<100000;++j)
                ctable[j] = Vec3f(0.3)+0.7*normalize(Vec3f(gel_rand(),gel_rand(),gel_rand()));
            ctable[3] = Vec3f(1,0,0);
            ctable[4] = Vec3f(0,1,0);
            ctable[5] = Vec3f(0,0,1);
            ctable[6] = Vec3f(1,0,1);
        }
        return ctable[i%100000];
    }
    }

    void color_graph_node_sets(AMGraph3D& g, const NodeSetVec& node_set_vec) {
        double w_max = 0;
        for(const auto& [w,ns]:node_set_vec)
            w_max = max(w,w_max);
        for(auto n : g.node_ids())
            g.node_color[n] = Vec3f(0.8);
        for(auto e : g.edge_ids())
            g.edge_color[e] = Vec3f(0.8);
        int cnt = 0;
        for(const auto& [w,ns]:node_set_vec) {
            //        cout << " w, sz:" << w << " , " << ns.size() << endl;
            Vec3f col  = sqr(get_color(cnt++));
            //        Vec3f col = Vec3f(w,0,1-w);// = GLGraphics::get_color(cnt);
            for(auto n: ns) {
                g.node_color[n] = col;
                for(auto m: ns)
                    if( n!= m) {
                        auto e = g.find_edge(n, m);
                        if (g.valid_edge_id(e))
                            g.edge_color[e] = col;
                    }
            }
        }
    }

    AttribVec<NodeID, double> leaf_distance(const AMGraph3D& g) {
        BreadthFirstSearch bfs(g);
        for(auto n: g.node_ids()) {
            if(g.neighbors(n).size()<=1)
                bfs.add_init_node(n,0);
        }
        while(bfs.Dijkstra_step());
        return bfs.dist;
    }


    AttribVecDouble smooth_dist(const AMGraph3D& g, const AttribVecDouble& _dist, int smooth_iter) {
        AttribVecDouble dist = _dist;
        AttribVecDouble dist_new = dist;
        double wgt = 0.4;
        for(int i=0;i<smooth_iter;++i) {
            for(auto n: g.node_ids()) {
                const auto& N =  g.neighbors(n);
                bool is_max=true;
                bool is_min=true;
                double d_nn_sum=0;
                for (auto nn: N){
                    if(dist[nn]<dist[n])
                        is_min=false;
                    else if(dist[nn]>dist[n])
                        is_max=false;
                    d_nn_sum += dist[nn];
                }
                if(!(is_max || is_min)) {
                    dist_new[n] *= (1.0-wgt);
                    dist_new[n] += d_nn_sum * wgt/N.size();
                }
            }
            dist=dist_new;
        }
        return dist;
    }


    AttribVecDouble projection(const AMGraph3D& g, const Vec3d& dir, int smooth_iter) {
        AttribVecDouble dist;
        for(auto n: g.node_ids()) {
            dist[n] = dot(dir, g.pos[n]);
        }
        return smooth_dist(g,dist,smooth_iter);
    }

    AttribVecDouble negate_dist(const AMGraph3D& g, const AttribVecDouble& dist_in) {
        AttribVecDouble dist_out;
        for(auto n: g.node_ids()) {
            dist_out[n] = -dist_in[n];
        }
        return dist_out;
    }


    void prune(Geometry::AMGraph3D& g) {
        vector<NodeID> garbage;
        for(auto n: g.node_ids()) {
            auto N = g.neighbors(n);
            if (N.size()==1)
            {
                auto m = N[0];
                auto M = g.neighbors(m);
                if(M.size()>2)
                    garbage.push_back(n);
            }
            if(N.size()==0)
                garbage.push_back(n);
        }
        for(auto n: garbage)
            g.remove_node(n);
        
        g = clean_graph(g);
    }





    AMGraph3D voxel_graph_from_mesh(Manifold& m, int res) {
        Vec3d p0,p7;
        bbox(m, p0, p7);
        Vec3d diag = p7-p0;
        double D = diag.min_coord();
        double l = diag.max_coord()/res;
        Vec3i dim = Vec3i(diag/l)+Vec3i(1);
        
        AMGraph3D g;
        OBBTree tree;
        build_OBBTree(m, tree);
        RGrid<NodeID> node_grid(dim,-1);
        for(Vec3i p: Range3D(dim))
        {
            Vec3d x = p0 + Vec3d(p)*l;
            double d = tree.compute_signed_distance(Vec3f(x));
            if(d<=0.0) {
                NodeID n =g.add_node(x);
                node_grid[p] = n;
                g.node_color[n] = Vec3f(-4.0*d/D,0,0);
            }
        }
        
        for(Vec3i p: Range3D(dim))
            if(node_grid[p] != -1)
                for(Vec3i nbor : Range3D(Vec3i(-1), Vec3i(2)))
                    if(nbor != Vec3i(0)) {
                        Vec3i pn = p + nbor;
                        if(node_grid.in_domain(pn) && node_grid[pn] != -1)
                            g.connect_nodes(node_grid[p], node_grid[pn]);
                    }
        return g;
    }


    pair<Vec3d, double> approximate_bounding_sphere(const AMGraph3D& g, const NodeSetUnordered& s) {
        vector<Vec3d> pts;
        if(s.empty())
            for(auto n: g.node_ids())
                pts.push_back(g.pos[n]);
        else
            for(auto n : s)
                pts.push_back(g.pos[n]);
        return approximate_bounding_sphere(pts);
    }



    template<typename IndexType, typename  FuncType>
    auto minimum(IndexType b, IndexType e, FuncType f)  {
        using ValueType = decltype(f(*b));
        IndexType min_idx = b;
        ValueType min_val = f(*b);
        for(IndexType i=++b; i != e; ++i) {
            ValueType val = f(*i);
            if (val < min_val) {
                min_val = val;
                min_idx = i;
            }
        }
        return make_pair(min_idx, min_val);
    }


    NodeSetVec k_means_node_clusters(AMGraph3D& g, int N, int MAX_ITER) {
        vector<NodeID> seeds(begin(g.node_ids()), end(g.node_ids()));
        srand(0);
        shuffle(begin(seeds), end(seeds), default_random_engine(rand()));
        seeds.resize(N);
        vector<Vec3d> seed_pos(N);
        for(int i=0;i<N;++i)
            seed_pos[i] = g.pos[seeds[i]];
        
        NodeSetVec clusters(N);
        for (int iter=0;iter<MAX_ITER;++iter) {
            KDTree<Vec3d, size_t> seed_tree;
            for(int i=0;i<N;++i)
                seed_tree.insert(seed_pos[i], i);
            seed_tree.build();
            clusters.clear();
            clusters.resize(N);
            for(auto n: g.node_ids())
            {
                double dist = DBL_MAX;
                Vec3d k;
                size_t idx;
                if(seed_tree.closest_point(g.pos[n], dist, k, idx))
                    clusters[idx].second.insert(n);
            }
            
            multimap<double, NodeID> raw_seeds;
            map<NodeID, Vec3d> cluster_center;
            for(auto [_, ns]: clusters) {
                auto nsc_vec = connected_components(g, ns);
                for(auto nsc: nsc_vec) {
                    Vec3d avg_pos = accumulate(begin(nsc), end(nsc), Vec3d(0.0), [&](const Vec3d& a, NodeID n){return a+g.pos[n];});
                    avg_pos /= nsc.size();
                    
                    size_t b_cnt = 0;
                    for(auto n: nsc)
                        for(auto m: g.neighbors(n))
                            if (nsc.count(m) == 0)
                                b_cnt += 1;
                    
                    auto [node_iter, min_dist] = minimum(begin(nsc), end(nsc), [&](NodeID n){return length(g.pos[n]-avg_pos);});
                    raw_seeds.insert(make_pair(b_cnt/double(nsc.size()), *node_iter));
                    cluster_center[*node_iter] = avg_pos;
                }
            }
            auto node_iter = raw_seeds.begin();
            for(int i=0;i<N;++i, ++node_iter) {
                NodeID n = node_iter->second;
                seeds[i] = n;
                seed_pos[i] = cluster_center[n];
            }
            
        }
        
        color_graph_node_sets(g, clusters);
        return clusters;
    }

    GraphDist::GraphDist(const AMGraph3D& g) {
        for(auto n : g.node_ids())
            if(g.in_use(n))
                for(auto m : g.neighbors(n))
                    if(n<m) {
                        Vec3d pn = g.pos[n];
                        Vec3d pm = g.pos[m];
                        Vec3d pc = (pn+pm)/2.0;
                        seg_tree.insert(pc, segments.size());
                        segments.push_back(LineSegment(pn,pm));
                        R = max(R, length(pc-pn));
                    }
        seg_tree.build();
    }

    double GraphDist::dist(const Vec3d& p) {
        double dist = 1e32;
        Vec3d k;
        size_t segment_idx;
        if(seg_tree.closest_point(p, dist, k, segment_idx)) {
            dist = sqrt(segments[segment_idx].sqr_distance(p).sqr_dist);
            dist += R;
            vector<Vec3d> keys;
            vector<size_t> vals;
            seg_tree.in_sphere(p, dist, keys, vals);
            for (size_t idx: vals)
                dist = min(dist, sqrt(segments[idx].sqr_distance(p).sqr_dist));
        }
        return dist;
    }

    pair<double,double> graph_H_dist(const Geometry::AMGraph3D& g0, const Geometry::AMGraph3D& g, size_t samples) {
        
        GraphDist gd0(g0);
        
        double total_length = 0;
        for(auto n : g.node_ids())
            for(auto m : g.neighbors(n))
                if(g.valid_node_id(n) && g.valid_node_id(m) && n<m) {
                    total_length += length(g.pos[m]-g.pos[n]);
                }
        int cnt = 0;
        double avg_dist = 0;
        double max_dist = 0.0;
        srand(0);
        for(auto n : g.node_ids())
            for(auto m : g.neighbors(n))
                if(g.valid_node_id(n) && g.valid_node_id(m) && n<m) {
                    double l = length(g.pos[m]-g.pos[n]);
                    int samples_per_edge = samples*(l/total_length) + 0.5;
                    for (int s = 0; s < samples_per_edge; ++s) {
                        double r = rand()/double(RAND_MAX);
                        Vec3d p = r * g.pos[m] + (1.0-r) * g.pos[n];
                        double d = gd0.dist(p);
                        avg_dist += d;
                        max_dist = max(max_dist,d);
                        cnt += 1;
                    }
                }
        avg_dist /= cnt;
        return make_pair(avg_dist, max_dist);
        
    }
}
