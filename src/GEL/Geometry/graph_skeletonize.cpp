//
//  graph_skeletonize.cpp
//  MeshEditE
//
//  Created by Andreas Bærentzen on 26/01/2020.
//  Copyright © 2020 J. Andreas Bærentzen. All rights reserved.
//

#include <future>
#include <thread>
#include <unordered_set>
#include <queue>
#include <list>
#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <GEL/Util/AttribVec.h>
#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/graph_util.h>

using namespace std;
using namespace CGLA;
using namespace Util;
using namespace Geometry;

namespace  Geometry {

    using NodeID = AMGraph::NodeID;
    using NodeSetUnordered = unordered_set<NodeID>;
    using NodeQueue = queue<NodeID>;

    #include <GEL/Geometry/graph_util.h>
    #include <GEL/Geometry/graph_skeletonize.h>


    void greedy_weighted_packing(const AMGraph3D& g, NodeSetVec& node_set_vec, bool normalize) {
       
        vector<pair<double,int>> node_set_index;
        
        if(normalize) {
            vector<double> opportunity_cost(node_set_vec.size(), 0.0);
            
            for(int i=0;i<node_set_vec.size();++i) {
                const auto& [w_i,ns_i] = node_set_vec[i];
                for (int j=i+1;j<node_set_vec.size();++j) {
                        const auto& [w_j,ns_j] = node_set_vec[j];
                        int matches = test_intersection(ns_i,ns_j);
                        if(matches >0){
                            opportunity_cost[i] += w_j;
                            opportunity_cost[j] += w_i;
                        }
                    }
                double weight = w_i/opportunity_cost[i];
                node_set_index.push_back(make_pair(weight, i));
            }
        }
        else for(int i=0;i<node_set_vec.size();++i) {
            const auto& [w_i,ns_i] = node_set_vec[i];
            node_set_index.push_back(make_pair(w_i, i));
        }
        
        sort(begin(node_set_index), end(node_set_index), greater<pair<double,int>>());
        NodeSetVec node_set_vec_new;
        AttribVec<NodeID, size_t> set_index(g.no_nodes(), -1);
        for(const auto& [norm_weight, ns_idx]: node_set_index) {
            const auto& [weight,node_set] = node_set_vec[ns_idx];
            bool immaculate = true;
            for(auto n: node_set)
                if (set_index[n] != -1) {
                    immaculate = false;
                    break;
                }
            if (immaculate) {
                node_set_vec_new.push_back(node_set_vec[ns_idx]);
                for(auto n: node_set)
                    set_index[n] = ns_idx;
            }
        }
        swap(node_set_vec_new,node_set_vec);
    }



    NodeSetVec maximize_node_set_vec(AMGraph3D& g, const NodeSetVec& _node_set_vec) {
        NodeSetVec node_set_vec = _node_set_vec;
        
        BreadthFirstSearch bfs(g);
        AttribVec<NodeID, int> nsv_membership(g.no_nodes(),-1);
        for (int nsv_cnt = 0; nsv_cnt < node_set_vec.size(); ++nsv_cnt) {
            const auto& nsv = node_set_vec[nsv_cnt];
            for(auto n: nsv.second) {
                bfs.add_init_node(n);
                nsv_membership[n] = nsv_cnt;
            }
        }

        while(bfs.Dijkstra_step());

        for(auto n: g.node_ids())
            if(nsv_membership[n]==-1 && bfs.pred[n] != AMGraph::InvalidNodeID) {
                auto m = n;
                vector<NodeID> path;
                while(nsv_membership[m]==-1) {
                    path.push_back(m);
                    m = bfs.pred[m];
                }
                auto nsv_number = nsv_membership[m];
                for(auto l: path)
                    nsv_membership[l] = nsv_number;
            }
        
        for(auto n: g.node_ids())
            if (nsv_membership[n] != -1) {
                node_set_vec[nsv_membership[n]].second.insert(n);
            }
        return node_set_vec;
    }



    pair<AMGraph3D, Util::AttribVec<AMGraph::NodeID,AMGraph::NodeID>>
    skeleton_from_node_set_vec(AMGraph3D& g, const NodeSetVec& _node_set_vec, bool merge_branch_nodes, int smooth_steps) {
        // First expand the node_set_vec so that all nodes are assigned.
        NodeSetVec node_set_vec = maximize_node_set_vec(g, _node_set_vec);
    //    color_graph_node_sets(g, node_set_vec);
    //    return make_pair(g, AttribVec<NodeID, NodeID> ());
            
        // Skeleton graph
        AMGraph3D skel;
        Util::AttribVec<AMGraph::NodeID, double> node_size;
        
        // Map from g nodes to skeleton nodes
        AttribVec<NodeID, NodeID> skel_node_map(g.no_nodes(),AMGraph::InvalidNodeID);
        
        // Map from skeleton node to its weight.
        AttribVec<NodeID, double> skel_node_weight;
        
        // Create a skeleton node for each node set.
        for(const auto& [w,ns]: node_set_vec)
            if(ns.size()>0)
            {
                const NodeID skel_node = skel.add_node(Vec3d(0));
                Vec3d avg_pos(0);
                for (auto n: ns) {
                    avg_pos += g.pos[n];
                    skel_node_map[n] = skel_node;
                }
                avg_pos /= ns.size();
                
                vector<double> lengths;
                for(auto n: ns)
                    lengths.push_back(length(g.pos[n]-avg_pos));
                nth_element(begin(lengths), begin(lengths)+lengths.size()/2, end(lengths));
                node_size[skel_node] = lengths[lengths.size()/2];
                skel.pos[skel_node] = avg_pos;
                skel_node_weight[skel_node] = (ns.size());
            }
        
        // If two graph nodes are connected and belong to different skeleton nodes,
        // we also connect their respective skeleton nodes.
        for(NodeID n0: g.node_ids())
            for(NodeID m : g.neighbors(n0)) {
                NodeID skel_node_n0 = skel_node_map[n0];
                NodeID skel_node_m = skel_node_map[m];
                if(skel_node_m != skel_node_n0) {
                    skel.connect_nodes(skel_node_n0, skel_node_m);
                }
            }
        
        // At this point, we return if the merging of branch nodes is not desired.
        if(!merge_branch_nodes)
            return make_pair(skel, skel_node_map);
        
        // If skeletal nodes s0, s1, and s2 form a clique, we add them to the cliques
        // vector of NodeSets.
        vector<NodeSet> cliques;
        for(NodeID s0: skel.node_ids()) {
            auto N_s0 = skel.neighbors(s0);
            for(NodeID s1: N_s0)
                for(NodeID s2: N_s0)
                    if(s1 != s2 && skel.find_edge(s1, s2) != AMGraph::InvalidEdgeID)
                        cliques.push_back({s0,s1,s2});
        }
        
    #define MULTI_PASS_SKELETONIZATION 0
        
    #if MULTI_PASS_SKELETONIZATION
        vector<pair<int,int>> adjacent_cliques;
    #endif
        // If two cliques intersect with more than a single node, we join them.
        for(int i = 0; i< cliques.size(); ++i)
            for(int j = 0; j< cliques.size(); ++j)
                if (i != j) {
                    if(test_intersection(cliques[i], cliques[j])>1) {
    #if MULTI_PASS_SKELETONIZATION
                        adjacent_cliques.push_back(make_pair(i, j));
    #else
                        cliques[i].insert(begin(cliques[j]),end(cliques[j]));
                        cliques[j].clear();
    #endif
                    }
                }
        
        // Now, we create a branch node connected to all of the nodes in the
        // merged clique
        vector<NodeID> branch_nodes;
        for(auto& ns: cliques)
    #if MULTI_PASS_SKELETONIZATION
            if(ns.size()>0)
    #endif
            {
                Vec3d avg_pos(0);
                double wsum = 0;
                double rad=0;
                for(auto n: ns) {
                    avg_pos += skel_node_weight[n]*skel.pos[n];
                    rad += skel_node_weight[n]*node_size[n];
                    wsum += skel_node_weight[n];
                }
                avg_pos /= wsum;
                rad /= wsum;
                auto n_branch = skel.add_node(avg_pos);
                branch_nodes.push_back(n_branch);
                skel.node_color[n_branch] = Vec3f(1,0,0);
                node_size[n_branch] = rad;
                skel_node_weight[n_branch] = wsum / ns.size();
                for(auto n: ns)
                    skel.connect_nodes(n_branch, n);
            }
        
        
        // Disconnect all of the nodes that are now connected to a
        // common branch node.
        for(auto n: branch_nodes) {
            const auto& N = skel.neighbors(n);
            for(auto nn: N)
                for(auto nm: N)
                    skel.disconnect_nodes(nn, nm);
        }
        
    #if MULTI_PASS_SKELETONIZATION
        for(auto np: adjacent_cliques) {
            skel.connect_nodes(branch_nodes[np.first], branch_nodes[np.second]);
        }
    #endif

        // Smooth gently
        for(int iter=0;iter< smooth_steps;++iter) {
            auto skel_new_pos = skel.pos;
            for(auto sn: skel.node_ids()) {
                skel_new_pos[sn] = Vec3d(0);
                double w_sum = 0;
                for (auto nsn: skel.neighbors(sn)) {
                    double w = sqrt(skel_node_weight[nsn])/skel.neighbors(nsn).size();
                    skel_new_pos[sn] += w*skel.pos[nsn];
                    w_sum += w;
                }
                skel_new_pos[sn] /= w_sum;
            }
            
            for(auto sn: skel.node_ids()) {
                double w = 1.0/skel.neighbors(sn).size();
                skel.pos[sn] = skel.pos[sn] * w + (1.0-w) * skel_new_pos[sn];
            }
        }

        // Finally, store the node size in the green channel of the node color
        // it is perhaps not the best idea, but this way we do not need another
        // way of storing the size....
        for(auto n: skel.node_ids())
            skel.node_color[n][1]=node_size[n];
        
        return make_pair(skel, skel_node_map);
    }

    AttribVec<NodeID, double> junction_distance(const AMGraph3D& g) {
        BreadthFirstSearch bfs(g);
        for(auto n: g.node_ids()) {
            if(g.neighbors(n).size()>2)
                bfs.add_init_node(n,0);
        }
        while(bfs.Dijkstra_step());
        return bfs.dist;
    }

    NodeSetVec skeletal_reweighting(AMGraph3D& g, const NodeSetVec& nsv_for_skel) {

        auto [skel,_] = skeleton_from_node_set_vec(g, nsv_for_skel, true, 0);
        auto leaf_dist = junction_distance(skel);
        NodeSetVec nsv;
        for(int i=0;i<nsv_for_skel.size();++i) {
            const auto& [w,ns] = nsv_for_skel[i];
            double l = leaf_dist[NodeID(i)];
            nsv.push_back(make_pair(sqrt(l+1)*w, ns));
        }
        return nsv;
    }

    NodeSetVec separating_node_sets(AMGraph3D& g, const AttribVec<NodeID, double>& dist, int shift) {
        
        BreadthFirstSearch bfs(g,dist);
        while (bfs.step());
        vector<pair<int,NodeID>> nodes_by_tin;
        
        for(auto n : g.node_ids())
            nodes_by_tin.push_back(make_pair(bfs.T_in[n], n));
    //    sort(begin(nodes_by_tin), end(nodes_by_tin));
        shuffle(begin(nodes_by_tin), end(nodes_by_tin), default_random_engine(rand()));

       
        vector<vector<NodeID>> separators;
        AttribVec<NodeID, int> separator_idx(g.no_nodes(), -1);
        int first_T0 = shift;
        for(const auto [T0,n0] : nodes_by_tin)
            if(/*T0 >= first_T0 &&*/ separator_idx[n0]==-1){
                int new_sep_idx = separators.size();
                separator_idx[n0] = new_sep_idx;
                queue<NodeID> nq;
                nq.push(n0);
                
                bool intersects_previous = false;
                vector<NodeID> sep({n0});
                while(!nq.empty() && !intersects_previous) {
                    NodeID n = nq.front();
                    nq.pop();
                    for(auto nn: g.neighbors(n)) {
                        if(separator_idx[nn]==-1) {
                            int T_in = bfs.T_in[nn];
                            int T_out = bfs.T_out[nn];
                            if(T_in <= T0 && T0 < T_out) {
                                separator_idx[nn] = new_sep_idx;
                                sep.push_back(nn);
                                nq.push(nn);
                            }
                        }
                        else if (separator_idx[nn] != new_sep_idx) {
                            intersects_previous = true;
                            break;
                        }
                    }
                }
                if(intersects_previous) {
                    for(auto n: sep)
                        separator_idx[n] = -1;
                }
                else
                    separators.push_back(sep);
            }
        
        NodeSetVec nsv_for_skel;
        for(const auto& nv: separators)
            if (nv.size() > 10) {
                NodeSet ns = NodeSet(begin(nv),end(nv));
                double c = vertex_separator_curvature(g, ns, bfs.T_out);
                nsv_for_skel.push_back(make_pair(1.0/(1e-5+c), ns));
            }
        
        return nsv_for_skel;
    //    return skeletal_reweighting(g,nsv_for_skel);
    }


    NodeSetVec front_separators(AMGraph3D& g, const vector<AttribVecDouble>& dvv)
    {
        auto process_dist = [](AMGraph3D& g, const AttribVecDouble& dist, int shift) -> NodeSetVec
        {
            auto node_set_vec = separating_node_sets(g, dist, shift);
            return node_set_vec;
        };
        
        size_t N = dvv.size();
        NodeSetVec node_set_vec_global;
        vector<future<NodeSetVec>> nsvfutures(N);
        
        for(int i=0;i<N;++i)
        nsvfutures[i] = async(launch::async, process_dist, ref(g), dvv[i], 0);
        
        for(int i=0;i<N;++i) {
            NodeSetVec nsv =nsvfutures[i].get();
            node_set_vec_global.insert(end(node_set_vec_global), begin(nsv), end(nsv));
        }
        
        greedy_weighted_packing(g, node_set_vec_global, true);
        color_graph_node_sets(g,node_set_vec_global);
        return node_set_vec_global;
    }

    int find_component(const AMGraph3D& g, NodeID n, const vector<NodeSetUnordered>& front_components) {
        int component = -1;
        for(auto m: g.neighbors(n))
            for(int i=0;i<front_components.size();++i)
        if(front_components[i].count(m)) {
            if(component == -1)
                component = i;
            else if (component != i) {
                component = -1;
                return component;
            }
        }
        return component;
    };


    template<typename T>
    void smooth_attribute(const AMGraph3D& g, AttribVec<NodeID, T>& attrib, const NodeSetUnordered& node_set,
                          int N_iter = 1, const AttribVec<NodeID, Vec3d>* _pos = 0) {
        double delta = 0.5;
        const AttribVec<NodeID, Vec3d>& pos = (_pos == 0) ? g.pos : *_pos;
        auto attrib_new = attrib;
        for(int iter=0;iter<N_iter;++iter) {
            for(auto n: node_set) {
                auto N = g.neighbors(n);
                attrib_new[n] = T(0);
                double w_sum = 0.0;
                for(auto m: N) {
                    double w = 1.0/(1e-30+length(pos[m]-pos[n]));
                    attrib_new[n] += w*attrib[m];
                    w_sum += w;
                }
                attrib_new[n] = ((1.0-delta)*attrib[n] + delta*attrib_new[n]/w_sum);
            }
            swap(attrib_new,attrib);
        }
    }


    void node_set_thinning(const AMGraph3D& g, NodeSetUnordered& separator,
                           vector<NodeSetUnordered>& front_components,
                           const AttribVecDouble& priority) {
        using DN_pair = pair<double, NodeID>;
        priority_queue<DN_pair> DNQ;
        for(auto n: separator)
            DNQ.push(make_pair(priority[n],n));
        
        bool did_work = false;
        do {
            did_work = false;
            priority_queue<DN_pair> DNQ_new;
            while (!DNQ.empty()) {
                auto dnp = DNQ.top();
                auto n = dnp.second;
                DNQ.pop();
                int component = find_component(g,n,front_components);
                if(component != -1) {
                    separator.erase(n);
                    front_components[component].insert(n);
                    did_work = true;
                }
                else DNQ_new.push(dnp);
            }
            swap(DNQ_new,DNQ);
        }
        while(did_work);
    }

 

    void optimize_separator(const AMGraph3D& g, NodeSetUnordered& separator,
                          vector<NodeSetUnordered>& front_components) {
        if(separator.size() > 0) {
            NodeID n0 = *begin(separator);
            auto nbors = neighbors(g, separator);
            separator.insert(begin(nbors), end(nbors));
            front_components = connected_components(g, neighbors(g, separator));

            BreadthFirstSearch bfs(g);
            for(auto n: g.node_ids())
                bfs.mask[n] = 0;
            for(auto n: separator)
                bfs.mask[n] = 1;

            bfs.add_init_node(n0);
            while(bfs.Dijkstra_step());
            
            node_set_thinning(g, separator, front_components, bfs.dist);
        }
    }


 
    void shrink_separator(const AMGraph3D& g,
                          NodeSetUnordered& separator,
                          vector<NodeSetUnordered>& front_components,
                          const Vec3d& sphere_centre, int opt_steps) {
        // Next, we thin out the separator until it becomes minimal (i.e. removing one more node
        // would make it cease to be a separator. We remove nodes iteratively and always remove the
        // last visited nodes first.
        const auto separator_orig = separator;
        const auto front_components_orig = front_components;
        
        auto smpos = g.pos;
        AttribVec<NodeID, double> center_dist;
        for(auto n: separator_orig)
            center_dist[n] = sqr_length(smpos[n]-sphere_centre);
        smooth_attribute(g, smpos, separator, sqrt(separator.size()));
        node_set_thinning(g, separator, front_components, center_dist);
        
        for(int iter=0;iter<opt_steps;++iter)
            optimize_separator(g, separator, front_components);
    }



    /** For a given graph, g,  and a given node n0, we compute a local separator.
     The algorithm proceeds in a way similar to Dijkstra, finding a set of nodes separator such that there is anoter set of nodes, front,
     connected to separator via edges and front consists of two connected components.
     thick_front indicates whether we want to add a layer of nodes to the front before checking the number of connected  components.
     persistence is how many iterations the front must have two connected components before we consider the interior
     a local separator.
     The final node set returned is then thinned to the minimal separator.
     */
    pair<double,NodeSet> local_separator(AMGraph3D& g, NodeID n0, double quality_noise_level, int optimization_steps) {
        
        auto front_size_ratio = [](const vector<NodeSetUnordered>& fc) {
            if(fc.size()<2)
                return  0.0;
            auto comp_sizes = minmax({fc[0].size(),fc[1].size()});
            return double(comp_sizes.first)/comp_sizes.second;
        };
        
        // Create the separator node set and the temporary node set (used during computation)
        // The tmp sets are needed because of persistence. Whenever we have had two connected components
        // in front for a number of iterations = persistence, we go back to the original separator.
        NodeSetUnordered Sigma({n0});

        // Create the front node set. Note that a leaf node is a separator by definition,
        // so if there is only one neighbor, we are done here.
        auto N = g.neighbors(n0);
        if(N.size()==0)
            return make_pair(0.0,NodeSet());
        if(N.size()==1)
            return make_pair(1.0, NodeSet({n0}));
        NodeSetUnordered F(begin(N), end(N));

        // We will need node sets for the connected components of the front.
        vector C_F = connected_components(g, F);
        
        // Create the initial sphere which is of radius zero centered at the input node.
        Vec3d centre = g.pos[n0];
        double radius = 0.0;
        
        NodeID last_n = AMGraph3D::InvalidNodeID; // Very last node added to separator.
        // Now, proceed by expanding a sphere
        while (C_F.size()==1 || front_size_ratio(C_F) < quality_noise_level)
        {
            // Find the node in front closest to the center
            const NodeID n = *min_element(begin(F), end(F), [&](NodeID a, NodeID b) {
                return sqr_length(g.pos[a]-centre)<sqr_length(g.pos[b]-centre);
            });
            
            // Update the sphere centre and radius to contain the new point.
            const Vec3d p_n = g.pos[n];
            double l = length(centre-p_n);
            if (l>radius) {
                radius = 0.5 * (radius + l);
                centre = p_n + radius * (centre-p_n)/(1e-30+length(centre-p_n));
            }

            // Now, remove n from F and put it in Sigma.
            // Add n's neighbours (not in Sigma) to F.
            last_n = n;
            F.erase(n);
            Sigma.insert(n);
            for(auto m: g.neighbors(n))
                if(Sigma.count(m)==0)
                    F.insert(m);
            
            // If the front is empty, we must have included an entire
            // connected component in "separator". Bail!
            if(F.size() == 0)
                return make_pair(0.0, NodeSet());
            
            C_F = connected_components(g, F);
        }
        double quality = front_size_ratio(C_F);
//        shrink_separator(g, Sigma, C_F, centre, optimization_steps);
//
//        // We have to check if the local separator is in fact split into two
//        // components. If so, get rid of it.
//        if(connected_components(g, Sigma).size()>1)
//            return make_pair(0.0, NodeSet());
        
        return make_pair(quality, order(Sigma));
    }

    using hrc = chrono::high_resolution_clock;

    NodeSetVec local_separators(AMGraph3D& g, bool sampling, double quality_noise_level, int optimization_steps) {
        
        // Because we are greedy: all cores belong to this task!
        const int CORES = thread::hardware_concurrency();
        
        // touched will help us keep track of how many separators use a given node.
        Util::AttribVec<NodeID, int> touched(g.no_nodes(), 0);
        
        // Create a random order vector of nodes.
        vector<NodeID> node_id_vec;
        for(auto n: g.node_ids())
            node_id_vec.push_back(n);
        srand(1);
        shuffle(begin(node_id_vec), end(node_id_vec), default_random_engine(rand()));
        
        auto t1 = hrc::now();
        // Each core will have its own vector of NodeSets in which to store
        // separators.
        vector<NodeSetVec> nsvv(CORES);
        int cnt = 0;
        auto create_separators = [&](int core) {
            auto& nsv = nsvv[core];
            for(auto n: node_id_vec) {
                double probability = 1.0/int_pow(2.0, touched[n]);
                if (n%CORES==core && (!sampling || rand() <= probability*RAND_MAX)) {
                    cnt += 1;
                    auto ns = local_separator(g, n, quality_noise_level, optimization_steps);
                    if(ns.second.size()>0) {
                        nsv.push_back(ns);
                        for(auto m: ns.second)
                            touched[m] += 1;
                    }
                }
            }
        };

        vector<thread> threads(CORES);
        for(int i=0;i<CORES;++i)
            threads[i] = thread(create_separators, i);

        for(int i=0;i<CORES;++i)
            threads[i].join();

        auto t2 = hrc::now();

        NodeSetVec node_set_vec_global;
        for(const auto& nsv: nsvv)
            for(const auto& ns: nsv)
                node_set_vec_global.push_back(ns);

        auto sep_bef = node_set_vec_global.size();
        greedy_weighted_packing(g, node_set_vec_global, true);
        auto sep_aft = node_set_vec_global.size();
        auto t3 = hrc::now();
        
        cout << "Computed " << cnt << " separators" << endl;
        cout << "Found " << sep_bef << " separators" << endl;
        cout << "Packed " << sep_aft << " separators" << endl;
        cout << "Finding separators: " << (t2-t1).count()*1e-9 << endl;
        cout << "Packing separators: " << (t3-t2).count()*1e-9 << endl;

        // Color the node sets selected by packing, so we can get a sense of the
        // selection.
        color_graph_node_sets(g, node_set_vec_global);

        return node_set_vec_global;
    }
}
