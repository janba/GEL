#ifndef dyncon_hpp
#define dyncon_hpp

#include <iostream>
#include <list>
#include <stack>
#include <unordered_map>
#include <set>
#include <map>
#include <queue>

namespace Geometry {
    template<typename T>
    class DynCon {
    private:
        template<typename WT>
        struct UndirectedEdge {
            // like a pair but (first, second) is the same as (second, first).
            WT first, second;
            // Location in adjacency list.
            typename std::list<WT>::iterator adjListFirst, adjListSecond;

            UndirectedEdge(const WT a, const WT b) {
                first = std::min(a, b);
                second = std::max(a, b);
            }

            bool operator==(const UndirectedEdge &e) const {
                return (first == e.first && second == e.second);
            }
        };

        template<typename WT>
        struct UndirectedEdgeHash {
            size_t operator()(const UndirectedEdge<WT> &e) const {
                uint32_t v = e.first, w = e.second;
                return std::hash<uint64_t>()((uint64_t) v << 32 | w);
            }
        };

        template<typename WT>
        struct VertexPairHash {
            size_t operator()(const std::pair<WT, WT> &e) const {
                uint32_t v = e.first, w = e.second;
                return std::hash<uint64_t>()((uint64_t) v << 32 | w);
            }
        };

        template<typename WT>
        struct SplayT {
            WT u, v, size; // endpoints and size of subtree
            bool adjT, adjNT, marked; // subtree has tree/nontree to consider, node itself should be considered
            size_t l, r, p; // left,right,parent

            SplayT(WT x, WT y) {
                u = x;
                v = y;
                size = 1;
                adjNT = false;
                marked = x < y;
                adjT = marked;
                l = r = p = -1;
            }
        };

        template<typename WT>
        class st_wrapper{
        private:
            std::vector<SplayT<WT>> v;

        public:
            size_t p(size_t t){return v[t].p;}
            size_t lc(size_t t){return v[t].l;}
            size_t rc(size_t t){return v[t].r;}

            bool has_p(size_t t){return v[t].p != -1;}
            bool has_l(size_t t){return v[t].l != -1;}
            bool has_r(size_t t){return v[t].r != -1;}

            size_t add(WT p, WT q){
                v.emplace_back(SplayT<WT>(p,q));
                return v.size()-1;
            }

            size_t find_representative(size_t t) {
                while (has_p(t)) {
                    //std::cout << "("<<t->u<<","<<t->v << ")-";
                    t = p(t);
                }
                return t;
            }

            int count(size_t t) {
                return t!=-1 ? v[t].size : 0;
            }

            bool hasAdjNT(size_t t) {
                return t != -1 && v[t].adjNT;
            }

            bool hasAdjT(size_t t) {
                return t != -1 && v[t].adjT;
            }

            void update_data(size_t t) {
                if (t!=-1) {
                    v[t].size = 1 + count(lc(t)) + count(rc(t));
                    v[t].adjNT = hasAdjNT(lc(t)) || hasAdjNT(rc(t)) || (v[t].u == v[t].v && v[t].marked);
                    v[t].adjT = hasAdjT(lc(t)) || hasAdjT(rc(t)) || (v[t].u < v[t].v && v[t].marked);
                }
            }

            void propagate_data(size_t t) {
                bool pastNT;
                bool pastT;
                while (t!=-1) {
                    pastNT = v[t].adjNT;
                    pastT = v[t].adjT;
                    update_data(t);
                    if ((pastNT != v[t].adjNT) || (pastT != v[t].adjT)) t = p(t);
                    else break;
                }
            }

            void rotateR(size_t t) {
                size_t l = lc(t);
                v[t].l = rc(l);
                if (has_r(l)) v[rc(l)].p = t;
                v[l].p = p(t);
                if (has_p(t)) {
                    if (lc(p(t)) == t) v[p(t)].l = l;
                    else v[p(t)].r = l;
                }
                v[l].r = t;
                v[t].p = l;
                update_data(t);
                update_data(l);
            }

            void rotateL(size_t t) {
                size_t r = rc(t);
                v[t].r = lc(r);
                if (has_l(r)) v[lc(r)].p = t;
                v[r].p = p(t);
                if (has_p(t)) {
                    if (lc(p(t)) == t) v[p(t)].l = r;
                    else v[p(t)].r = r;
                }
                v[r].l = t;
                v[t].p = r;
                update_data(t);
                update_data(r);
            }

            void splay(size_t t) {
                while (has_p(t)) {
                    if (!has_p(p(t))) {
                        if (lc(p(t)) == t) rotateR(p(t));
                        else rotateL(p(t));
                    } else if (lc(p(t)) == t && lc(p(p(t))) == p(t)) {
                        rotateR(p(p(t)));
                        rotateR(p(t));
                    } else if (rc(p(t)) == t && rc(p(p(t))) == p(t)) {
                        rotateL(p(p(t)));
                        rotateL(p(t));
                    } else if (rc(p(t)) == t && lc(p(p(t))) == p(t)) {
                        rotateL(p(t));
                        rotateR(p(t));
                    } else {
                        rotateR(p(t));
                        rotateL(p(t));
                    }
                }
            }

            void join(size_t& t, size_t l, size_t r) {
                if(r==-1) {t = l; return;}
                if(l==-1) {t = r; return;}

                size_t max = l;
                while (has_r(max)) max = rc(max);
                splay(max);
                v[max].r = r;
                v[r].p = max;
                update_data(max);
                t = max;
            }

            void split(size_t t, size_t &l, size_t &r) {
                splay(t);
                l = lc(t);
                v[t].l = -1;
                r = t;
                if (l!=-1) {
                    v[l].p = -1;
                    update_data(l);
                }
                update_data(r);
            }

            void remove_first(size_t t) {
                if (has_r(t)) v[rc(t)].p = -1;
                v[t].l = v[t].r = v[t].p = -1;
            }

            void mark(size_t t, bool val){
                v[t].marked = val;
                propagate_data(t);
            }

            size_t reroot(size_t t) {
                size_t A, B;
                split(t, A, B);
                join(B,B, A);
                return B;
            }

            void verify_children(size_t t) {
                if (t==-1) return;
                if (has_p(t) && !(lc(p(t)) == t || rc(p(t)) == t))
                    std::cout << "ERROR P " << "(" << v[t].u << "," << v[t].v << ")" << "," << "(" << v[p(t)].u << "," << v[p(t)].v
                              << ")" << std::endl;
                if (has_l(t) && (p(lc(t)) != t))
                    std::cout << "ERROR L " << "(" << v[t].u << "," << v[t].v << ")" << "," << "(" << v[lc(t)].u << "," << v[lc(t)].v
                              << ")" << std::endl;
                if (has_r(t) && (p(rc(t)) != t))
                    std::cout << "ERROR R " << "(" << v[t].u << "," << v[t].v << ")" << "," << "(" << v[rc(t)].u << "," << v[rc(t)].v
                              << ")" << std::endl;
                if (v[t].size != 1 + count(lc(t)) + count(rc(t))) std::cout << "ERROR SIZE" << std::endl;
                verify_children(lc(t));
                verify_children(rc(t));
            }

            SplayT<WT> access(size_t t){
                return v[t];
            }
        };

        template<typename WT>
        class EulerTourForest {
        private:
            std::map<std::pair<WT, WT>, size_t> v_map;
            st_wrapper<WT>& stw;
            //std::unordered_map<std::pair<WT, WT>, SplayT<WT> *, VertexPairHash<WT>> v_map;

        public:
            EulerTourForest(DynCon::st_wrapper<WT>& _stw) : stw(_stw){}

            size_t add(WT u, WT v){
                size_t t = stw.add(u,v);
                v_map[std::pair<WT, WT>(u, v)] = t;
                return t;
            }

            size_t find_tree(WT u, WT v) {
                auto search = v_map.find(std::pair<WT, WT>(u, v));
                if (search == v_map.end()) return -1;
                else return search->second;
            }

            size_t find_tree(WT x) {
                return find_tree(x, x);
            }

            size_t link(WT u, WT v) {
                size_t t_u = find_tree(u);
                size_t t_v = find_tree(v);

                if (t_u == -1) t_u = add(u,u);
                if (t_v == -1) t_v = add(v,v);

                stw.splay(t_u);
                stw.splay(t_v);
                if (stw.p(t_u) == t_v || (stw.has_p(t_u) && stw.p(stw.p(u)) == t_v)) return t_v;

                t_u = stw.reroot(t_u);
                t_v = stw.reroot(t_v);

                auto forward = add(u,v);
                auto retreat = add(v, u);

                stw.join(t_u, t_u, forward);
                stw.join(t_v, t_v, retreat);
                stw.join(t_u, t_u, t_v);

                return t_u;
            }

            // Cuts the edge between u and v. A and B will be the roots of the new trees after the cut (not the first element in ETT sequence).
            void cut(WT u, WT v, size_t &A, size_t &B) {
                auto s_advance = v_map.find(std::pair<WT, WT>(u, v));
                auto s_retreat = v_map.find(std::pair<WT, WT>(v, u));
                if (s_advance == v_map.end() || s_retreat == v_map.end()) {
                    A = B = -1;
                    return;
                }
                size_t J, K, L;
                size_t t_advance = s_advance->second;
                size_t t_retreat = s_retreat->second;

                stw.split(t_advance, J, L);
                L = stw.rc(L);
                stw.remove_first(t_advance);
                v_map.erase(s_advance);

                if (stw.find_representative(t_retreat) == J) {
                    stw.split(t_retreat, J, K);
                    K = stw.rc(K);
                } else {
                    stw.split(t_retreat, K, L);
                    L = stw.rc(L);
                }

                stw.remove_first(t_retreat);
                v_map.erase(s_retreat);

                A = K;
                stw.join(B, J, L);
            }

            bool is_connected(WT u, WT v) {
                size_t t_u = find_tree(u);
                size_t t_v = find_tree(v);

                if (!(t_u!=-1 && t_v!=-1)) return false;

                stw.splay(t_u);
                stw.splay(t_v);

                return (stw.p(t_u) == t_v || (stw.has_p(t_u) && stw.p(stw.p(t_u)) == t_v));
            }

            bool edge_exists(WT u, WT v) {
                return find_tree(u, v) != -1;
            }

            void mark(WT x, WT y, bool val) {
                auto t = find_tree(x, y);
                if (t == -1) {t = add(x,y);}

                stw.mark(t,val);
            }

            void mark(WT x, bool val) { mark(x, x, val); }
        };

        EulerTourForest<T> *forest;

        // An incident list for each vertex. Uses adjacency lists since we then only have to store one int.
        std::map<T, std::list<T> *> adjacencyLists;

        std::unordered_map<UndirectedEdge<T>, UndirectedEdge<T> *, UndirectedEdgeHash<T>> edgeSet;

        // Adds a non-tree edge with the given specifications
        void addNonTree(T v, T w, UndirectedEdge<T> *&e) {
            if (v > w) std::swap(v, w);

            // Mark that the corresponding nodes in forest has adjacent non-tree edges
            forest->mark(v, true);
            forest->mark(w, true);

            // Add e to adjacency lists.
            if (adjacencyLists.count(v) == 0) {
                adjacencyLists[v] = new std::list<T>();
            }
            if (adjacencyLists.count(w) == 0) {
                adjacencyLists[w] = new std::list<T>();
            }
            auto listV = adjacencyLists.find(v)->second;
            auto listW = adjacencyLists.find(w)->second;

            listV->push_front(w);
            listW->push_front(v);

            // Store iterators to positions in adjacency lists for O(1) removal
            e->adjListFirst = listV->begin();
            e->adjListSecond = listW->begin();
        }

        void disconnect_nontree(UndirectedEdge<T> *&e) {
            int v = *e->adjListSecond;
            int w = *e->adjListFirst;
            auto listV = adjacencyLists.find(v)->second;
            auto listW = adjacencyLists.find(w)->second;
            listV->erase(e->adjListFirst);
            listW->erase(e->adjListSecond);
            if (listV->empty()) forest->mark(v, false);
            if (listW->empty()) forest->mark(w, false);
        }

        st_wrapper<T> stw;

    public:
        DynCon() {
            forest = new EulerTourForest<T>(stw);
            adjacencyLists = {};
            edgeSet = {};
        }

        ~DynCon() {
            delete forest;

            for (auto l: adjacencyLists) {
                delete l.second;
            }

            for (auto e: edgeSet) {
                delete e.second;
            }
        }

        // Insert the edge going from v to w
        void insert(T v, T w) {
            // Swap v and w so that an edge (v, w) is the same as (w, v).
            if (v > w) std::swap(v, w);

            // Edge already exists - early return
            if (edgeSet.find(UndirectedEdge<T>(v, w)) != edgeSet.end()) {
                return;
            }

            auto e = new UndirectedEdge<T>(v, w);
            edgeSet[*e] = e;

            // if v & w disconnected in F_0 -> insert as tree edge in F_0.
            if (!is_connected(v, w)) {
                forest->link(v, w);
            } else { // Add as non-tree otherwise
                addNonTree(v, w, e);
            }
        }

        // Returns false if no replacement edge found
        bool remove(T v, T w) {
            if (v > w) std::swap(v, w);

            // Early return if edge doesn't exist
            auto search = edgeSet.find(UndirectedEdge<T>(v, w));
            if (search == edgeSet.end()) return true;

            auto e = search->second;

            // Remove from edgeset
            edgeSet.erase(*e);

            // If edge is non-tree, removal cannot break connectivity
            if (!forest->edge_exists(v, w)) {
                // But we must update adjacency lists
                disconnect_nontree(e);
                delete e;
                return true;
            }

            UndirectedEdge<T> *replacement = nullptr;

            size_t V, W;

            forest->cut(v, w, V, W);

            // Ensure that |V| < |W|. We only need to use V.
            if (stw.access(V).size > stw.access(W).size) std::swap(V, W);

            // If there are no adjacent non-tree edges in V no need for searching
            if (!stw.hasAdjNT(V)) {
                delete e;
                return false;
            }

            // Iterative level order traversal of V and non-tree edges
            std::queue<size_t> queue;
            size_t current = V;
            queue.push(current);

            while (!replacement && !queue.empty()) { // Process
                current = queue.front();
                queue.pop();

                if (stw.has_l(current) && stw.hasAdjNT(stw.lc(current))) queue.push(stw.lc(current));
                if (stw.has_r(current) && stw.hasAdjNT(stw.rc(current))) queue.push(stw.rc(current));

                if (stw.access(current).marked &&
                    stw.access(current).u == stw.access(current).v) { // Node represents vertex with adjacent non-tree edges
                    auto list = adjacencyLists.find(stw.access(current).u)->second;
                    auto iter = list->begin();
                    while (iter != list->end()) { // Destructive iteration of adjacent non-tree edges
                        T candidate = *iter;
                        iter++;
                        auto c_e = edgeSet.find(UndirectedEdge<T>(stw.access(current).u, candidate))->second;
                        if (stw.find_representative(forest->find_tree(candidate)) != V) { // Non-tree edge goes to W
                            disconnect_nontree(c_e);
                            replacement = c_e;
                            forest->link(c_e->first, c_e->second);
                            break;
                        }
                    }
                }
            }
            delete e;

            if (replacement) {
                return true;
            }
            return false;
        }

        // Returns true if there exists a path going between v and w.
        bool is_connected(T v, T w) {
            // Check connected in F_0.
            return forest->is_connected(v, w);
        }

        // Returns representative of component of v
        T get_representative(T v) {
            size_t rep = forest->find_tree(v);
            if (rep == -1) return -1;
            rep = stw.find_representative(rep);
            return stw.access(rep).u;
        }

        T get_size(T v) {
            // Returns size of component of v
            // Structure actually stores |V|+2*|E| nodes
            // Means (size+2)/3 vertices are in structure
            size_t rep = forest->find_tree(v);
            if (rep == -1) return 0;
            return (stw.access(stw.find_representative(rep)).size + 2) / 3;
        }

        void print_tree(T v) {
            in_order(find_representative(forest->find_tree(v)));
        }
    };
}

#endif //dyncon_hpp
