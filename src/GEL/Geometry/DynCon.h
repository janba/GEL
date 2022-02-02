#ifndef dyncon_hpp
#define dyncon_hpp

#include "DynCon.h"
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
        struct SplayT {
            WT u, v, size; // endpoints and size of subtree
            bool adjT, adjNT, marked; // subtree has tree/nontree to consider, node itself should be considered
            SplayT<WT> *l, *r, *p; // left,right,parent

            SplayT(WT x, WT y) {
                u = x;
                v = y;
                size = 1;
                adjNT = false;
                marked = x < y;
                adjT = marked;
                l = r = p = nullptr;
            }
        };

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
        static void in_order(SplayT<WT> *t) {
            if (!t) return;
            in_order(t->l);
            std::cout << "(" << t->u << "," << t->v << ")";
            in_order(t->r);
        }

        template<typename WT>
        static void post_order(SplayT<WT> *t) {
            if (!t) return;
            post_order(t->l);
            post_order(t->r);
            std::cout << "(" << t->u << "," << t->v << ")";
        }

        template<typename WT>
        static SplayT<WT> *find_representative(SplayT<WT> *t) {
            while (t->p != nullptr) {
                //std::cout << "("<<t->u<<","<<t->v << ")-";
                t = t->p;
            }
            return t;
        }

        template<typename WT>
        static int count(SplayT<WT> *t) {
            return t ? t->size : 0;
        }

        template<typename WT>
        static bool hasAdjNT(SplayT<WT> *t) {
            return t != nullptr && t->adjNT;
        }

        template<typename WT>
        static bool hasAdjT(SplayT<WT> *t) {
            return t != nullptr && t->adjT;
        }

        template<typename WT>
        static void update_data(SplayT<WT> *t) {
            if (t) {
                t->size = 1 + count(t->l) + count(t->r);
                t->adjNT = hasAdjNT(t->l) || hasAdjNT(t->r) || (t->u == t->v && t->marked);
                t->adjT = hasAdjT(t->l) || hasAdjT(t->r) || (t->u < t->v && t->marked);
            }
        }

        template<typename WT>
        static void propagate_data(SplayT<WT> *t) {
            bool pastNT;
            bool pastT;
            while (t) {
                pastNT = t->adjNT;
                pastT = t->adjT;
                update_data(t);
                if ((pastNT != t->adjNT) || (pastT != t->adjT)) t = t->p;
                else break;
            }
        }

        template<typename WT>
        static void rotateR(SplayT<WT> *t) {
            SplayT<WT> *l = t->l;
            t->l = l->r;
            if (l->r) l->r->p = t;
            l->p = t->p;
            if (t->p) {
                if (t->p->l == t) t->p->l = l;
                else t->p->r = l;
            }
            l->r = t;
            t->p = l;
            update_data(t);
            update_data(l);
        }

        template<typename WT>
        static void rotateL(SplayT<WT> *t) {
            SplayT<WT> *r = t->r;
            t->r = r->l;
            if (r->l) r->l->p = t;
            r->p = t->p;
            if (t->p) {
                if (t->p->l == t) t->p->l = r;
                else t->p->r = r;
            }
            r->l = t;
            t->p = r;
            update_data(t);
            update_data(r);
        }

        template<typename WT>
        static void splay(SplayT<WT> *t) {
            while (t->p) {
                if (!t->p->p) {
                    if (t == t->p->l) rotateR(t->p);
                    else rotateL(t->p);
                } else if (t == t->p->l && t->p == t->p->p->l) {
                    rotateR(t->p->p);
                    rotateR(t->p);
                } else if (t == t->p->r && t->p == t->p->p->r) {
                    rotateL(t->p->p);
                    rotateL(t->p);
                } else if (t == t->p->r && t->p == t->p->p->l) {
                    rotateL(t->p);
                    rotateR(t->p);
                } else {
                    rotateR(t->p);
                    rotateL(t->p);
                }
            }
        }

        template<typename WT>
        static void join(SplayT<WT> *&t, SplayT<WT> *l, SplayT<WT> *r) {
            if (!l || !r) {
                t = l ? l : r;
                return;
            }
            SplayT<WT> *max = l;
            while (max->r) max = max->r;
            splay(max);
            max->r = r;
            r->p = max;
            t = max;
            update_data(t);
        }

        template<typename WT>
        static void split(SplayT<WT> *t, SplayT<WT> *&l, SplayT<WT> *&r) {
            splay(t);
            l = t->l;
            t->l = nullptr;
            r = t;
            if (l) {
                l->p = nullptr;
                update_data(l);
            }
            update_data(r);
        }

        template<typename WT>
        static void remove_first(SplayT<WT> *t) {
            if (t->r) t->r->p = nullptr;
            t->l = t->r = t->p = nullptr;
            delete t;
        }

        template<typename WT>
        static SplayT<WT> *reroot(SplayT<WT> *t) {
            SplayT<WT> *A, *B;
            split(t, A, B);
            join(B, B, A);
            return B;
        }

        template<typename WT>
        static void verify_children(SplayT<WT> *t) {
            if (!t) return;
            if (t->p && !(t->p->l == t || t->p->r == t))
                std::cout << "ERROR P " << "(" << t->u << "," << t->v << ")" << "," << "(" << t->p->u << "," << t->p->v
                          << ")" << std::endl;
            if (t->l && (t->l->p != t))
                std::cout << "ERROR L " << "(" << t->u << "," << t->v << ")" << "," << "(" << t->l->u << "," << t->l->v
                          << ")" << std::endl;
            if (t->r && (t->r->p != t))
                std::cout << "ERROR R " << "(" << t->u << "," << t->v << ")" << "," << "(" << t->r->u << "," << t->r->v
                          << ")" << std::endl;
            if (t->size != 1 + count(t->l) + count(t->r)) std::cout << "ERROR SIZE" << std::endl;
            verify_children(t->l);
            verify_children(t->r);
        }

        template<typename WT>
        class EulerTourForest {
        private:
            std::map<std::pair<WT, WT>, SplayT<WT> *> v_map;
            //std::unordered_map<std::pair<WT, WT>, SplayT<WT> *, VertexPairHash<WT>> v_map;

        public:
            EulerTourForest() = default;

            ~EulerTourForest() {
                for (auto t: v_map) {
                    delete t.second;
                }
            }

            SplayT<WT> *find_tree(WT u, WT v) {
                auto search = v_map.find(std::pair<WT, WT>(u, v));
                if (search == v_map.end()) return nullptr;
                else return search->second;
            }

            SplayT<WT> *find_tree(WT x) {
                return find_tree(x, x);
            }

            SplayT<WT> *link(WT u, WT v) {
                SplayT<WT> *t_u = find_tree(u);
                SplayT<WT> *t_v = find_tree(v);

                if (!t_u) {
                    t_u = new SplayT<WT>(u, u);
                    v_map[std::pair<WT, WT>(u, u)] = t_u;
                }
                if (!t_v) {
                    t_v = new SplayT<WT>(v, v);
                    v_map[std::pair<WT, WT>(v, v)] = t_v;
                }

                splay(t_u);
                splay(t_v);
                if (t_u->p == t_v || (t_u->p && t_u->p->p == t_v)) return t_v;

                t_u = reroot(t_u);
                t_v = reroot(t_v);

                auto forward = new SplayT<WT>(u, v);
                v_map[std::pair<WT, WT>(u, v)] = forward;
                auto retreat = new SplayT<WT>(v, u);
                v_map[std::pair<WT, WT>(v, u)] = retreat;

                join(t_u, t_u, forward);
                join(t_v, t_v, retreat);
                join(t_u, t_u, t_v);

                return t_u;
            }

            // Cuts the edge between u and v. A and B will be the roots of the new trees after the cut (not the first element in ETT sequence).
            void cut(WT u, WT v, SplayT<WT> *&A, SplayT<WT> *&B) {
                auto s_advance = v_map.find(std::pair<WT, WT>(u, v));
                auto s_retreat = v_map.find(std::pair<WT, WT>(v, u));
                if (s_advance == v_map.end() || s_retreat == v_map.end()) {
                    A = B = nullptr;
                    return;
                }
                SplayT<WT> *J, *K, *L;
                SplayT<WT> *t_advance = s_advance->second;
                SplayT<WT> *t_retreat = s_retreat->second;

                split(t_advance, J, L);
                L = L->r;
                remove_first(t_advance);
                v_map.erase(s_advance);

                if (find_representative(t_retreat) == J) {
                    split(t_retreat, J, K);
                    K = K->r;
                } else {
                    split(t_retreat, K, L);
                    L = L->r;
                }

                remove_first(t_retreat);
                v_map.erase(s_retreat);

                A = K;
                join(B, J, L);
            }

            bool is_connected(WT u, WT v) {
                SplayT<WT> *t_u = find_tree(u);
                SplayT<WT> *t_v = find_tree(v);

                if (!(t_u && t_v)) return false;

                splay(t_u);
                splay(t_v);

                return (t_u->p == t_v || (t_u->p && t_u->p->p == t_v));
            }

            bool edge_exists(WT u, WT v) {
                return find_tree(u, v) != nullptr;
            }

            void mark(WT x, WT y, bool val) {
                auto t = find_tree(x, y);
                if (!t) {
                    t = new SplayT<WT>(x, y);
                    v_map[std::pair<WT, WT>(x, y)] = t;
                }

                t->marked = val;
                propagate_data(t);
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

    public:
        DynCon() {
            forest = new EulerTourForest<T>;
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

            SplayT<T> *V, *W;

            forest->cut(v, w, V, W);

            // Ensure that |V| < |W|. We only need to use V.
            if (V->size > W->size) std::swap(V, W);

            // If there are no adjacent non-tree edges in V no need for searching
            if (!V->adjNT) {
                delete e;
                return false;
            }

            // Iterative level order traversal of V and non-tree edges
            std::queue<SplayT<T> *> queue;
            SplayT<T> *current = V;
            queue.push(current);

            while (!replacement && !queue.empty()) { // Process
                current = queue.front();
                queue.pop();

                if (current->l && current->l->adjNT) queue.push(current->l);
                if (current->r && current->r->adjNT) queue.push(current->r);

                if (current->marked &&
                    current->u == current->v) { // Node represents vertex with adjacent non-tree edges
                    auto list = adjacencyLists.find(current->u)->second;
                    auto iter = list->begin();
                    while (iter != list->end()) { // Destructive iteration of adjacent non-tree edges
                        T candidate = *iter;
                        iter++;
                        auto c_e = edgeSet.find(UndirectedEdge<T>(current->u, candidate))->second;
                        if (find_representative(forest->find_tree(candidate)) != V) { // Non-tree edge goes to W
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
            SplayT<T> *rep = forest->find_tree(v);
            if (!rep) return -1;
            rep = find_representative(rep);
            return rep->u;
        }

        T get_size(T v) {
            // Returns size of component of v
            // Structure actually stores |V|+2*|E| nodes
            // Means (size+2)/3 vertices are in structure
            SplayT<T> *rep = forest->find_tree(v);
            if (!rep) return 0;
            return (find_representative(rep)->size + 2) / 3;
        }

        void print_tree(T v) {
            in_order(find_representative(forest->find_tree(v)));
        }
    };
}

#endif //dyncon_hpp
