//
// Created by etoga on 11/6/23.
//

#ifndef C_ETF_H
#define C_ETF_H

#include <vector>
#include <map>


using uint = unsigned int;

inline
uint pcg_hash(uint input){
    uint state = input * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

struct Treap{
    uint prio, size = 1u;
    uint l=0, r=0, p=0;

    Treap(uint x, uint y): prio(pcg_hash(x<<sizeof(uint)/2 ^ y)){};
    Treap(): prio(0),size(0){};
};

struct ETF{
private:
    uint root = 0;
    std::vector<Treap> arena = {Treap()};

    inline
    uint& par(uint i){
        return arena[i].p;
    }

    inline
    uint& lc(uint i){
        return arena[i].l;
    }

    inline
    uint& rc(uint i){
        return arena[i].r;
    }

    inline
    bool is_lc(uint i){
        return (par(i) && lc(par(i)) == i);
    }

    inline
    uint prio(uint i){
        return arena[i].prio;
    }

    inline
    uint count(uint i){
        return arena[i].size;
    }

    inline
    void update_count(uint i){
        if(i) arena[i].size = 1 + count(lc(i)) + count(rc(i));
    }

    void join(uint& t, uint l, uint r){
        if(!l || !r) t = l ? l:r;
        else if (prio(l) > prio(r)){
            join(rc(l), rc(l), r);
            par(rc(l)) = l;
            t = l;
        } else {
            join(lc(r), l, lc(r));
            par(lc(r)) = r;
            t = r;
        }
        update_count(t);
    }

    // Split s.t. l < t <= r
    void split(uint t, uint& l, uint& r){
        l = lc(t);
        if(l) par(l) = 0;
        r = t;
        arena[t].size -= count(l);
        lc(t) = 0;
        while(par(t)){
            if(is_lc(t)){
                if(t == l){
                    lc(par(t)) = r;
                    if(r) par(r) = par(t);
                    r = par(t);
                    par(t) = 0;
                } else {
                    r = par(t);
                }
                t = r;
            } else {
                if(t == r){
                    rc(par(t)) = l;
                    if(l) par(l) = par(t);
                    l = par(t);
                    par(t) = 0;
                } else {
                    l = par(t);
                }
                t = l;
            }
            update_count(t);
        }
        update_count(t);
    }



public:
    void reserve(uint size){
        arena.reserve(size);
    }

    uint representative(uint t) {
        while(par(t)){
            t = par(t);
        }
        return t;
    }

    uint accumulate(){
        uint size = arena.size();
        arena.emplace_back(size,size);
        if(size>1) {
            uint _;
            join(_, root, size);
            root = representative(root);
        } else {
            root = 1;
        }
        return size;
    }

    // Takes indices in the ETF arena
    // Returns indices of the new edges in the ETF arena
    std::pair<uint,uint> insert(const uint u, const uint v){
        uint advance = arena.size();
        arena.emplace_back(u,v);
        uint retreat = arena.size();
        arena.emplace_back(v,u);

        // Do the split
        uint A,B,J=0,K=0,L=0;
        split(u,A,B);
        if(representative(v)==A){ // v left of u
            split(v,J,K);
            L = B;
            A = advance;
            B = retreat;
        } else { // u left of v
            split(v,K,L);
            J = A;
            A = retreat;
            B = advance;
        }
        join(K,K,A);
        join(J,J,B);
        join(J,J,L);
        return {advance,retreat};
    }

    bool connected(const uint u, const uint v){
        return representative(u)==representative(v);
    }

    bool connected(const uint u, const uint v, uint& out_tree, uint& out_to_tree) {
        out_tree = representative(u);
        out_to_tree = representative(v);
        return out_tree == out_to_tree;
    }
};

#endif //C_ETF_H
