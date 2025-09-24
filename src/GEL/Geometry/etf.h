//
// Created by etoga on 11/6/23.
//

#ifndef GEOMETRY_ETF_H
#define GEOMETRY_ETF_H

#include <vector>

namespace Geometry
{
using uint = std::uint32_t;

namespace detail
{
    constexpr
    uint pcg_hash(const uint input)
    {
        const uint state = input * 747796405u + 2891336453u;
        const uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
        return (word >> 22u) ^ word;
    }
}


/// Treap - data structure for a randomized binary search tree
/// @sa ETF
struct Treap {
    uint priority;
    uint size = 1u;
    uint left_child = 0;
    uint right_child = 0;
    uint parent = 0;

    Treap(uint x, uint y): priority(detail::pcg_hash(x << sizeof(uint) / 2 ^ y)) {}
    Treap(): priority(0), size(0) {}
};

/// Euler Tour Forest
struct ETF {
private:
    uint root = 0;
    std::vector<Treap> arena = {Treap()};

    [[nodiscard]]
    const uint& parent(const uint i) const
    {
        return arena[i].parent;
    }

    uint& parent(const uint i)
    {
        return arena[i].parent;
    }

    [[nodiscard]]
    const uint& left_child(const uint i) const
    {
        return arena[i].left_child;
    }

    uint& left_child(const uint i)
    {
        return arena[i].left_child;
    }

    [[nodiscard]]
    const uint& right_child(const uint i) const
    {
        return arena[i].right_child;
    }

    uint& right_child(const uint i)
    {
        return arena[i].right_child;
    }

    [[nodiscard]]
    bool is_lc(const uint i) const
    {
        return (parent(i) && left_child(parent(i)) == i);
    }

    [[nodiscard]]
    uint priority(const uint i) const
    {
        return arena[i].priority;
    }

    [[nodiscard]]
    uint count(const uint i) const
    {
        return arena[i].size;
    }

    void update_count(uint i)
    {
        if (i) arena[i].size = 1 + count(left_child(i)) + count(right_child(i));
    }

    void join(uint& t, const uint l, const uint r)
    {
        if (!l || !r) t = l ? l : r;
        else if (priority(l) > priority(r)) {
            join(right_child(l), right_child(l), r);
            parent(right_child(l)) = l;
            t = l;
        } else {
            join(left_child(r), l, left_child(r));
            parent(left_child(r)) = r;
            t = r;
        }
        update_count(t);
    }

    // Split s.t. l < t <= r
    void split(uint t, uint& l, uint& r)
    {
        l = left_child(t);
        if (l) parent(l) = 0;
        r = t;
        arena[t].size -= count(l);
        left_child(t) = 0;
        while (parent(t)) {
            if (is_lc(t)) {
                if (t == l) {
                    left_child(parent(t)) = r;
                    if (r) parent(r) = parent(t);
                    r = parent(t);
                    parent(t) = 0;
                } else {
                    r = parent(t);
                }
                t = r;
            } else {
                if (t == r) {
                    right_child(parent(t)) = l;
                    if (l) parent(l) = parent(t);
                    l = parent(t);
                    parent(t) = 0;
                } else {
                    l = parent(t);
                }
                t = l;
            }
            update_count(t);
        }
        update_count(t);
    }

public:
    void reserve(const uint size)
    {
        arena.reserve(size);
    }

    [[nodiscard]]
    uint representative(uint t) const
    {
        while (parent(t)) {
            t = parent(t);
        }
        return t;
    }

    uint accumulate()
    {
        uint size = arena.size();
        arena.emplace_back(size, size);
        if (size > 1) {
            uint _;
            join(_, root, size);
            root = representative(root);
        } else {
            root = 1;
        }
        return size;
    }

    /// Takes indices in the ETF arena, returns indices of the new edges in the ETF arena
    std::pair<uint, uint> insert(const uint u, const uint v)
    {
        uint advance = arena.size();
        arena.emplace_back(u, v);
        uint retreat = arena.size();
        arena.emplace_back(v, u);

        // Do the split
        uint A, B, J = 0, K = 0, L = 0;
        split(u, A, B);
        if (representative(v) == A) {
            // v left of u
            split(v, J, K);
            L = B;
            A = advance;
            B = retreat;
        } else {
            // u left of v
            split(v, K, L);
            J = A;
            A = retreat;
            B = advance;
        }
        join(K, K, A);
        join(J, J, B);
        join(J, J, L);
        return {advance, retreat};
    }

    [[nodiscard]]
    bool connected(const uint u, const uint v) const
    {
        return representative(u) == representative(v);
    }

    bool connected(const uint u, const uint v, uint& out_tree, uint& out_to_tree) const
    {
        out_tree = representative(u);
        out_to_tree = representative(v);
        return out_tree == out_to_tree;
    }
};
}

#endif // GEOMETRY_ETF_H
