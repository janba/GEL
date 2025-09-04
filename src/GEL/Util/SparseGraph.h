//
// Created by Cem Akarsubasi on 8/29/25.
//

#ifndef GEL_SPARSE_GRAPH_H
#define GEL_SPARSE_GRAPH_H

#include <ranges>
#include <limits>
#include <GEL/Util/AssociativeContainers.h>
#include <GEL/Util/Range.h>
#include <queue>
#include <vector>

namespace Util::detail
{

/// @brief A sparse graph implementation heavily based on the one shown
/// in github.com/ashvardanian/less_slow.cpp (Apache 2.0 License)
///
/// Used mainly for experimentation as this can be significantly
/// more efficient than AMGraph when the graph is actually sparse.
///
/// There is an implicit requirement that the built in data structure
/// allow for returning multiple values from std::equal_range despite
/// this seeminly not making any sense for containers besides multimap
/// and multiset. Apparently this requirement is met by the absl containers
/// and the std containers but not phmap::flat_hash_set/map. The BTree based
/// version should still have decent performance though.
///
/// @tparam EdgeValue type stored for each edge
template <std::default_initializable EdgeValue>
class SparseGraph {
public:
    using NodeID = size_t;
    //using EdgeID = size_t;

    struct NodeIDPair {
        NodeID from;
        NodeID to;

        std::strong_ordering operator<=>(NodeIDPair other) const noexcept
        {
            return std::tie(from, to) <=> std::tie(other.from, other.to);
        }

        std::weak_ordering operator<=>(NodeID other) const noexcept { return from <=> other; }
    };

    struct edge_t {
        NodeID from;
        NodeID to;
        EdgeValue key;
    };

    struct graph_size_t {
        NodeID nodes;
        NodeID edges;
    };

    struct equal_t {
        using is_transparent = std::true_type;

        bool operator()(edge_t const& lhs, edge_t const& rhs) const noexcept
        {
            return lhs.from == rhs.from && lhs.to == rhs.to;
        }

        bool operator()(NodeID lhs, edge_t const& rhs) const noexcept { return lhs == rhs.from; }
        bool operator()(edge_t const& lhs, NodeID rhs) const noexcept { return lhs.from == rhs; }

        bool operator()(edge_t const& lhs, NodeIDPair const& rhs) const noexcept
        {
            return lhs.from == rhs.from && lhs.to == rhs.to;
        }

        bool operator()(NodeIDPair const& lhs, edge_t const& rhs) const noexcept
        {
            return lhs.from == rhs.from && lhs.to == rhs.to;
        }
    };

    using compare_t = std::less<>;

    struct hash_t {
        using is_transparent = std::true_type;
        std::size_t operator()(NodeID from) const noexcept { return std::hash<NodeID>{}(from); }
        std::size_t operator()(edge_t const& edge) const noexcept { return std::hash<NodeID>{}(edge.from); }

        std::size_t operator()(NodeIDPair const& pair) const noexcept
        {
            return std::hash<NodeID>{}(pair.from);
        }
    };

    /// ID type for edges
    //using EdgeID = size_t;

    /// Special ID value for invalid node
    static constexpr NodeID InvalidNodeID = std::numeric_limits<NodeID>::max();

    /// Special ID value for invalid edge
    static constexpr std::nullopt_t InvalidEdgeID = std::nullopt;


private:
    using Key = NodeIDPair;

    /// The array containing the actual nodes
    //Util::HashSet<edge_t, hash_t, equal_t> m_edges;
    Util::BTreeMap<Key, EdgeValue, compare_t> m_edges;
    size_t no_nodes_created = 0;

    /// Number of edges.

public:
    /// Clear the graph
    void clear()
    {
        m_edges.clear();
        no_nodes_created = 0;
    }

    void reserve(size_t capacity)
    {
        m_edges.reserve(capacity * 2);
    }

    graph_size_t size() const noexcept
    {
        graph_size_t size;
        size.edges = m_edges.size();
        size.nodes = 0;
        for (auto const& edge : m_edges) size.nodes = std::max(size.nodes, edge.from);
        return size;
    }

    /** Return number of nodes. Note that nodes are never removed from edge_map, so the number of nodes is not decremented.
    If you use the derived class AMGraph3D and call cleanup, the number of nodes will  be the actual number. */
    // FIXME: just return the actual number
    [[nodiscard]]
    size_t no_nodes() const { return no_nodes_created; }

    /** Return number of edges created. Note that if nodes were disconnected thereby removing an edge,
     this number is not decremented. If you use the derived class AMGraph3D and call cleanup, the number
     of edges will  be the actual number. */
    [[nodiscard]]
    size_t no_edges() const { return m_edges.size() / 2; }

    /// Return whether the node is valid (i.e. in the graph)
    // FIXME: just return the actual validity
    // bool valid_node_id(NodeID n) const { return n < no_nodes();}

    /// Return whether an edge is valid (i.e. in the graph)
    // // FIXME: just return the actual validity
    // bool valid_edge_id(EdgeID e) const { return e < no_edges_created;}

    /// Returns true if the graph contains no nodes, false otherwise
    [[nodiscard]]
    bool empty() const { return no_nodes() == 0; }

    /// Add a node to the graph
    NodeID add_node()
    {
        auto id = no_nodes_created;
        no_nodes_created += 1;
        return id;
    }

    void connect_nodes(NodeID from, NodeID to, EdgeValue value) noexcept(false)
    {
        if (from == to) return; // Skip self-loop

        m_edges.emplace(Key {from, to}, value);
        m_edges.emplace(Key {to, from}, value);
    }

    std::optional<EdgeValue> find_edge(NodeID from, NodeID to) const noexcept
    {
        if (auto it = m_edges.find({from, to}); it != m_edges.end()) return it->second; //it->key;
        return std::nullopt;
    }

    void remove_edge(NodeID from, NodeID to) noexcept
    {
        m_edges.erase(Key (from, to));
        m_edges.erase(Key (to, from));
    }

    void compact() noexcept(false)
    {
        // Erasing does not trigger a rehash, so we do it manually:
        m_edges.rehash(0);
    }

    std::ranges::input_range auto all_edges_lazy() const noexcept
    {
        return std::ranges::subrange(m_edges.cbegin(), m_edges.cend())
        | std::views::filter([](const auto& edge) -> bool { return edge.first.from < edge.first.to; });
        //| std::views::filter([](const Key& edge) -> bool { return edge.from < edge.to; });
    }

    std::ranges::input_range auto edges_lazy(NodeID from) const noexcept
    {
        auto [begin, end] = m_edges.equal_range(from);
        return std::ranges::subrange(begin, end);
    }

    std::ranges::input_range auto neighbors_lazy(NodeID from) const noexcept
    {
        return edges_lazy(from) | std::views::keys | std::views::transform([](const Key& edge) { return edge.to; });
    }

    [[nodiscard]] std::ranges::input_range
    auto shared_neighbors_lazy(const NodeID n0, const NodeID n1) const
    {
        auto neighbors0 = neighbors_lazy(n0);
        //auto neighbors1 = neighbors_lazy(n1);
        return neighbors0
            | std::views::filter([this, n1](const NodeID& n0_n) { return m_edges.contains(NodeIDPair{n0_n, n1}); });
    }

    template <typename visitor_type_>
    void for_edges(NodeID from, visitor_type_ visitor) const noexcept
    {
        for (auto& edge : edges_lazy(from)) { visitor(edge.from, edge.to, edge.weight); }
    }

    /// The range returned can be used in range based for loops over all node ids
    [[nodiscard]]
    Util::Range node_ids() const { return {0, no_nodes_created};}
    //
    // /// The range returned can be used in range based for loops over all node ids
    // const Util::Range edge_ids() const { return Util::Range(0,no_edges_created);}

    /// Return the number of edges incident on a given node.
    [[nodiscard]]
    size_t valence(NodeID n) const { return std::ranges::distance(edges_lazy(n)); }
};

// FIXME: single generic connected components implementation
template <typename EdgeKey>
std::vector<HashSet<typename SparseGraph<EdgeKey>::NodeID>>
connected_components(const SparseGraph<EdgeKey>& g, const HashSet<typename SparseGraph<EdgeKey>::NodeID>& s) {
    using NodeID = SparseGraph<EdgeKey>::NodeID;
    using NodeSet = HashSet<typename SparseGraph<EdgeKey>::NodeID>;
    using NodeQueue = std::queue<NodeID>;

    NodeSet s_visited;
    std::vector<NodeSet> component_vec;
    for(auto nf0 : s) {
        if(s_visited.count(nf0)==0)
        {
            NodeQueue Q;
            Q.push(nf0);
            s_visited.insert(nf0);
            NodeSet component;
            while(!Q.empty())
            {
                NodeID nf = Q.front();
                Q.pop();
                component.insert(nf);
                for(auto nnf: g.neighbors_lazy(nf))
                    if(s.count(nnf) >0 && s_visited.count(nnf)==0) {
                        Q.push(nnf);
                        s_visited.insert(nnf);
                    }
            }
            component_vec.push_back(component);
        }
    }
    return component_vec;
}
}

#endif //GEL_FLATGRAPH_H
