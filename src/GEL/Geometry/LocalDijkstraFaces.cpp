/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/Geometry/LocalDijkstraFaces.h>

#include <algorithm>
#include <cstddef>
#include <queue>
#include <unordered_set>
#include <vector>

namespace Geometry {
namespace {

using NodeID = AMGraph::NodeID;

struct CycleHash {
    size_t operator()(const std::vector<NodeID>& cycle) const noexcept
    {
        size_t h = cycle.size();
        for (const NodeID id : cycle)
            h ^= id + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

void ancestors_of(NodeID node, NodeID root, const std::vector<NodeID>& parent, int max_edges, std::vector<NodeID>& path)
{
    path.clear();
    for (int guard = 0; guard <= max_edges && node != AMGraph::InvalidNodeID; ++guard) {
        path.push_back(node);
        if (node == root)
            break;
        node = parent[node];
    }
}

void cycle_through(const std::vector<NodeID>& up_u, const std::vector<NodeID>& up_v, std::vector<NodeID>& cycle)
{
    cycle.clear();
    size_t at_u = 0;
    size_t at_v = 0;
    bool found = false;
    for (size_t j = 0; j < up_v.size() && !found; ++j) {
        for (size_t i = 0; i < up_u.size(); ++i) {
            if (up_u[i] != up_v[j])
                continue;
            at_u = i;
            at_v = j;
            found = true;
            break;
        }
    }
    if (!found)
        return;

    cycle.reserve(at_u + at_v + 1);
    for (size_t i = 0; i <= at_u; ++i)
        cycle.push_back(up_u[i]);
    size_t j = at_v;
    while (j > 0) {
        --j;
        cycle.push_back(up_v[j]);
    }
}

std::vector<DijkstraFace> faces_from_trees(const AMGraph3D& graph, int max_edges)
{
    const size_t node_count = graph.no_nodes();
    std::unordered_set<std::vector<NodeID>, CycleHash> seen;
    std::vector<DijkstraFace> faces;

    std::vector<int> depth(node_count, -1);
    std::vector<NodeID> parent(node_count, AMGraph::InvalidNodeID);
    std::vector<NodeID> ball;
    std::vector<NodeID> up_node;
    std::vector<NodeID> up_neighbor;
    std::vector<NodeID> cycle;
    std::vector<NodeID> key;
    up_node.reserve(static_cast<size_t>(max_edges) + 1);
    up_neighbor.reserve(static_cast<size_t>(max_edges) + 1);
    key.reserve(static_cast<size_t>(max_edges) * 2 + 1);
    std::queue<NodeID> queue;

    for (const NodeID root : graph.node_ids()) {
        std::fill(depth.begin(), depth.end(), -1);
        ball.clear();
        depth[root] = 0;
        parent[root] = AMGraph::InvalidNodeID;
        queue.push(root);
        ball.push_back(root);

        while (!queue.empty()) {
            const NodeID node = queue.front();
            queue.pop();
            if (depth[node] >= max_edges)
                continue;
            for (const NodeID neighbor : graph.neighbors_lazy(node)) {
                if (depth[neighbor] >= 0)
                    continue;
                depth[neighbor] = depth[node] + 1;
                parent[neighbor] = node;
                queue.push(neighbor);
                ball.push_back(neighbor);
            }
        }

        for (const NodeID node : ball) {
            for (const NodeID neighbor : graph.neighbors_lazy(node)) {
                if (neighbor <= node || depth[neighbor] < 0)
                    continue;
                if (parent[node] == neighbor || parent[neighbor] == node)
                    continue;
                ancestors_of(node, root, parent, max_edges, up_node);
                ancestors_of(neighbor, root, parent, max_edges, up_neighbor);
                cycle_through(up_node, up_neighbor, cycle);
                if (cycle.size() < 3)
                    continue;
                key.assign(cycle.begin(), cycle.end());
                std::sort(key.begin(), key.end());
                if (seen.find(key) != seen.end())
                    continue;
                seen.insert(key);
                faces.push_back(DijkstraFace{std::move(cycle), root});
            }
        }
    }
    return faces;
}

} // namespace

std::vector<DijkstraFace> local_dijkstra_faces(const AMGraph3D& graph, int max_edges)
{
    if (max_edges < 1 || graph.no_nodes() == 0)
        return {};
    return faces_from_trees(graph, max_edges);
}

} // namespace Geometry
