/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <GEL/Geometry/JohnsonCycles.h>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace Geometry {
namespace {

using NodeID = AMGraph::NodeID;

/// Johnson's search. `lock[v]` is the longest path, counted in vertices and
/// including v, that may still enter v. A failed search leaves the lock at the
/// length just tried, so a longer path does not retry v. A cycle raises the
/// lock, and the B-list raises vertices that reached v and failed, one step
/// tighter for each edge of that trail. A walk of two vertices is the
/// undirected edge just used, read backwards, and is not a cycle.
class CycleSearch {
public:
    explicit CycleSearch(const AMGraph& graph)
        : node_count(graph.no_nodes())
    {
        adjacency.resize(node_count);
        for (const NodeID node : graph.node_ids()) {
            const auto neighbors = graph.neighbors_lazy(node);
            adjacency[node].assign(neighbors.begin(), neighbors.end());
        }
    }

    std::vector<ElementaryCycle> run(int max_edges)
    {
        std::vector<ElementaryCycle> cycles;
        if (max_edges < 3 || node_count == 0)
            return cycles;

        bound = max_edges;
        lock.assign(node_count, bound);
        in_path.assign(node_count, 0);
        blist.assign(node_count, {});
        stamp.assign(node_count, 0);
        epoch = 1;
        touched.clear();
        path.clear();
        frames.clear();

        for (NodeID origin_node = 0; origin_node < node_count; ++origin_node) {
            if (neighbors_above(origin_node) < 2)
                continue;
            origin = origin_node;
            circuit_from_origin(cycles);
            reset_touched();
        }
        return cycles;
    }

private:
    struct Frame {
        NodeID vertex = 0;
        size_t next = 0;
        int backward = 0;
    };

    size_t first_kept(NodeID vertex) const
    {
        const auto& row = adjacency[vertex];
        const auto begin = std::lower_bound(row.begin(), row.end(), origin);
        return static_cast<size_t>(begin - row.begin());
    }

    int neighbors_above(NodeID vertex) const
    {
        int count = 0;
        for (const NodeID neighbor : adjacency[vertex]) {
            if (neighbor > vertex && ++count == 2)
                break;
        }
        return count;
    }

    void touch(NodeID vertex)
    {
        if (stamp[vertex] == epoch)
            return;
        stamp[vertex] = epoch;
        touched.push_back(vertex);
    }

    void reset_touched()
    {
        for (const NodeID vertex : touched) {
            lock[vertex] = bound;
            in_path[vertex] = 0;
            blist[vertex].clear();
        }
        touched.clear();
        if (++epoch == 0) {
            std::fill(stamp.begin(), stamp.end(), 0);
            epoch = 1;
        }
    }

    void remember_failure(NodeID vertex)
    {
        const auto& row = adjacency[vertex];
        for (size_t index = first_kept(vertex); index < row.size(); ++index) {
            const NodeID neighbor = row[index];
            touch(neighbor);
            auto& bucket = blist[neighbor];
            if (std::find(bucket.begin(), bucket.end(), vertex) == bucket.end())
                bucket.push_back(vertex);
        }
    }

    void relax(NodeID vertex, int backward)
    {
        std::vector<std::pair<int, NodeID>> pending;
        pending.emplace_back(backward, vertex);
        while (!pending.empty()) {
            const auto item = pending.back();
            pending.pop_back();
            if (item.first > bound)
                continue;
            touch(item.second);
            const int raised = bound - item.first + 1;
            if (!(lock[item.second] < raised))
                continue;
            lock[item.second] = raised;
            for (const NodeID predecessor : blist[item.second]) {
                if (!in_path[predecessor])
                    pending.emplace_back(item.first + 1, predecessor);
            }
        }
    }

    void circuit_from_origin(std::vector<ElementaryCycle>& cycles)
    {
        path.push_back(origin);
        in_path[origin] = 1;
        lock[origin] = 0;
        touch(origin);
        frames.push_back(Frame{origin, first_kept(origin), bound});

        while (!frames.empty()) {
            Frame frame = frames.back();
            bool advanced = false;
            const auto& row = adjacency[frame.vertex];
            while (frame.next < row.size()) {
                const NodeID neighbor = row[frame.next++];
                if (neighbor == origin) {
                    if (path.size() >= 3 && path[1] < path.back())
                        cycles.push_back(ElementaryCycle{path, origin});
                    if (path.size() >= 3)
                        frame.backward = 1;
                    continue;
                }
                if (path.size() < static_cast<size_t>(lock[neighbor]) && !in_path[neighbor]) {
                    frames.back() = frame;
                    path.push_back(neighbor);
                    in_path[neighbor] = 1;
                    lock[neighbor] = static_cast<int>(path.size());
                    touch(neighbor);
                    frames.push_back(Frame{neighbor, first_kept(neighbor), bound});
                    advanced = true;
                    break;
                }
            }
            if (advanced)
                continue;

            frames.pop_back();
            path.pop_back();
            in_path[frame.vertex] = 0;
            if (!frames.empty())
                frames.back().backward = std::min(frames.back().backward, frame.backward);
            if (frame.backward < bound)
                relax(frame.vertex, frame.backward);
            else
                remember_failure(frame.vertex);
        }
    }

    size_t node_count = 0;
    int bound = 0;
    NodeID origin = 0;
    unsigned epoch = 1;
    std::vector<std::vector<NodeID>> adjacency;
    std::vector<int> lock;
    std::vector<char> in_path;
    std::vector<std::vector<NodeID>> blist;
    std::vector<unsigned> stamp;
    std::vector<NodeID> touched;
    std::vector<NodeID> path;
    std::vector<Frame> frames;
};

} // namespace

std::vector<ElementaryCycle> elementary_cycles(const AMGraph& graph, int max_edges)
{
    return CycleSearch(graph).run(max_edges);
}

} // namespace Geometry
