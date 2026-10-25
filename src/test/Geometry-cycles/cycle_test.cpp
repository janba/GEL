#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest.h>

#include <GEL/Geometry/Graph.h>
#include <GEL/Geometry/JohnsonCycles.h>

#include <algorithm>
#include <random>
#include <set>
#include <vector>

using Geometry::AMGraph;
using Geometry::ElementaryCycle;
using NodeID = AMGraph::NodeID;

namespace {

AMGraph graph_with_edges(size_t node_count, const std::vector<std::pair<NodeID, NodeID>>& edges)
{
    AMGraph graph;
    for (size_t index = 0; index < node_count; ++index)
        graph.add_node();
    for (const auto& edge : edges)
        graph.connect_nodes(edge.first, edge.second);
    return graph;
}

std::vector<NodeID> canonical_cycle(std::vector<NodeID> cycle)
{
    const auto least = std::min_element(cycle.begin(), cycle.end());
    std::rotate(cycle.begin(), least, cycle.end());
    if (cycle.size() > 2 && cycle.back() < cycle[1])
        std::reverse(cycle.begin() + 1, cycle.end());
    return cycle;
}

std::set<std::vector<NodeID>> oracle_cycles(const AMGraph& graph, int max_edges)
{
    std::set<std::vector<NodeID>> found;
    if (max_edges < 3)
        return found;
    const size_t node_count = graph.no_nodes();
    std::vector<char> used(node_count, 0);
    std::vector<NodeID> path;

    const auto walk = [&](auto&& self, NodeID start, NodeID vertex, int edges) -> void {
        for (const NodeID neighbor : graph.neighbors(vertex)) {
            if (neighbor == start) {
                if (edges + 1 >= 3 && edges + 1 <= max_edges)
                    found.insert(canonical_cycle(path));
                continue;
            }
            if (used[neighbor] || edges + 1 >= max_edges)
                continue;
            used[neighbor] = 1;
            path.push_back(neighbor);
            self(self, start, neighbor, edges + 1);
            path.pop_back();
            used[neighbor] = 0;
        }
    };

    for (NodeID start = 0; start < node_count; ++start) {
        std::fill(used.begin(), used.end(), 0);
        used[start] = 1;
        path = {start};
        walk(walk, start, start, 0);
    }
    return found;
}

std::set<std::vector<NodeID>> listed_cycles(const std::vector<ElementaryCycle>& cycles)
{
    std::set<std::vector<NodeID>> found;
    for (const ElementaryCycle& cycle : cycles)
        found.insert(canonical_cycle(cycle.vertices));
    return found;
}

void expect_same_cycles(const AMGraph& graph, int max_edges)
{
    const std::vector<ElementaryCycle> cycles = Geometry::elementary_cycles(graph, max_edges);
    const std::set<std::vector<NodeID>> found = listed_cycles(cycles);
    CHECK(found.size() == cycles.size());
    CHECK(found == oracle_cycles(graph, max_edges));
    for (const ElementaryCycle& cycle : cycles) {
        REQUIRE(cycle.vertices.size() >= 3);
        CHECK(cycle.vertices.size() <= static_cast<size_t>(max_edges));
        CHECK(cycle.start == cycle.vertices.front());
        CHECK(cycle.start == *std::min_element(cycle.vertices.begin(), cycle.vertices.end()));
        CHECK(cycle.vertices[1] < cycle.vertices.back());
        std::set<NodeID> unique(cycle.vertices.begin(), cycle.vertices.end());
        CHECK(unique.size() == cycle.vertices.size());
        for (size_t index = 0; index < cycle.vertices.size(); ++index) {
            const NodeID from = cycle.vertices[index];
            const NodeID to = cycle.vertices[(index + 1) % cycle.vertices.size()];
            CHECK(graph.find_edge(from, to) != AMGraph::InvalidEdgeID);
        }
    }
}

}

TEST_CASE("no cycles below length 3")
{
    const AMGraph empty;
    CHECK(Geometry::elementary_cycles(empty, 10).empty());

    AMGraph loop;
    loop.add_node();
    loop.connect_nodes(0, 0);
    CHECK(Geometry::elementary_cycles(loop, 10).empty());

    const AMGraph triangle = graph_with_edges(3, {{0, 1}, {1, 2}, {2, 0}});
    CHECK(Geometry::elementary_cycles(triangle, 2).empty());
    const auto one = Geometry::elementary_cycles(triangle, 3);
    REQUIRE(one.size() == 1);
    CHECK(one[0].vertices == std::vector<NodeID>{0, 1, 2});
    CHECK(one[0].start == 0);
}

TEST_CASE("square with a diagonal keeps both triangles and the 4-cycle")
{
    const AMGraph square = graph_with_edges(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}});
    expect_same_cycles(square, 3);
    expect_same_cycles(square, 4);
    const auto bounded = listed_cycles(Geometry::elementary_cycles(square, 3));
    const auto full = listed_cycles(Geometry::elementary_cycles(square, 4));
    CHECK(bounded.size() == 2);
    CHECK(full.size() == 3);
}

TEST_CASE("complete graphs match the elementary-cycle counts")
{
    std::vector<std::pair<NodeID, NodeID>> edges;
    for (NodeID from = 0; from < 4; ++from)
        for (NodeID to = from + 1; to < 4; ++to)
            edges.emplace_back(from, to);
    const AMGraph complete = graph_with_edges(4, edges);
    CHECK(Geometry::elementary_cycles(complete, 3).size() == 4);
    CHECK(Geometry::elementary_cycles(complete, 4).size() == 7);
    expect_same_cycles(complete, 4);

    edges.clear();
    for (NodeID from = 0; from < 6; ++from)
        for (NodeID to = from + 1; to < 6; ++to)
            edges.emplace_back(from, to);
    const AMGraph larger = graph_with_edges(6, edges);
    for (int bound = 3; bound <= 6; ++bound)
        expect_same_cycles(larger, bound);
}

TEST_CASE("a failed long path does not hide a shorter cycle")
{
    const AMGraph graph = graph_with_edges(7, {
        {0, 1}, {1, 2}, {2, 6}, {6, 3},
        {0, 4}, {4, 3}, {3, 5}, {5, 0},
    });
    for (int bound = 3; bound <= 7; ++bound)
        expect_same_cycles(graph, bound);
}

TEST_CASE("bridge and grid")
{
    const AMGraph bridge = graph_with_edges(6, {
        {0, 1}, {1, 2}, {2, 0},
        {2, 3},
        {3, 4}, {4, 5}, {5, 3},
    });
    expect_same_cycles(bridge, 6);

    std::vector<std::pair<NodeID, NodeID>> grid;
    constexpr NodeID width = 4;
    for (NodeID y = 0; y < width; ++y) {
        for (NodeID x = 0; x < width; ++x) {
            const NodeID node = y * width + x;
            if (x + 1 < width)
                grid.emplace_back(node, node + 1);
            if (y + 1 < width)
                grid.emplace_back(node, node + width);
        }
    }
    const AMGraph lattice = graph_with_edges(width * width, grid);
    for (int bound = 4; bound <= 8; ++bound)
        expect_same_cycles(lattice, bound);
}

TEST_CASE("random graphs agree with exhaustive search")
{
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> coin(0.0, 1.0);
    for (int trial = 0; trial < 20; ++trial) {
        constexpr size_t node_count = 8;
        std::vector<std::pair<NodeID, NodeID>> edges;
        for (NodeID from = 0; from < node_count; ++from) {
            for (NodeID to = from + 1; to < node_count; ++to) {
                if (coin(rng) < 0.45)
                    edges.emplace_back(from, to);
            }
        }
        const AMGraph graph = graph_with_edges(node_count, edges);
        for (int bound = 3; bound <= 8; ++bound)
            expect_same_cycles(graph, bound);
    }
}
