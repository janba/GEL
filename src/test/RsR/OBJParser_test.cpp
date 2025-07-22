//
// Created by Cem Akarsubasi on 5/21/25.
//

#include <GEL/Util/RawObj.h>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>

using namespace Util;
using namespace Util::Combinators;

static constexpr auto FILE_CAPITAL_A = "../../../../data/PointClouds/Capital_A.obj";
static constexpr auto FILE_BUNNY_SIMPLE = "../../../../data/bunny.obj";
static constexpr auto FILE_BUNNY_COMPLEX = "../../../../data/PointClouds/bun_complete.obj";

TEST_CASE("ignore_spaces")
{
    constexpr auto test = "    a";
    constexpr auto test2 = "     ";

    auto view = std::string_view(test);
    auto view2 = std::string_view(test2);
    SUBCASE("simple") {
        CHECK_EQ(ignore_spaces(view), true);
        CHECK_EQ(view, "a");
    }
    SUBCASE("ignore all") {
        CHECK_EQ(ignore_spaces(view2), true);
        CHECK_EQ(view2, "");
    }
}

TEST_CASE("parse_string")
{
    constexpr auto test = "hello";
    auto view = std::string_view(test);

    SUBCASE("simple") {
        CHECK_EQ(parse_string("hel", view), true);
        CHECK_EQ(view, "lo");
    }
}

TEST_CASE("parse_float")
{
    constexpr auto test = "1.5";
    auto view = std::string_view(test);
    SUBCASE("simple") {
        CHECK_EQ(parse_float(view), 1.5);
        CHECK_EQ(parse_float(view), std::nullopt);
    }

    constexpr auto test2 = " 1.5";
    auto view2 = std::string_view(test2);
    SUBCASE("leading space") {
        CHECK_EQ(parse_float(view2), std::nullopt);
        ignore_spaces(view2);
        CHECK_EQ(parse_float(view2), 1.5);
    }
}

TEST_CASE("parse triplet")
{
    constexpr auto test = "v 1.0 0.0 -1.0";
    constexpr auto test2 = "vn 1.0 0.0 -1.0";
    auto view = std::string_view(test);
    auto view2 = std::string_view(test2);
    SUBCASE("simple")
    {
        CHECK_EQ(parse_prefix_then_float_triplet("v", view), CGLA::Vec3d(1.0, 0.0, -1.0));
        CHECK_EQ(parse_prefix_then_float_triplet("v", view2), std::nullopt);
    }
}

bool float_eq(const double a, const double b, const float eps = 1.0e-6)
{
    return std::abs(a - b) < eps;
}

// TEST_CASE("read raw obj")
// {
//     auto file_path = "../../../../data/bunny.obj";
//     auto file_path_str = std::string(file_path);
//     auto robj = read_raw_obj(file_path);
//     auto pc = read_obj(file_path_str);
//     CHECK_EQ(robj.vertices.size(), pc.vertices.size());
//     CHECK_EQ(robj.normals.size(), pc.normals.size());
//
//     std::cout << pc.vertices.size() << std::endl;
// }