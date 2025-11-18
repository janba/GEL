//
// Created by arkeo on 4/25/25.
//

#include <vector>
#include <functional>
#include <thread>
#include <iostream>
#include <numeric>
#include <algorithm>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include <nanobench.h>

#include <GEL/Util/ParallelAdapters.h>

using Util::detail::ThreadPool;
using namespace Util::detail::Parallel;
using namespace Util::detail;

// bool test_empty()
// {
//     ThreadPool pool(2);
//     const auto empty = map(pool, std::vector<int>{}, [](auto x) { return x; });
//     return empty.empty();
// }

static constexpr auto problem_size = 10000;

inline float non_trivial_function(const int x)
{
    float total = 0.0;
    for (int i = 0; i < 500; ++i) {
        total += 1.0f / static_cast<float>(i+1) * std::cos(static_cast<float>(x)/static_cast<float>(i+1));
    }
    return total;
}

bool test_oracle()
{
    constexpr auto size = 10;
    ThreadPool pool(3);
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::vector<float> vout2(size);
    std::iota(v.begin(), v.end(), 0);
    std::ranges::transform(v.cbegin(), v.cend(), vout.begin(), non_trivial_function);
    const auto right = map(pool, v, std::move(vout2), non_trivial_function);
    CHECK_EQ(vout, right);
    return true;
}

// TEST_CASE("Destructor doesn't block")
// {
//     ThreadPool pool(2);
// }
//
// TEST_CASE("Inert waitAll")
// {
//     ThreadPool pool(2);
//     pool.waitAll();
// }

// TODO: maybe reenable them in the future
// TEST_CASE("deadlock_test1")
// {
//     for (int i = 0; i < 100000; ++i) {
//         CHECK(test_empty());
//         //std::cout << "DT1  " << i << "\n";
//     }
// }
//
// TEST_CASE("deadlock_test2")
// {
//     for (int i = 0; i < 100000; ++i) {
//         test_oracle();
//         //std::cout << "DT2 " << i << "\n";
//     }
// }



TEST_CASE("parallel_map_performance_singlethread")
{
    constexpr auto size = problem_size;
    ThreadPool pool(1);
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("parallel_map_singlethreaded", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            map(pool, v, std::move(vout), non_trivial_function)
        );
    });
}

TEST_CASE("parallel_map_performance_2_threads")
{
    constexpr auto size = problem_size;
    ThreadPool pool(2);
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("parallel_map_2_threads", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            map(pool, v, std::move(vout), non_trivial_function)
        );
    });
}

TEST_CASE("parallel_map_performance_4_threads")
{
    constexpr auto size = problem_size;
    ThreadPool pool(4);
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("parallel_map_4_threads", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            map(pool, v, std::move(vout), non_trivial_function)
        );
    });
}

TEST_CASE("parallel_map_performance_8_threads")
{
    constexpr auto size = problem_size;
    ThreadPool pool(8);
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("parallel_map_8_threads", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            map(pool, v, std::move(vout), non_trivial_function)
        );
    });
}

TEST_CASE("transform_performance")
{
    constexpr auto size = problem_size;
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("std::transform", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            std::transform(v.cbegin(), v.cend(), vout.begin(), non_trivial_function)
        );
    });
}

TEST_CASE("ranges::transform_performance")
{
    constexpr auto size = problem_size;
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    ankerl::nanobench::Bench().run("std::ranges::transform", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            std::ranges::transform(v.cbegin(), v.cend(), vout.begin(), non_trivial_function)
        );
    });
}

TEST_CASE("for_loop_performance")
{
    constexpr auto size = problem_size;
    std::vector<int> v(size);
    std::vector<float> vout(size);
    std::iota(v.begin(), v.end(), 0);
    auto for_loop = [&]()
    {
        for (auto i = 0; i < vout.size(); ++i)
        {
            vout[i] = non_trivial_function(i);
        }
        return 0;
    };

    ankerl::nanobench::Bench().run("for_loop", [&]()
    {
        ankerl::nanobench::doNotOptimizeAway(
            for_loop()
        );
    });
}

auto static_func(int in) -> double { return static_cast<double>(in); };

struct IncompatibleStruct {};

template <typename Ty>
auto identity(Ty in) -> Ty { return in;}

TEST_CASE("parallel_adapter_callability")
{
    ThreadPool pool(1);
    std::vector<int> in(0);
    std::vector<double> out(0);
    const auto& in_ref = in;
    auto& out_ref = out;
    auto in_rref = std::move(in);
    auto out_rref = std::move(out);
    auto lambda = [&](auto x) { return x; };
    // auto r0 = map(pool, in, lambda);
    // auto r1 = map(pool, in_ref, lambda);
    // auto r2 = map(pool, in_rref, lambda);
    // auto r3 = map(pool, std::move(in), lambda); // pointless but allowed
    auto& r4 = map(pool, in, out_ref, lambda);     // return is the same object as out_ref
    auto& r5 = map(pool, in_ref, out_ref, lambda); // return is the same object as out_ref
    auto& r6 = map(pool, in_rref, out_ref, lambda);// return is the same object as out_ref

    auto r7 = map(pool, in, std::move(out), lambda);     // preferred passing. out is moved to return
    auto r8 = map(pool, in_ref, std::move(out), lambda); // moved into the return value
    auto r9 = map(pool, in_rref, std::move(out), lambda);
    auto& r10 = map(pool, in, out_rref, lambda);
    auto& r11 = map(pool, in_ref, out_rref, lambda);
    auto& r12 = map(pool, in_rref, out_rref, lambda);

    auto& lambda_ref = lambda;
    auto&& lambda_rref = lambda;
    auto r13 = map(pool, in, out_rref, lambda_ref);

    std::function<double(int)> lambda_wrap = lambda;
    auto& lambda_wrap_ref = lambda_wrap;
    auto&& lambda_wrap_rref = std::move(lambda);
    auto r14 = map(pool, in, out_rref, lambda_wrap);
    auto r15 = map(pool, in, out_rref, lambda_wrap_ref);
    auto r16 = map(pool, in, out_rref, lambda_wrap_rref);
    auto r17 = map(pool, in, out_rref, static_func);

    // auto r18 = map(pool, in_ref, identity<double>);
}

TEST_CASE("parallel_map_correctness")
{
    ThreadPool pool(3);
    const auto v = [] {
        std::vector<int> v(10);
        std::iota(v.begin(), v.end(), 0);
        return v;
    }();
    const auto times_two = [](auto x){ return x * 2; };
    std::vector<int> out;

    auto& out2 = map(pool, v, out, times_two);
    CHECK(&out2 == &out);
    for (auto id = 0; id < v.size(); ++id) {
        CHECK(out.at(id) == 2 * id);
    }
}

TEST_CASE("parallel_enumerate_map_correctness")
{
    ThreadPool pool(3);
    const auto v = [] {
        std::vector<int> v(10);
        std::iota(v.begin(), v.end(), 0);
        return v;
    }();
    const auto times_two = [](auto idx, auto x){ return idx * x; };
    std::vector<int> out;

    auto& out2 = enumerate_map(pool, v, out, times_two);
    CHECK(&out2 == &out);
    for (auto id = 0; id < v.size(); ++id) {
        CHECK(out.at(id) == id * id);
    }
}

TEST_CASE("parallel_filter_correctness")
{
    ThreadPool pool(3);
    const auto v = [] {
        std::vector<int> v(10);
        std::iota(v.begin(), v.end(), 0);
        return v;
    }();
    const auto is_even = [](auto x){ return x % 2 == 0; };
    std::vector<int> out;

    auto& out2 = filter(pool, v, out, is_even);
    CHECK(&out2 == &out);
    CHECK(out.size() == 5);
    for (auto id = 0; id < out.size(); ++id) {
        CHECK(out.at(id) == id * 2);
    }
}

TEST_CASE("parallel_map_filter_correctness")
{
    // TODO: There might be a race condition that causes a deadlock here
    ThreadPool pool(3);
    const auto v = [] {
        std::vector<int> v(10);
        std::iota(v.begin(), v.end(), 0);
        return v;
    }();
    const auto times_two_is_even = [](int x) -> std::optional<int> {
        if (x % 2 == 0) { return std::make_optional(x * 2); }
        else { return std::nullopt; }
    };
    std::vector<int> out;

    auto& out2 = Parallel::map_filter(pool, v, out, times_two_is_even);
    CHECK(&out2 == &out);
    CHECK(out.size() == 5);
    for (auto id = 0; id < out.size(); ++id) {
        CHECK(out.at(id) == id * 4);
    }
}

TEST_CASE("parallel_enumerate_map_filter_correctness")
{
    ThreadPool pool(3);
    const auto v = [] {
        std::vector<int> v(10);
        std::iota(v.begin(), v.end(), 0);
        return v;
    }();
    const auto times_two_is_even = [](size_t idx, int x) -> std::optional<int> {
        if (x % 2 == 0) { return std::make_optional(x * 2 + idx); }
        else { return std::nullopt; }
    };
    std::vector<int> out;

    auto& out2 = Parallel::enumerate_map_filter(pool, v, out, times_two_is_even);
    CHECK(&out2 == &out);
    CHECK(out.size() == 5);
    for (auto id = 0; id < out.size(); ++id) {
        CHECK(out.at(id) == id * 6);
    }
}

TEST_CASE("parallel_enumerate_map2_filter2_correctness")
{
    ThreadPool pool(3);
    const auto v1 = [] {
        std::vector<int> v(10);
        std::iota(v.begin(), v.end(), 0);
        return v;
    }();
    const auto v2 = [] {
        std::vector<int> v(10);
        std::iota(v.begin(), v.end(), 0);
        return v;
    }();
    const auto times_two_is_even = [](size_t idx, int x, int y) -> std::optional<std::pair<int, int>> {
        if (x % 2 == 0) { return std::make_optional(std::make_pair(x, y)); }
        else { return std::nullopt; }
    };
    std::vector<int> out1;
    std::vector<int> out2;
    auto& out1_ref = out1;
    auto& out2_ref = out2;

    // TODO: map_filter should check if the passed function actually returns an optional
    Parallel::enumerate_map2_filter2(pool, v1, v2, out1, out2, times_two_is_even);
    //CHECK(&out2 == &out);
    CHECK(out1.size() == 5);
    CHECK(out2.size() == 5);

    for (auto id = 0; id < out1.size(); ++id) {
        CHECK(out1.at(id) == id * 2);
    }
    for (auto id = 0; id < out2.size(); ++id) {
        CHECK(out2.at(id) == id * 2);
    }
}