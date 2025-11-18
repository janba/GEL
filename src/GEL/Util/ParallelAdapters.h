//
// Created by Cem Akarsubasi on 5/5/25.
//
// TODO: examples

#ifndef GEL_UTIL_PARALLELADAPTERS_H
#define GEL_UTIL_PARALLELADAPTERS_H
#pragma once

#include <GEL/Util/ThreadPool.h>
#include <optional>
#include <barrier>
#ifdef _MSC_VER
#include <array> // std::size
#include <limits> // std::numeric_limits
#endif

namespace Util::detail {
    /// Concepts for constraining the inputs of parallel adapters
    namespace Concepts {

        /// A range with random access that knows its size in constant time
        template<typename C> concept SizedRange = std::ranges::sized_range<C> && std::ranges::random_access_range<C>;

        /// A range with random access that knows its size in constant time that can be resized
        template<typename R> concept ResizableRange = SizedRange<R> && requires(R& r, size_t size) {
            { r.reserve(size) };
            { r.resize(size) };
        };

        /// Returns the value inside a range
        template<SizedRange Range> using Item = std::ranges::range_value_t<Range>;

        /// Maintains the value classification of T
        template<typename T> using Forward = decltype(std::forward<T>(std::declval<T>()));

        template<typename F, typename A, typename R> concept UnaryFunction = requires(A a, F f) {
            { std::is_invocable<F, decltype(a)>() };
            { f(a) } -> std::convertible_to<R>;
        };

        template<typename F, typename A, typename B, typename R> concept BinaryFunction = requires(A a, B b, F f) {
            { std::is_invocable<F, decltype(a), decltype(b)>() };
            { f(a, b) } -> std::convertible_to<R>;
        };

        template<typename F, typename A, typename B, typename C, typename R> concept TernaryFunction
                = requires(A a, B b, C c, F f) {
                      { std::is_invocable<F, decltype(a), decltype(b), decltype(c)>() };
                      { f(a, b, c) } -> std::convertible_to<R>;
                  };

        /// @code
        /// Functor: 'a -> 'b where I: SizedRange<Item = 'a>
        /// @endcode
        template<typename F, typename I> concept ForeachFunctor = SizedRange<I> && UnaryFunction<F, Item<I>, void>;

        /// @code
        /// Functor: 'a -> 'b where I: SizedRange<Item = 'a>, O: SizedRange<Item = 'b>
        /// @endcode
        template<typename F, typename I, typename O> concept MapFunctor
                = SizedRange<I> && SizedRange<O> && UnaryFunction<F, Item<I>, Item<O>>;

        /// @code
        /// Functor: size_t -> 'a -> 'b where I: SizedRange<Item = 'a>, O: SizedRange<Item = 'b>
        /// @endcode
        template<typename F, typename I, typename O> concept EnumMapper
                = SizedRange<I> && SizedRange<O> && BinaryFunction<F, size_t, Item<I>, Item<O>>;

        /// @code
        /// Functor: 'a -> 'b -> 'c where I1: SizedRange<Item = 'a>, I2: SizedRange<Item = 'b>
        /// @endcode
        template<typename F, typename I1, typename I2> concept Foreach2Functor
                = SizedRange<I1> && SizedRange<I2> && BinaryFunction<F, Item<I1>, Item<I2>, void>;

        /// @code
        /// Functor: size_t -> 'a -> 'b -> 'c where I1: SizedRange<Item = 'a>, I2: SizedRange<Item = 'b>
        /// @endcode
        template<typename F, typename I1, typename I2> concept EnumForeach2Functor
                = SizedRange<I1> && SizedRange<I2> && TernaryFunction<F, size_t, Item<I1>, Item<I2>, void>;

        /// @code
        /// Functor: 'a -> 'b -> 'c where I1: SizedRange<Item = 'a>, I2: SizedRange<Item = 'b>,  O: SizedRange<Item =
        /// 'b>
        /// @endcode
        template<typename F, typename I1, typename I2, typename O> concept Map2Functor
                = SizedRange<I1> && SizedRange<I2> && SizedRange<O> && BinaryFunction<F, Item<I1>, Item<I2>, Item<O>>;

        /// @code
        /// Functor: 'a -> 'b -> 'c where I1: SizedRange<Item = 'a>, I2: SizedRange<Item = 'b>,  O: SizedRange<Item =
        /// 'b>
        /// @endcode
        template<typename F, typename I1, typename I2, typename O> concept EnumMap2Functor
                = SizedRange<I1> && SizedRange<I2> && SizedRange<O>
                  && TernaryFunction<F, size_t, Item<I1>, Item<I2>, Item<O>>;

        /// @code
        /// Functor: 'a -> 'b where I: SizedRange<Item = 'a>, O: SizedRange<Item = 'b>
        /// @endcode
        template<typename F, typename I, typename O> concept FilterFunctor
                = SizedRange<I> && ResizableRange<O> && UnaryFunction<F, Item<I>, bool>
                  && std::same_as<Item<I>, Item<O>>;

        /// @code
        /// Functor: 'a -> 'b where I: SizedRange<Item = 'a>, O: SizedRange<Item = 'b>
        /// @endcode
        template<typename F, typename I, typename O> concept MapFilterFunctor
                = SizedRange<I> && ResizableRange<O> && UnaryFunction<F, Item<I>, std::optional<Item<O>>>;

        /// @code
        /// Functor: 'a -> 'b where I: SizedRange<Item = 'a>, O: SizedRange<Item = 'b>
        /// @endcode
        template<typename F, typename I, typename O> concept EnumMapFilterFunctor
                = SizedRange<I> && ResizableRange<O> && BinaryFunction<F, size_t, Item<I>, std::optional<Item<O>>>;

        /// @code
        /// Functor: 'a -> 'b where I: SizedRange<Item = 'a>, O: SizedRange<Item = 'b>
        /// @endcode
        template<typename F, typename I1, typename I2, typename O> concept EnumMap2FilterFunctor
                = SizedRange<I1> && SizedRange<I2> && ResizableRange<O>
                  && TernaryFunction<F, size_t, Item<I1>, Item<I2>, std::optional<Item<O>>>;

        /// @code
        /// Functor: 'a -> 'b -> Option Pair 'c * 'd
        ///     where I1: SizedRange<Item = 'a>,
        ///           I2: SizedRange<Item = 'b>,
        ///           O1: ResizableRange<Item = 'c>
        ///           O2: ResizableRange<Item = 'd>
        /// @endcode
        template<typename F, typename I1, typename I2, typename O1, typename O2> concept Map2Filter2Functor
                = SizedRange<I1> && SizedRange<I2> && SizedRange<O1> && SizedRange<O2>
                  && BinaryFunction<F, Item<I1>, Item<I2>, std::optional<std::pair<Item<O1>, Item<O2>>>>;

        /// @code
        /// Functor: size_t -> 'a -> 'b -> Option Pair 'c * 'd
        ///     where I1: SizedRange<Item = 'a>,
        ///           I2: SizedRange<Item = 'b>,
        ///           O1: ResizableRange<Item = 'c>
        ///           O2: ResizableRange<Item = 'd>
        /// @endcode
        template<typename F, typename I1, typename I2, typename O1, typename O2> concept EnumMap2Filter2Functor
                = SizedRange<I1> && SizedRange<I2> && SizedRange<O1> && SizedRange<O2>
                  && TernaryFunction<F, size_t, Item<I1>, Item<I2>, std::optional<std::pair<Item<O1>, Item<O2>>>>;

        template<typename F, typename State, typename I> concept FoldFunctor
                = SizedRange<I> && BinaryFunction<F, State, Item<I>, State>;
    } // namespace Concepts

    /// Helper functions for parallel algorithms
    /// @internal
    namespace ParallelUtil {
        /// @brief variadic function to find the shortest iterator
        /// @tparam Iterator The first iterator
        /// @tparam Iterators Other iterators
        /// @param list1 The first iterator
        /// @param lists Further iterators
        /// @return The length of the shortest iterator
        template<typename Iterator, typename... Iterators>
        constexpr size_t smallest_size(const Iterator& list1, const Iterators&... lists)
        {
            size_t min_size = std::numeric_limits<size_t>::max();
            for (std::array<size_t, sizeof...(lists) + 1> sizes = {list1.size(), std::size(lists...)};
                 auto size: sizes) {
                if (size < min_size) { min_size = size; }
            }
            return min_size;
        }

        /// @brief Ceiling division
        /// @tparam IntegerType1 Integral type
        /// @tparam IntegerType2 Integral type
        /// @param lhs lhs
        /// @param rhs rhs
        /// @return ceil(lhs / div)
        template<typename IntegerType1, typename IntegerType2>
        requires std::is_integral_v<IntegerType1> && std::is_integral_v<IntegerType2>
        constexpr IntegerType1 div_ceil(IntegerType1 lhs, IntegerType2 rhs)
        {
            return lhs / rhs + (lhs % rhs != 0);
        }

        /// @brief Wraps up outputs to a value or a reference pair depending on the value classification
        /// @return either @code std::pair(T1&, T2&) or std::pair(T1, T2) @endcode
        template<typename T1, typename T2> auto make_pair_wrapper(T1&& t1, T2&& t2)
        {
            if constexpr (std::is_lvalue_reference_v<T1> && std::is_lvalue_reference_v<T2>) {
                return std::make_pair(std::ref(t1), std::ref(t2));
            }
            else {
                return std::make_pair(std::forward<T1>(t1), std::forward<T2>(t2));
            }
        }

        // 64 bytes on x86-64 │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ __cacheline_aligned │ ...
        static constexpr std::size_t CACHE_LINE_SIZE = 64; // 64 bytes is as good a guess as any

        struct alignas(CACHE_LINE_SIZE) CacheLineCounter {
            size_t inner;
        };

        /// @brief Moves chunks and resizes the output array for
        /// @param target output array
        /// @param counters the counter vector containing chunk sizes
        /// @param max_chunk_size gap size between the chunks
        template <typename Target>
        void fix_chunks(Target& target, const std::vector<CacheLineCounter>& counters, size_t max_chunk_size)
        {
            size_t total_size = 0;
            for (const auto& counter: counters) { total_size += counter.inner; }

            auto end_of_chunk1 = target.begin() + counters.at(0).inner;
            for (size_t chunk = 1; chunk < counters.size(); ++chunk) {
                auto this_chunk_size = counters[chunk].inner;
                auto chunk_begin = target.begin() + chunk * max_chunk_size;
                auto chunk_end = chunk_begin + this_chunk_size;
                std::copy(chunk_begin, chunk_end, end_of_chunk1);
                end_of_chunk1 += this_chunk_size;
            }
            target.resize(total_size);
        }
    } // namespace ParallelUtil

    /// Parallel algorithms
    namespace Parallel {
        using namespace Concepts;

        /// @brief Runs the given closure for each element of the input iterator
        /// @tparam F a unary function type of the form T -> void
        /// @tparam InputIt An input iterator type that yields elements of T
        /// @param pool The threadpool to use
        /// @param it the input iterator
        /// @param f the function to map over
        ///
        template<SizedRange InputIt, ForeachFunctor<InputIt> F>
        auto foreach (IExecutor& pool, const InputIt& it, F && f) -> void
        {
            const auto pool_size = pool.size();
            const auto work_size = it.size();
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);
            if (work_size == 0) { return; }
            for (size_t i = 0; i < pool_size; ++i) {
                pool.addTask([i, reduced_size, work_size, &it, &f] {
                    const auto max_size = std::min((i + 1) * reduced_size, work_size);
                    for (auto j = i * reduced_size; j < max_size; ++j) { f(it[j]); }
                });
            }
            pool.waitAll();
        }

        /// @brief Runs the given closure for each element of two input iterators
        /// @tparam F a binary function type of the form (T, U) -> void
        /// @tparam InputIt1 An input iterator type that yields elements of T
        /// @tparam InputIt2 An input iterator type that yields elements of U
        /// @param pool The threadpool to use
        /// @param it1 the first input iterator
        /// @param it2 the second input iterator
        /// @param f the function to map over
        template<SizedRange InputIt1, SizedRange InputIt2, Foreach2Functor<InputIt1, InputIt2> F>
        auto foreach2(IExecutor& pool, InputIt1& it1, InputIt2& it2, F&& f) -> void
        {
            const auto pool_size = pool.size();
            const auto work_size = ParallelUtil::smallest_size(it1, it2);
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);
            if (work_size == 0) { return; }
            for (size_t i = 0; i < pool_size; ++i) {
                pool.addTask([i, reduced_size, work_size, &it1, &it2, &f] {
                    const auto max_size = std::min((i + 1) * reduced_size, work_size);
                    for (auto j = i * reduced_size; j < max_size; ++j) { f(it1[j], it2[j]); }
                });
            }
            pool.waitAll();
        }

        /// @brief Runs the given closure for each element of two input iterators
        /// @tparam F a binary function type of the form (T, U) -> void
        /// @tparam InputIt1 An input iterator type that yields elements of T
        /// @tparam InputIt2 An input iterator type that yields elements of U
        /// @param pool The threadpool to use
        /// @param it1 the first input iterator
        /// @param it2 the second input iterator
        /// @param f the function to map over
        template<SizedRange InputIt1, SizedRange InputIt2, EnumForeach2Functor<InputIt1, InputIt1> F>
        auto enumerate_foreach2(IExecutor& pool, InputIt1& it1, InputIt2& it2, F&& f) -> void
        {
            const auto pool_size = pool.size();
            const auto work_size = ParallelUtil::smallest_size(it1, it2);
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);
            if (work_size == 0) { return; }
            for (size_t i = 0; i < pool_size; ++i) {
                pool.addTask([i, reduced_size, work_size, &it1, &it2, &f] {
                    const auto max_size = std::min((i + 1) * reduced_size, work_size);
                    for (auto j = i * reduced_size; j < max_size; ++j) { f(j, it1[j], it2[j]); }
                });
            }
            pool.waitAll();
        }

        /// @brief Map over an iterator.
        /// @tparam F a unary function type of the form T -> U
        /// @tparam InputIt An input iterator type that yields elements of T
        /// @tparam OutputIt An output iterator type that takes elements of type U
        /// @param pool The threadpool to use
        /// @param it the input iterator
        /// @param out the output iterator
        /// @param f the function to map over
        /// @return the forwarded output iterator
        ///
        template<SizedRange InputIt, SizedRange OutputIt, MapFunctor<InputIt, OutputIt> F>
        auto map(IExecutor& pool, const InputIt& it, OutputIt&& out, F&& f) -> Forward<OutputIt>
        {
            const auto pool_size = pool.size();
            const auto work_size = it.size();
            const auto reduced_size = work_size / pool_size + (work_size % pool_size != 0);
            if (work_size == 0) { return std::forward<OutputIt>(out); }
            if (out.size() != work_size) {
                out.reserve(work_size);
                out.resize(work_size);
            }

            for (size_t i = 0; i < pool_size; ++i) {
                pool.addTask([&out, i, reduced_size, work_size, &it, &f] {
                    auto max_size = std::min((i + 1) * reduced_size, work_size);
                    for (auto j = i * reduced_size; j < max_size; ++j) { out[j] = f(it[j]); }
                });
            }
            pool.waitAll();
            return std::forward<OutputIt>(out);
        }

        /// @brief Map over an iterator with the indexes
        /// @tparam F a binary function type of the form (size_t, T) -> U
        /// @tparam InputIt An input iterator type that yields elements of T
        /// @tparam OutputIt An output iterator type that takes elements of type U
        /// @param pool The threadpool to use
        /// @param it the input iterator
        /// @param out the output iterator
        /// @param f the function to map over
        /// @return the forwarded output iterator.
        template<SizedRange InputIt, SizedRange OutputIt, EnumMapper<InputIt, OutputIt> F>
        auto enumerate_map(IExecutor& pool, const InputIt& it, OutputIt&& out, F&& f) -> decltype(out)
        {
            const auto pool_size = pool.size();
            const auto work_size = it.size();
            const auto reduced_size = work_size / pool_size + (work_size % pool_size != 0);
            if (work_size == 0) { return std::forward<OutputIt>(out); }
            if (out.size() != work_size) {
                out.reserve(work_size);
                out.resize(work_size);
            }

            for (size_t i = 0; i < pool_size; ++i) {
                pool.addTask([&out, i, reduced_size, work_size, &it, &f] {
                    auto max_size = std::min((i + 1) * reduced_size, work_size);
                    for (auto j = i * reduced_size; j < max_size; ++j) { out[j] = f(j, it[j]); }
                });
            }
            pool.waitAll();
            return std::forward<OutputIt>(out);
        }

        /// @brief Maps over two input iterators at the same time. If the iterators are of different length,
        /// iterates up to the end of the shorter iterator.
        /// @tparam F a unary function type of the form (T, U) -> V
        /// @tparam InputIt1 An input iterator type that yields elements of T
        /// @tparam InputIt2 An input iterator type that yields elements of U
        /// @tparam OutputIt An output iterator type that takes elements of type V
        /// @param pool The threadpool to use
        /// @param it1 the first input iterator
        /// @param it2 the second input iterator
        /// @param out the output iterator
        /// @param f the function to map over
        /// @return the output iterator by value
        ///
        template<SizedRange InputIt1, SizedRange InputIt2, SizedRange OutputIt,
                 Map2Functor<InputIt1, InputIt2, OutputIt> F>
        auto map2(IExecutor& pool, const InputIt1& it1, const InputIt2& it2, OutputIt&& out, F&& f) -> decltype(out)
        {
            const auto pool_size = pool.size();
            const auto work_size = ParallelUtil::smallest_size(it1, it2);
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);
            if (work_size == 0) { return std::forward<OutputIt>(out); }
            if (out.size() != work_size) {
                out.reserve(work_size);
                out.resize(work_size);
            }

            for (size_t i = 0; i < pool_size; ++i) {
                pool.addTask([&out, i, reduced_size, work_size, &it1, &it2, &f] {
                    auto max_size = std::min((i + 1) * reduced_size, work_size);
                    for (auto j = i * reduced_size; j < max_size; ++j) { out[j] = f(it1[j], it2[j]); }
                });
            }
            pool.waitAll();
            return std::forward<OutputIt>(out);
        }

        /// @brief Map over an iterator with the indexes
        /// @tparam F a ternary function type of the form (size_t, T, U) -> V
        /// @tparam InputIt1 An input iterator type that yields elements of T
        /// @tparam InputIt2 An input iterator type that yields elements of U
        /// @tparam OutputIt An output iterator type that takes elements of type V
        /// @param pool The threadpool to use
        /// @param it1 the first input iterator
        /// @param it2 the second input iterator
        /// @param out the output iterator
        /// @param f the function to map over
        /// @return the forwarded output iterator.
        ///
        template<SizedRange InputIt1, typename InputIt2, SizedRange OutputIt,
                 EnumMap2Functor<InputIt1, InputIt2, OutputIt> F>
        auto enumerate_map2(IExecutor& pool, const InputIt1& it1, const InputIt2& it2, OutputIt&& out, F&& f)
                -> decltype(out)
        {
            const auto pool_size = pool.size();
            const auto work_size = ParallelUtil::smallest_size(it1, it2);
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);
            if (work_size == 0) { return std::forward<OutputIt>(out); }
            if (out.size() != work_size) {
                out.reserve(work_size);
                out.resize(work_size);
            }

            for (size_t i = 0; i < pool_size; ++i) {
                pool.addTask([&out, i, reduced_size, work_size, &it1, &it2, &f] {
                    auto max_size = std::min((i + 1) * reduced_size, work_size);
                    for (auto j = i * reduced_size; j < max_size; ++j) { out[j] = f(j, it1[j], it2[j]); }
                });
            }
            pool.waitAll();
            return std::forward<OutputIt>(out);
        }

        /// @brief perform a filter
        /// @tparam F a unary function predicate of the form T -> bool
        /// @tparam InputIt An input iterator type that yields elements of T
        /// @tparam OutputIt An output iterator type that takes elements of T
        /// @param pool The threadpool to use
        /// @param it the input iterator
        /// @param out the output iterator
        /// @param p the predicate
        /// @return the output iterator by value
        ///
        template<SizedRange InputIt, SizedRange OutputIt, FilterFunctor<InputIt, OutputIt> F>
        auto filter(IExecutor& pool, const InputIt& it, OutputIt&& out, F&& p) -> decltype(out)
        {
            // strategy: we allocate out to the same size as it, every thread has its own counter, and when they are all
            // done, we perform a serial memcpy. memcpy is memory-bound, meaning there is not a lot to gain from
            // parallelization.
            const auto pool_size = pool.size();
            const auto work_size = it.size();
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);
            if (work_size == 0) { return std::forward<OutputIt>(out); }
            if (out.size() != work_size) {
                out.reserve(work_size);
                // safety post-condition: we need to manually resize out to the correct size at the end
                out.resize(work_size);
            }
            using ParallelUtil::CacheLineCounter;
            std::vector<CacheLineCounter> counters(pool_size, CacheLineCounter(0));
            for (size_t thread_id = 0; thread_id < pool_size; ++thread_id) {
                auto& thread_counter = counters[thread_id];
                thread_counter.inner = 0;
                pool.addTask([&out, thread_id, reduced_size, work_size, &it, &p, &thread_counter] {
                    auto max_size = std::min((thread_id + 1) * reduced_size, work_size);
                    auto begin = thread_id * reduced_size;
                    for (auto j = begin; j < max_size; ++j) {
                        const auto& value = it[j];
                        if (p(value)) {
                            out[begin + thread_counter.inner] = value;
                            thread_counter.inner++;
                        }
                    }
                });
            }
            pool.waitAll();

            ParallelUtil::fix_chunks(out, counters, reduced_size);

            return std::forward<OutputIt>(out);
        }

        /// @brief perform a map and a filter operation simultaneously
        /// @tparam F a unary function type of the form T -> std::optional U
        /// @tparam InputIt An input iterator type that yields elements of T
        /// @tparam OutputIt An output iterator type that takes elements of type U
        /// @param pool The threadpool to use
        /// @param it the input iterator
        /// @param out the output iterator
        /// @param f the function to map-filter over
        /// @return the output iterator by value
        ///
        template<SizedRange InputIt, SizedRange OutputIt, MapFilterFunctor<InputIt, OutputIt> F>
        auto map_filter(IExecutor& pool, const InputIt& it, OutputIt&& out, F&& f) -> decltype(out)
        {
            // strategy: we allocate out to the same size as it, every thread has its own counter, and when they are all
            // done, we perform a serial memcpy. memcpy is memory-bound, meaning there is not a lot to gain from
            // parallelization.
            const auto pool_size = pool.size();
            const auto work_size = it.size();
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);
            if (work_size == 0) { return std::forward<OutputIt>(out); }
            if (out.size() != work_size) {
                out.reserve(work_size);
                // safety post-condition: we need to manually resize out to the correct size at the end
                out.resize(work_size);
            }
            using ParallelUtil::CacheLineCounter;
            std::vector<CacheLineCounter> counters(pool_size, CacheLineCounter(0));
            for (auto thread_id = 0; thread_id < pool_size; ++thread_id) {
                auto& thread_counter = counters[thread_id].inner;
                thread_counter = 0;
                pool.addTask([&, thread_id, reduced_size, work_size] {
                    const auto max_size = std::min((thread_id + 1) * reduced_size, work_size);
                    const auto begin = thread_id * reduced_size;
                    for (auto j = begin; j < max_size; ++j) {
                        const auto& input_value = it[j];
                        if (const auto&& result = f(input_value); result.has_value()) {
                            out[begin + thread_counter] = *result;
                            thread_counter++;
                        }
                    }
                });
            }
            pool.waitAll();

            ParallelUtil::fix_chunks(out, counters, reduced_size);

            return std::forward<OutputIt>(out);
        }

        /// @brief perform a map and a filter operation simultaneously over an enumerated iterator
        /// @tparam F a binary function type of the form (size_t, T) -> std::optional U
        /// @tparam InputIt An input iterator type that yields elements of T
        /// @tparam OutputIt An output iterator type that takes elements of type U
        /// @param pool The threadpool to use
        /// @param it the input iterator
        /// @param out the output iterator
        /// @param f the function to map-filter over
        /// @return the output iterator by value
        ///
        template<SizedRange InputIt, SizedRange OutputIt, EnumMapFilterFunctor<InputIt, OutputIt> F>
        auto enumerate_map_filter(IExecutor& pool, const InputIt& it, OutputIt&& out, F&& f) -> decltype(out)
        {
            // strategy: we allocate out to the same size as it, every thread has its own counter, and when they are all
            // done, we perform a serial memcpy. memcpy is memory-bound, meaning there is not a lot to gain from
            // parallelization.
            const auto pool_size = pool.size();
            const auto work_size = it.size();
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);
            if (work_size == 0) { return std::forward<OutputIt>(out); }
            if (out.size() != work_size) {
                out.reserve(work_size);
                // safety post-condition: we need to manually resize out to the correct size at the end
                out.resize(work_size);
            }
            using ParallelUtil::CacheLineCounter;
            std::vector<CacheLineCounter> counters(pool_size, CacheLineCounter(0));
            for (size_t thread_id = 0; thread_id < pool_size; ++thread_id) {
                auto& thread_counter = counters[thread_id];
                thread_counter.inner = 0;
                pool.addTask([&out, thread_id, reduced_size, work_size, &it, &f, &thread_counter] {
                    auto max_size = std::min((thread_id + 1) * reduced_size, work_size);
                    auto begin = thread_id * reduced_size;
                    for (auto j = begin; j < max_size; ++j) {
                        const auto& input_value = it[j];
                        if (const auto&& result = f(j, input_value); result.has_value()) {
                            auto inner = std::move(result.value());
                            out[begin + thread_counter.inner] = inner;
                            thread_counter.inner++;
                        }
                    }
                });
            }
            pool.waitAll();

            ParallelUtil::fix_chunks(out, counters, reduced_size);

            return std::forward<OutputIt>(out);
        }

        /// @brief perform a map and a filter operation simultaneously over two input iterators
        /// @tparam F a binary function type of the form (T, U) -> std::optional V
        /// @tparam InputIt1 An input iterator type that yields elements of T
        /// @tparam InputIt2 An input iterator type that yields elements of U
        /// @tparam OutputIt An output iterator type that takes elements of type V
        /// @param pool The threadpool to use
        /// @param it1 the first input iterator
        /// @param it2 the second input iterator
        /// @param out the output iterator
        /// @param f the function to map-filter over
        /// @return the output iterator by value
        ///
        template<SizedRange InputIt1, SizedRange InputIt2, ResizableRange OutputIt,
                 EnumMap2FilterFunctor<InputIt1, InputIt2, OutputIt> F>
        auto map2_filter(IExecutor& pool, const InputIt1& it1, const InputIt2& it2, OutputIt&& out, F&& f)
                -> decltype(out)
        {
            // strategy: we allocate out to the same size as it, every thread has its own counter, and when they are all
            // done, we perform a serial memcpy. memcpy is memory-bound, meaning there is not a lot to gain from
            // parallelization.
            const auto pool_size = pool.size();
            const auto work_size = ParallelUtil::smallest_size(it1, it2);
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);
            if (work_size == 0) { return std::forward<OutputIt>(out); }
            if (out.size() != work_size) {
                out.reserve(work_size);
                // safety post-condition: we need to manually resize out to the correct size at the end
                out.resize(work_size);
            }
            using ParallelUtil::CacheLineCounter;
            std::vector<CacheLineCounter> counters(pool_size, CacheLineCounter(0));
            for (size_t thread_id = 0; thread_id < pool_size; ++thread_id) {
                auto& thread_counter = counters[thread_id];
                thread_counter.inner = 0;
                pool.addTask([&out, thread_id, reduced_size, work_size, &it1, &f, &thread_counter, &it2] {
                    auto max_size = std::min((thread_id + 1) * reduced_size, work_size);
                    auto begin = thread_id * reduced_size;
                    for (auto j = begin; j < max_size; ++j) {
                        const auto& input_value1 = it1[j];
                        const auto& input_value2 = it2[j];
                        const auto&& result = f(input_value1, input_value2);
                        if (result.has_value()) {
                            auto&& inner = result.value();
                            out[begin + thread_counter.inner] = inner;
                            thread_counter.inner++;
                        }
                    }
                });
            }
            pool.waitAll();

            ParallelUtil::fix_chunks(out, counters, reduced_size);

            return std::forward<OutputIt>(out);
        }

        /// @brief perform a map operation over two input iterators and a filter operation through two output iterators
        /// @tparam F a binary function type of the form (T, U) -> std::optional std::pair V, W
        /// @tparam InputIt1 An input iterator type that yields elements of T
        /// @tparam OutputIt1 An output iterator type that takes elements of type U
        /// @param pool The threadpool to use
        /// @param it1 the first input iterator
        /// @param it2 the second input iterator
        /// @param out1 the first output iterator
        /// @param out2 the second output iterator
        /// @param f the function to map-filter over
        /// @return the pair of output iterators by value
        ///
        template<SizedRange InputIt1, SizedRange InputIt2, ResizableRange OutputIt1, ResizableRange OutputIt2,
                 Map2Filter2Functor<InputIt1, InputIt2, OutputIt1, OutputIt2> F>
        auto map2_filter2(IExecutor& pool, const InputIt1& it1, const InputIt2& it2, OutputIt1&& out1, OutputIt2&& out2,
                          F&& f) -> std::pair<decltype(out1), decltype(out2)>
        {
            // strategy: we allocate out to the same size as it, every thread has its own counter, and when they are all
            // done, we perform a serial memcpy. memcpy is memory-bound, meaning there is not a lot to gain from
            // parallelization.
            const auto pool_size = pool.size();
            const auto work_size = ParallelUtil::smallest_size(it1, it2);
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);
            if (work_size == 0) { return std::make_pair(std::forward<OutputIt1>(out1), std::forward<OutputIt2>(out2)); }
            if (out1.size() != work_size) {
                out1.reserve(work_size);
                // safety post-condition: we need to manually resize out to the correct size at the end
                out1.resize(work_size);
            }
            if (out2.size() != work_size) {
                out2.reserve(work_size);
                // safety post-condition: we need to manually resize out to the correct size at the end
                out2.resize(work_size);
            }
            using ParallelUtil::CacheLineCounter;
            std::vector<CacheLineCounter> counters(pool_size, CacheLineCounter(0));
            for (size_t thread_id = 0; thread_id < pool_size; ++thread_id) {
                auto& thread_counter = counters[thread_id];
                thread_counter.inner = 0;
                pool.addTask([&out1, thread_id, reduced_size, work_size, &it1, &f, &thread_counter, &it2, out2] {
                    auto max_size = std::min((thread_id + 1) * reduced_size, work_size);
                    auto begin = thread_id * reduced_size;
                    for (auto j = begin; j < max_size; ++j) {
                        const auto& input_value1 = it1[j];
                        const auto& input_value2 = it2[j];
                        const auto result = f(input_value1, input_value2);
                        if (result.has_value()) {
                            auto [inner1, inner2] = result.value();
                            out1[begin + thread_counter.inner] = std::move(inner1);
                            out2[begin + thread_counter.inner] = std::move(inner2);
                            thread_counter.inner++;
                        }
                    }
                });
            }
            pool.waitAll();


            ParallelUtil::fix_chunks(out1, counters, reduced_size);
            ParallelUtil::fix_chunks(out2, counters, reduced_size);

            return std::make_pair(std::forward<OutputIt1>(out1), std::forward<OutputIt2>(out2));
        }

        /// @brief perform a map operation over two input iterators and a filter operation through two output iterators
        /// @tparam F a ternary function type of the form (size_t, T, U) -> std::optional std::pair V, W
        /// @tparam I1 An input iterator type that yields elements of T
        /// @tparam O1 An output iterator type that takes elements of type U
        /// @param pool The threadpool to use
        /// @param it1 the first input iterator
        /// @param it2 the second input iterator
        /// @param out1 the first output iterator
        /// @param out2 the second output iterator
        /// @param f the function to map-filter over
        /// @return the pair of output iterators by value
        template<SizedRange I1, SizedRange I2, ResizableRange O1, ResizableRange O2,
                 EnumMap2Filter2Functor<I1, I2, O1, O2> F>
        auto enumerate_map2_filter2(IExecutor& pool, I1&& it1, I2&& it2, O1&& out1, O2&& out2, F&& f)
                -> std::pair<Forward<O1>, Forward<O2>>
        {
            // strategy: we allocate out to the same size as it, every thread has its own counter, and when they are all
            // done, we perform a serial memcpy. memcpy is memory-bound, meaning there is not a lot to gain from
            // parallelization.
            const auto pool_size = pool.size();
            const auto work_size = ParallelUtil::smallest_size(it1, it2);
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);
            if (work_size == 0) {
                return ParallelUtil::make_pair_wrapper(std::forward<O1&&>(out1), std::forward<O2&&>(out2));
            }
            if (out1.size() != work_size) {
                out1.reserve(work_size);
                // safety post-condition: we need to manually resize out to the correct size at the end
                out1.resize(work_size);
            }
            if (out2.size() != work_size) {
                out2.reserve(work_size);
                // safety post-condition: we need to manually resize out to the correct size at the end
                out2.resize(work_size);
            }
            using ParallelUtil::CacheLineCounter;
            std::vector<CacheLineCounter> counters(pool_size, CacheLineCounter(0));
            for (size_t thread_id = 0; thread_id < pool_size; ++thread_id) {
                auto& thread_counter = counters[thread_id];
                thread_counter.inner = 0;
                pool.addTask([&out1, thread_id, reduced_size, work_size, &it1, &f, &thread_counter, &it2, &out2] {
                    auto max_size = std::min((thread_id + 1) * reduced_size, work_size);
                    auto begin = thread_id * reduced_size;
                    for (auto j = begin; j < max_size; ++j) {
                        const auto& input_value1 = it1[j];
                        const auto& input_value2 = it2[j];
                        auto result = f(j, input_value1, input_value2);
                        if (result.has_value()) {
                            auto [inner1, inner2] = std::move(result.value());
                            out1[begin + thread_counter.inner] = inner1;
                            out2[begin + thread_counter.inner] = inner2;
                            thread_counter.inner++;
                        }
                    }
                });
            }
            pool.waitAll();

            ParallelUtil::fix_chunks(out1, counters, reduced_size);
            ParallelUtil::fix_chunks(out2, counters, reduced_size);

            return ParallelUtil::make_pair_wrapper(std::forward<O1>(out1), std::forward<O2>(out2));
        }

        namespace detail {
            /// Serial fold
            template<typename State, typename InputIt, FoldFunctor<State, InputIt> F>
            auto fold(InputIt&& it, State&& initial, F&& f) -> State
            {
                State next = initial;
                for (size_t i = 0; i < it.size(); ++i) { next = f(next, it[i]); }
                return next;
            }
        } // namespace detail

        // 'state -> 'a collection -> ('state -> 'a -> 'state) -> 'state
        /// @brief Perform a parallel left fold
        /// @details
        /// @code Performs f(..., f(I[2], f([I[1], f(I[0], initial)))...) @endcode
        /// @tparam F a binary function of type State -> Item I -> State
        /// @tparam State the state type
        /// @tparam I Input range
        /// @param pool The threadpool to use
        /// @param it The input iterator
        /// @param initial Initial value of the state
        /// @param identity The left identity with respect to f. For example the left identity with respect to addition
        /// is the 0 matrix and the left identity with respect to multiplication is the identity matrix.
        /// @param f hello
        /// @return Folded state
        template<typename State, SizedRange I, FoldFunctor<State, I> F>
        auto fold(IExecutor& pool, I&& it, State&& initial, State&& identity, F&& f) -> State
        {
            const auto pool_size = pool.size();
            const auto work_size = it.size();
            const auto reduced_size = ParallelUtil::div_ceil(work_size, pool_size);

            std::vector<State> final(pool_size, identity);

            for (size_t tid = 0; tid < pool_size; ++tid) {
                auto max_size = std::min((tid + 1) * reduced_size, work_size);
                auto begin = tid * reduced_size;

                pool.addTask([&] {
                    State&& next = identity;
                    for (size_t i = begin; i < max_size; ++i) { next = f(next, it[i]); }
                    final[tid] = next;
                });
            }
            pool.waitAll();
            return fold(final, initial, std::forward<F>(f));
        }
    } // namespace Parallel
} // namespace Util

#endif