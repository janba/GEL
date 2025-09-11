//
// Created by Cem Akarsubasi on 9/9/25.
//
// Useful utilities for working with ranges

#ifndef GEL_RANGETOOLS_H
#define GEL_RANGETOOLS_H

#include <ranges>

// TODO: Cleaner encapsulation of types included here

/// Useful utilities for working with ranges
///
/// @details C++20 ranges library is unfortunately fairly incomplete, lacking a large chunk of useful view adapters
/// that are needed for functional style code.
namespace Util::Ranges
{
template <std::ranges::input_range R1, std::ranges::input_range R2>
struct zip_view_t : std::ranges::view_interface<zip_view_t<R1, R2>> {
    using value_type = std::pair<std::ranges::range_value_t<R1>, std::ranges::range_value_t<R2>>;
    using reference = std::pair<std::ranges::range_reference_t<R1>, std::ranges::range_reference_t<R2>>;

private:
public:
    class sentinel_t;
    std::ranges::iterator_t<R1> m_begin1;
    std::ranges::iterator_t<R2> m_begin2;
    std::ranges::sentinel_t<R1> m_end1;
    std::ranges::sentinel_t<R2> m_end2;

    class iterator_t {
    public:
        std::ranges::iterator_t<R1> ptr1;
        std::ranges::iterator_t<R2> ptr2;

    public:
        using difference_type = ptrdiff_t;
        using value_type = std::pair<typename decltype(ptr1)::value_type, typename decltype(ptr2)::value_type>;
        constexpr iterator_t() = default;
        constexpr explicit iterator_t(auto _ptr1, auto _ptr2) : ptr1(_ptr1), ptr2(_ptr2) {}

        iterator_t& operator++()
        {
            ++ptr1;
            ++ptr2;
            return *this;
        }

        iterator_t operator++(int)
        {
            iterator_t retval = *this;
            ++(*this);
            return retval;
        }

        bool operator==(iterator_t other) const
        {
            return ptr1 == other.ptr1 || ptr2 == other.ptr2;
        }

        bool operator==(sentinel_t other) const
        {
            return ptr1 == other.m_end1 || ptr2 == other.m_end2;
        }

        bool operator!=(iterator_t other) const { return !(*this == other); }

        reference operator*() const
        {
            return {*ptr1, *ptr2};
        }
    };

    class sentinel_t {
    public:
        std::ranges::sentinel_t<R1> m_end1;
        std::ranges::sentinel_t<R2> m_end2;
        friend class iterator_t;
    };

public:
    zip_view_t() = default;

    zip_view_t(R1&& range1, R2&& range2) :
        m_begin1(std::ranges::begin(range1)),
        m_begin2(std::ranges::begin(range2)),
        m_end1(std::ranges::end(range1)),
        m_end2(std::ranges::end(range2)) {}

    [[nodiscard]]
    auto begin() const -> iterator_t
    {
        return iterator_t(m_begin1, m_begin2);
    }

    [[nodiscard]]
    auto end() const -> sentinel_t
    {
        return sentinel_t(m_end1, m_end2);
    }
};

namespace detail
{
    using ExampleIter = std::vector<int>;
    static_assert(std::forward_iterator<zip_view_t<ExampleIter, ExampleIter>::iterator_t>);
    static_assert(
        std::sentinel_for<zip_view_t<ExampleIter, ExampleIter>::sentinel_t, zip_view_t<
                              ExampleIter, ExampleIter>::iterator_t>);
}

/// Helper function to create a zip_view_t
///
/// @related zip_view_t
struct Zip /*: : std::ranges::range_adaptor_closure<Zip> */ {
    template <std::ranges::input_range R1, std::ranges::input_range R2>
    [[nodiscard]]
    constexpr auto operator()(R1&& range1, R2&& range2) const -> zip_view_t<R1, R2>
    {
        return zip_view_t<R1, R2>{std::forward<R1>(range1), std::forward<R2>(range2)};
    }
};
// FIXME: replace with std::ranges::zip when upgrading to C++23
static constexpr Zip zip;

/// Cartesian product of two ranges
template <std::ranges::input_range R1, std::ranges::input_range R2>
std::ranges::input_range auto cartesian_product(R1&& range1, R2&& range2)
{
    //using T1 = std::ranges::range_value_t<R1>;
    //using T2 = std::ranges::range_value_t<R2>;

    auto rv1 = range1 | std::views::all;
    auto rv2 = range2 | std::views::all;

    auto size1 = std::ranges::size(rv1);
    auto size2 = std::ranges::size(rv2);

    auto wide1 = rv1 | std::views::transform([size2](auto e1) {
        std::views::iota(0, size2) | std::views::transform([e1](auto _i) {
            return e1;
        });
    }) | std::views::join;

    auto wide2 = std::views::iota(0, size1) | std::views::transform([rv2](auto _i) {
        return rv2;
    }) | std::views::join;

    //return wide1 | zip(wide2);
    return zip(wide1, wide2);
}

template <std::ranges::input_range R>
std::ranges::input_range auto repeat_range(R&& range)
{
    return std::views::iota(0ULL)
    | std::views::transform([range = std::forward<R>(range)](auto _unused) { return range; })
    | std::views::join;
}

template <std::copy_constructible V>
std::ranges::input_range auto repeat(V&& v)
{
    auto range = std::views::iota(0ULL) | std::views::transform([v = std::forward<V>(v)](auto _unused)  { return v; });
    return range;
}

//template <std::ranges::input_range R>
std::ranges::input_range auto shifted_wrapping(std::ranges::input_range auto&& range, size_t m)
{
    auto ranges = std::views::iota(0ULL) | std::views::transform([&](auto _unused) { return std::views::all(range); });
    auto joined = std::views::join(ranges);
    auto n = std::ranges::distance(range);
    m %= n;
    return joined | std::views::drop(m) | std::views::take(n);
}

template <std::ranges::viewable_range Range>
std::ranges::viewable_range auto slide2_wraparound(Range&& range)
{
    // TODO: figure this out
    //auto length = std::ranges::distance(range) + 1;
    //auto repeating1 = repeat_range(range);
    //auto repeating2 = repeat_range(range);
    //return zip(repeating1 | std::views::take(length), repeating2 | std::views::drop(1) | std::views::take(length));
    using Elem = std::remove_cvref_t<std::ranges::range_value_t<Range>>;
    const auto size = range.size();
    auto ptr = range.begin();
    const auto end = range.end();
    return std::views::iota(0UL, (size > 1) ? size : 0UL)
    | std::views::transform([=, &range]([[maybe_unused]] const auto& _i) mutable -> std::tuple<Elem const&, Elem const&> {
        const auto& current = ptr;
        ++ptr;
        const auto& next = (ptr == end) ? range.begin() : ptr;
        return std::tie(*current, *next);
    });
}

template<std::ranges::random_access_range Range, std::ranges::input_range Selector>
auto select_lazy(const Range& range, Selector&& selector)
{
    return selector | std::views::transform([&range](auto idx) {return range[idx];});
}

template <typename G, typename I>
concept GeneratorClosure =
    std::copy_constructible<I> && // Need to return I by value
    std::copy_constructible<G> &&
    requires(G g)
    {
        { g() } -> std::same_as<std::optional<I>>;
    };

/// An atrocious clone of Rust's Iterator interface
//template <std::copy_constructible Item, GeneratorClosure<Item> G>
template <typename Item, GeneratorClosure<Item> Generator>
struct Iterator {
    // The generator has to be stored so if the Iterator is reused, we can copy construct it.
    Generator generator;
    explicit Iterator(Generator&& g = Generator()) : generator(std::forward<Generator>(g)) {}

    struct sentinel_t;

    struct iterator_t {
        Generator generator;
        // iterator equality does not make much sense for such a generator, but C++ requires us to track this
        // meaning we have to store an index
        size_t index = 0;
        std::optional<Item> item = std::nullopt;
        using difference_type = ptrdiff_t;
        using value_type = Item;
        iterator_t() = default;
        iterator_t(const iterator_t& rhs) : generator(rhs.generator), index(rhs.index), item(rhs.item) {}
        iterator_t& operator=(const iterator_t& rhs) = default;
        iterator_t(iterator_t&& rhs) noexcept = default; //: generator(rhs.generator), index(rhs.index), item(rhs.item) {}
        iterator_t& operator=(iterator_t&& rhs) = default;
        explicit iterator_t(const Generator& generator) : generator(generator), item(std::nullopt)
        {
            // Since advancing the iterator and fetching elements are different operations
            // we need to initialize the item
            item = this->generator();
        }

        iterator_t& operator++()
        {
            index++;
            item = generator();
            return *this;
        }

        iterator_t operator++(int)
        {
            iterator_t retval = *this;
            ++(*this);
            return retval;
        }

        bool operator==(iterator_t other) const
        {
            return index == other.index;
        }

        bool operator==([[maybe_unused]] sentinel_t other) const
        {
            return !item.has_value();
        }

        bool operator!=(iterator_t other) const { return !(*this == other); }

        value_type operator*() const
        {
            return *item;
        }
    };
    // A sentinel can be entirely empty, since whether we hit the end depends on whether the closure returns further
    // items
    struct sentinel_t {};

    auto begin() -> iterator_t
    {
        return iterator_t(generator);
    }

    auto end() -> sentinel_t
    {
        return sentinel_t{};
    }

    static_assert(std::input_iterator<iterator_t>);
    static_assert(std::sentinel_for<sentinel_t, iterator_t>);
};

namespace detail
{
    struct DummyGenerator {
        std::optional<int> operator()()
        {
            return std::nullopt;
        }
    };
    struct FibonacciGenerator {
        std::uint64_t n0 = 0;
        std::uint64_t n1 = 1;

        std::optional<std::uint64_t> operator()()
        {
            const auto temp = n0 + n1;
            n0 = n1;
            n1 = temp;
            return n0 + n1;
        }
    };
    inline void test()
    {
        auto it = Iterator<std::uint64_t, FibonacciGenerator>(FibonacciGenerator());
        for (auto fib : it | std::views::take(50)) {
            std::cout << fib << '\n';
        }
    }
    using ExampleIter = std::vector<int>;
    static_assert(std::forward_iterator<Iterator<int, DummyGenerator>::iterator_t>);
    static_assert(
        std::sentinel_for<Iterator<int, DummyGenerator>::sentinel_t, Iterator<int, DummyGenerator>::iterator_t>);
}

/// Create an infinite generator using the provided initial state and fold function (State -> State)
template <std::copy_constructible State, typename Functor> requires std::same_as<std::invoke_result_t<Functor, State>, std::remove_cvref_t<State>>
std::ranges::viewable_range auto make_generator(State&& initial, Functor&& func)
{
    // For some reason, without the explicit remove_cvref_t, this closure returns a reference to "current" rather
    // than returning it by value
    auto generator = [state = initial, &func]([[maybe_unused]] int _unused) mutable -> std::remove_cvref_t<State> {
        auto current = state;
        state = func(current);
        return current;
    };
    return std::views::iota(0) | std::views::transform(generator);
}

/// Create an infinite generator using the provided initial state and fold function (State -> optional State)
// template <std::copy_constructible State, typename Functor> //requires std::same_as<std::invoke_result_t<Functor, State>, std::remove_cvref_t<State>>
// std::ranges::range auto make_generator_finite(State&& initial, Functor&& func)
// {
//     // using StateVal = std::remove_cvref_t<State>;
//     // struct Generator {
//     //     StateVal initial;
//     //     StateVal current;
//     //
//     //     std::optional<std::remove_cvref_t<State>> operator()()
//     //     {
//     //         if (current.has_value())
//     //             current = func(current.value());
//     //         return current;
//     //     }
//     // };
//     // auto generator = [state = initial, &func]() mutable -> std::optional<std::remove_cvref_t<State>> {
//     //     std::optional<std::remove_cvref_t<State>> current = state;
//     //     if (current.has_value())
//     //         state = func(current.value());
//     //     return current;
//     // };
//     // return Iterator<std::remove_cvref_t<State>, decltype(generator)>(std::move(generator));
// }

// struct iterator_t {
//     using difference_type = ptrdiff_t;
//     using value_type = Item;
//     iterator_t() = default;
//     iterator_t(const iterator_t& rhs) = default;
//     iterator_t& operator=(const iterator_t& rhs) = default;
//     iterator_t(iterator_t&& rhs) noexcept = default;
//     iterator_t& operator=(iterator_t&& rhs) = default;
//     explicit iterator_t(...): ... {}
//
//     iterator_t& operator++()
//     {
//         index++;
//         return *this;
//     }
//
//     iterator_t operator++(int)
//     {
//         iterator_t retval = *this;
//         ++(*this);
//         return retval;
//     }
//
//     bool operator==(iterator_t other) const
//     {
//         return /**/;
//     }
// };

}


#endif //GEL_RANGETOOLS_H
