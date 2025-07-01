//
// Created by Cem Akarsubasi on 5/28/25.
//
// Assertion macros

#ifndef GEL_ASSERT_H
#define GEL_ASSERT_H

#include <sstream>
#include <stdexcept>
#include <format>
#include <iterator>

/// The assert macro in <cassert> is inert in release mode. Certain graph algorithms
/// are very slow necessitating testing on release mode which is why these macros are
/// provided.
namespace Util::Assert::detail
{
template <typename T>
concept has_to_str_operator =
    requires(T t, std::stringstream& os)
    {
        { os << t };
    };

template <typename Left, typename Right>
concept has_eq =
    requires(Left left, Right right)
    {
        { left == right } -> std::same_as<bool>;
        { left != right } -> std::same_as<bool>;
    };

template <typename T>
std::string string_repr(T&& t)
    requires has_to_str_operator<T>
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

template <typename T>
std::string string_repr(T&& t)
{
    return "?";
}

template <typename... Args>
[[noreturn]]
#ifdef _MSC_VER
__declspec(noinline)
#else
__attribute__((noinline))
#endif
void gel_assert_failure(const char* file, const int line, auto func, const char* expr,
                        std::format_string<Args...>&& fmt,
                        Args&&... args)
{
    std::ostream_iterator<char> os(std::cerr);
    std::format_to(os, "Assertion failed at {}:{} in `{}`\nBecause '{}' is false\n", file, line, func, expr);
    std::format_to(os, fmt, std::forward<Args>(args)...);
    std::format_to(os, "\n");
    throw std::runtime_error("GEL assertion failure");
}

template <typename... Args>
[[noreturn]]
#ifdef _MSC_VER
__declspec(noinline)
#else
__attribute__((noinline))
#endif
void gel_assert_false_failure(const char* file, const int line, auto func, const char* expr,
                              std::format_string<Args...>&& fmt, Args&&... args)
{
    std::ostream_iterator<char> os(std::cerr);
    std::format_to(os, "Assertion failed at {}:{} in `{}`\nBecause '{}' is true\n", file, line, func, expr);
    std::format_to(os, fmt, std::forward<Args>(args)...);
    std::format_to(os, "\n");
    throw std::runtime_error("GEL assertion failure");
}

template <typename... Args>
[[noreturn]]
#ifdef _MSC_VER
__declspec(noinline)
#else
__attribute__((noinline))
#endif
void gel_assert_eq_failure(
    const char* file,
    const int line,
    const char* func,
    const char* left,
    const char* right,
    auto left_value,
    auto right_value,
    std::format_string<Args...>&& fmt,
    Args&&... args)
{
    std::ostream_iterator<char> os(std::cerr);
    std::format_to(os, "Assertion failed at {}:{} in `{}`\nBecause '{}' {{{{{}}}}} != '{}' {{{{{}}}}}\n",
                   file,
                   line,
                   func,
                   left,
                   detail::string_repr(left_value).c_str(), right, detail::string_repr(right_value).c_str());

    std::format_to(os, fmt, std::forward<Args>(args)...);
    std::format_to(os, "\n");
    throw std::runtime_error("GEL assertion failure");
}

template <typename... Args>
[[noreturn]]
#ifdef _MSC_VER
__declspec(noinline)
#else
__attribute__((noinline))
#endif
void gel_assert_neq_failure(
    const char* file,
    const int line,
    const char* func,
    const char* left,
    const char* right,
    auto left_value,
    auto right_value,
    std::format_string<Args...>&& fmt,
    Args&&... args)
{
    std::ostream_iterator<char> os(std::cerr);
    std::format_to(os, "Assertion failed at {}:{} in `{}`\nBecause '{}' {{{{{}}}}} == '{}' {{{{{}}}}}\n",
                   file,
                   line,
                   func,
                   left,
                   detail::string_repr(left_value).c_str(), right, detail::string_repr(right_value).c_str());
    std::format_to(os, fmt, std::forward<Args>(args)...);
    std::format_to(os, "\n");
    throw std::runtime_error("GEL assertion failure");
}

template <typename Predicate, typename... Args>
auto gel_assert_impl(Predicate&& p, auto expr, auto file, auto line, auto func, std::format_string<Args...>&& fmt = "",
                     Args&&... args) -> void
{
    [[unlikely]]
    if (!p) {
        gel_assert_failure(file, line, func, expr, std::forward<decltype(fmt)>(fmt), std::forward<Args>(args)...);
    }
}

template <typename Predicate, typename... Args>
auto gel_assert_false_impl(Predicate&& p, auto expr, auto file, auto line, auto func,
                           std::format_string<Args...>&& fmt = "",
                           Args&&... args) -> void
{
    [[unlikely]]
    if (p) {
        gel_assert_false_failure(file, line, func, expr, std::forward<decltype(fmt)>(fmt), std::forward<Args>(args)...);
    }
}

template <typename Left, typename Right, typename... Args>
auto gel_assert_eq_impl(Left&& l, Right&& r, auto expr_l, auto expr_r, auto file, auto line, auto func,
                        std::format_string<Args...>&& fmt = "",
                        Args&&... args) -> void requires detail::has_eq<Right, Left>
{
    [[unlikely]]
    if (l != r) {
        gel_assert_eq_failure(file, line, func, expr_l, expr_r, l, r, std::forward<decltype(fmt)>(fmt),
                              std::forward<Args>(args)...);
    }
}

template <typename Left, typename Right, typename... Args>
auto gel_assert_neq_impl(Left&& l, Right&& r, auto expr_l, auto expr_r, auto file, auto line, auto func,
                         std::format_string<Args...>&& fmt = "",
                         Args&&... args) -> void requires detail::has_eq<Right, Left>
{
    [[unlikely]]
    if (l == r) {
        gel_assert_neq_failure(file, line, func, expr_l, expr_r, l, r, std::forward<decltype(fmt)>(fmt),
                               std::forward<Args>(args)...);
    }
}
}

/// @brief Assertion macro
/// @details If the predicate is false, throws an std::runtime_error. Can be given a format string.
/// @code
/// auto i = 1;
/// auto j = 1;
/// GEL_ASSERT(i == j, "math is broken because {} is not equal to itself", i);
/// @endcode
#define GEL_ASSERT(pred, ...) do { ::Util::Assert::detail::gel_assert_impl((pred), #pred, __FILE__, __LINE__, __func__ __VA_OPT__(,) __VA_ARGS__); } while(false)

/// @brief Assertion macro
/// @details If the predicate is true, throws an std::runtime_error. Can be given a format string.
/// @code
/// auto i = 1;
/// auto j = 2;
/// GEL_ASSERT_FALSE(i == j, "math is broken because {} is equal to {}", i, j);
/// @endcode
#define GEL_ASSERT_FALSE(pred, ...) do { ::Util::Assert::detail::gel_assert_false_impl((pred), #pred, __FILE__, __LINE__, __func__ __VA_OPT__(,) __VA_ARGS__); } while(false)


/// @brief Assertion macro
/// @details If left is not equal to right, throws an std::runtime_error. Can be given a format string.
/// @code
/// auto i = 1;
/// auto j = 1;
/// GEL_ASSERT_EQ(i, j, "math is broken");
/// @endcode
#define GEL_ASSERT_EQ(left, right, ...) do { ::Util::Assert::detail::gel_assert_eq_impl((left), (right), #left, #right, __FILE__, __LINE__, __func__ __VA_OPT__(,) __VA_ARGS__); } while(false)

/// @brief Assertion macro
/// @details If left is equal to right, throws an std::runtime_error. Can be given a format string.
/// @code
/// auto i = 1;
/// auto j = 2;
/// GEL_ASSERT_NEQ(1, 2, "math is broken");
/// @endcode
#define GEL_ASSERT_NEQ(left, right, ...) do { ::Util::Assert::detail::gel_assert_neq_impl((left), (right), #left, #right, __FILE__, __LINE__, __func__ __VA_OPT__(,) __VA_ARGS__); } while(false)

#endif //GEL_ASSERT_H
