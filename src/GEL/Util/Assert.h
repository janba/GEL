//
// Created by Cem Akarsubasi on 5/28/25.
//
// Assertion macros

#ifndef GEL_ASSERT_H
#define GEL_ASSERT_H

#include <cstdio>
#include <sstream>
#include <stdexcept>

/// The assert macro in <cassert> is inert in release mode. Certain graph algorithms
/// are very slow necessitating testing on release mode which is why these macros are
/// provided.
namespace Util::Assert::detail {
    template<typename T> concept has_to_str_operator = requires(T t, std::stringstream& os) {
        { os << t };
    };

    template<typename Left, typename Right> concept has_eq = requires(Left left, Right right) {
        { left == right } -> std::same_as<bool>;
        { left != right } -> std::same_as<bool>;
    };

    template<typename T> std::string string_repr(T&& t)
    requires has_to_str_operator<T>
    {
        std::stringstream ss;
        ss << t;
        return ss.str();
    }

    template<typename T> std::string string_repr(T&& t)
    {
        return "?";
    }

#ifdef __GNUG__
#define IGNORE_WFORMAT _Pragma("GCC diagnostic ignored \"-Wformat-security\"")
#else
#define IGNORE_WFORMAT
#endif

#ifdef _MSC_VER
    __declspec(noinline)
#else
    __attribute__((noinline))
#endif
    void
    gel_assert_failure(const char* file, const int line, auto func, const char* expr, const char* fmt, auto... args)
    {
        fprintf(stderr, "Assertion failed at %s:%d in `%s`\nBecause '%s' is false\n", file, line, func, expr);
        IGNORE_WFORMAT
        fprintf(stderr, fmt, args...);
        fprintf(stderr, "\n");
        throw std::runtime_error("GEL assertion failure");
    }

#ifdef _MSC_VER
    __declspec(noinline)
#else
    __attribute__((noinline))
#endif
    void
    gel_assert_false_failure(const char* file, const int line, auto func, const char* expr, const char* fmt,
                             auto... args)
    {
        fprintf(stderr, "Assertion failed at %s:%d in `%s`\nBecause '%s' is true\n", file, line, func, expr);
        IGNORE_WFORMAT
        fprintf(stderr, fmt, args...);
        fprintf(stderr, "\n");
        throw std::runtime_error("GEL assertion failure");
    }

#ifdef _MSC_VER
    __declspec(noinline)
#else
    __attribute__((noinline))
#endif
    void
    gel_assert_eq_failure(const char* file, const int line, const char* func, const char* left, const char* right,
                          auto left_value, auto right_value, const char* fmt, auto... args)
    {
        fprintf(stderr, "Assertion failed at %s:%d in `%s`\nBecause '%s' {{%s}} != '%s' {{%s}}\n", file, line, func,
                left, string_repr(left_value).c_str(), right, string_repr(right_value).c_str());
        IGNORE_WFORMAT
        fprintf(stderr, fmt, args...);
        fprintf(stderr, "\n");
        throw std::runtime_error("GEL assertion failure");
    }

#ifdef _MSC_VER
    __declspec(noinline)
#else
    __attribute__((noinline))
#endif
    void
    gel_assert_neq_failure(const char* file, const int line, const char* func, const char* left, const char* right,
                           auto left_value, auto right_value, const char* fmt, auto... args)
    {
        fprintf(stderr, "Assertion failed at %s:%d in `%s`\nBecause '%s' {{%s}} == '%s' {{%s}}\n", file, line, func,
                left, string_repr(left_value).c_str(), right, string_repr(right_value).c_str());
        IGNORE_WFORMAT
        fprintf(stderr, fmt, args...);
        fprintf(stderr, "\n");
        throw std::runtime_error("GEL assertion failure");
    }

    template<typename Predicate, typename... Args>
    auto gel_assert_impl(Predicate&& p, auto expr, auto file, auto line, auto func, Args... args) -> void
    {
        [[unlikely]]
        if (!p) {
            gel_assert_failure(file, line, func, expr, args...);
        }
    }

    template<typename Predicate, typename... Args>
    auto gel_assert_false_impl(Predicate&& p, auto expr, auto file, auto line, auto func, Args... args) -> void
    {
        [[unlikely]]
        if (p) {
            gel_assert_false_failure(file, line, func, expr, args...);
        }
    }

    template<typename Left, typename Right, typename... Args>
    auto gel_assert_eq_impl(Left&& l, Right&& r, auto expr_l, auto expr_r, auto file, auto line, auto func,
                            Args... args) -> void
    requires has_eq<Right, Left>
    {
        [[unlikely]]
        if (l != r) {
            gel_assert_eq_failure(file, line, func, expr_l, expr_r, l, r, args...);
        }
    }

    template<typename Left, typename Right, typename... Args>
    auto gel_assert_neq_impl(Left&& l, Right&& r, auto expr_l, auto expr_r, auto file, auto line, auto func,
                             Args... args) -> void
    requires has_eq<Right, Left>
    {
        [[unlikely]]
        if (l == r) {
            gel_assert_neq_failure(file, line, func, expr_l, expr_r, l, r, args...);
        }
    }
} // namespace Util::Assert::detail

/// @brief Assertion macro
/// @details If pred is false, throws an std::runtime_error. Can be given a C-style format string.
/// @code
/// auto i = 1;
/// auto j = 2;
/// GEL_ASSERT_FALSE(i == j, "math is broken because %d is equal to %d", i, j);
/// @endcode
#define GEL_ASSERT(pred, ...)                                                                                          \
    do {                                                                                                               \
        ::Util::Assert::detail::gel_assert_impl((pred), #pred, __FILE__, __LINE__, __func__, "" __VA_ARGS__);          \
    }                                                                                                                  \
    while (false)

/// @brief Assertion macro
/// @details If pred is true, throws an std::runtime_error. Can be given a C-style format string.
/// @code
/// GEL_ASSERT(1 == 1, "math is broken");
/// @endcode
#define GEL_ASSERT_FALSE(pred, ...)                                                                                    \
    do {                                                                                                               \
        ::Util::Assert::detail::gel_assert_false_impl((pred), #pred, __FILE__, __LINE__, __func__, "" __VA_ARGS__);    \
    }                                                                                                                  \
    while (false)

/// @brief Assertion macro
/// @details If left is not equal to right, throws an std::runtime_error. Can be given a C-style format string.
/// @code
/// auto i = 1;
/// auto j = 1;
/// GEL_ASSERT_EQ(i, j, "math is broken");
/// @endcode
#define GEL_ASSERT_EQ(left, right, ...)                                                                                \
    do {                                                                                                               \
        ::Util::Assert::detail::gel_assert_eq_impl((left), (right), #left, #right, __FILE__, __LINE__, __func__,       \
                                                   "" __VA_ARGS__);                                                    \
    }                                                                                                                  \
    while (false)

/// @brief Assertion macro
/// @details If left is equal to right, throws an std::runtime_error. Can be given a C-style format string.
/// @code
/// auto i = 1;
/// auto j = 2;
/// GEL_ASSERT_NEQ(1, 2, "math is broken");
/// @endcode
#define GEL_ASSERT_NEQ(left, right, ...)                                                                               \
    do {                                                                                                               \
        ::Util::Assert::detail::gel_assert_neq_impl((left), (right), #left, #right, __FILE__, __LINE__, __func__,      \
                                                    "" __VA_ARGS__);                                                   \
    }                                                                                                                  \
    while (false)

#endif // GEL_ASSERT_H
