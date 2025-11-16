//
// Created by Cem Akarsubasi on 5/13/25.
//

#ifndef INPLACE_VECTOR_H
#define INPLACE_VECTOR_H
#include <iterator>
#include <ranges>

#include <GEL/Util/Assert.h>

// Should be a private type
namespace Util::detail {

    /// Variable-sized container that stores its data vector inside it
template<typename T, size_t N = 8>
class InplaceVector {
public:
    using value_type = T;
    using size_type = size_t;
    using reference = T&;
    using const_reference = const T&;
    using pointer = T*;
    using const_pointer = const T*;
    using iterator = T*;
    using const_iterator = const T*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

private:
    size_type m_size = 0;
    T m_data[N];
public:

    // constructors and destructors

    constexpr InplaceVector() noexcept: m_size(0)
    {
        // Always false
        for (size_type i = 0; i < m_size; ++i) {
            new (&m_data[i]) T();
        }
    }

    explicit constexpr InplaceVector(const size_type size) : m_size(0)
    {
        GEL_ASSERT(size <= N);
        for (size_type i = 0; i < size; ++i) {
            new (&m_data[i]) T();
        }
    }

    explicit constexpr InplaceVector(const size_type size, const T& value) : m_size(0)
    {
        GEL_ASSERT(size <= N);
        for (size_type i = 0; i < size; ++i) {
            new (&m_data[i]) T(value);
        }
    }

    template <std::input_iterator InputIt>
    constexpr InplaceVector(InputIt first, InputIt last)
    {
        for (; first != last; ++first) {
            push_back(*first);
        }
    }

    constexpr InplaceVector(const InplaceVector& other) {
        for (size_t i = 0; i < other.m_size; ++i) {
            m_data[i] = other.m_data[i];
        }
        m_size = other.m_size;

    }

    constexpr InplaceVector(InplaceVector&& other) noexcept
    {
        m_size = other.m_size;
        for (size_t i = 0; i < other.m_size; ++i) {
            m_data[i] = std::move(other.m_data[i]);
        }
        other.m_size = 0;
    }

    constexpr ~InplaceVector()
    {
        if constexpr (!std::is_trivially_destructible_v<T>) {
            for (size_type i = 0; i < m_size; ++i) {
                m_data[i].~T();
            }
        }
    }

    constexpr InplaceVector& operator=(const InplaceVector& other)
    {
        if (this != &other) {
            for (size_type i = 0; i < other.m_size; ++i) {
                m_data[i] = other.m_data[i];
            }
            m_size = other.m_size;
        }
        return *this;
    }

    constexpr InplaceVector& operator=(InplaceVector&& other) noexcept
    {
        if (this != &other) {
            for (size_type i = 0; i < other.m_size; ++i) {
                m_data[i] = std::move(other.m_data[i]);
            }
            m_size = other.m_size;
            other.m_size = 0;
        }
        return *this;
    }

    // Element access

    constexpr reference operator[](size_t index) {
        return *(m_data+index);
    }

    constexpr const_reference operator[](size_t index) const {
        return *(m_data+index);
    }

    constexpr reference at(size_t index) {
        GEL_ASSERT(index < m_size);
        return *(m_data+index);
    }

    [[nodiscard]] constexpr const_reference at(size_t index) const {
        GEL_ASSERT(index < m_size);
        return *(m_data+index);
    }

    constexpr reference front()
    {
        return *m_data;
    }

    [[nodiscard]] constexpr const_reference front() const
    {
        return *m_data;
    }

    constexpr reference back()
    {
        return *(m_data+m_size-1);
    }

    [[nodiscard]] constexpr const_reference back() const
    {
        return *(m_data+m_size-1);
    }

    constexpr pointer data() noexcept
    {
        return m_data;
    }

    [[nodiscard]] constexpr const_pointer data() const noexcept
    {
        return m_data;
    }

    // Iterators

    constexpr iterator begin()
    {
        return m_data;
    }

    constexpr iterator end()
    {
        return m_data + m_size;
    }

    [[nodiscard]]
    constexpr const_iterator begin() const
    {
        return m_data;
    }

    [[nodiscard]]
    constexpr const_iterator end() const
    {
        return m_data + m_size;
    }

    [[nodiscard]] constexpr const_iterator cbegin() const
    {
        return m_data;
    }

    [[nodiscard]] constexpr const_iterator cend() const
    {
        return m_data + m_size;
    }

    constexpr reverse_iterator rbegin()
    {
        return std::reverse_iterator(end());
    }

    constexpr reverse_iterator rend()
    {
        return std::reverse_iterator(begin());
    }

    [[nodiscard]] constexpr const_reverse_iterator crbegin() const
    {
        return std::reverse_iterator(cend());
    }

    [[nodiscard]] constexpr const_reverse_iterator crend() const
    {
        return std::reverse_iterator(cbegin());
    }

    // Capacity

    [[nodiscard]] constexpr bool empty() const
    {
        return (m_size == 0);
    }

    [[nodiscard]] constexpr size_type size() const
    {
        return m_size;
    }

    [[nodiscard]] constexpr size_type max_size() const
    {
        return N;
    }

    constexpr void reserve(const size_type new_cap)
    {
        // we literally can't do anything
        GEL_ASSERT(new_cap <= N);
    }

    [[nodiscard]] constexpr size_type capacity() const
    {
        return N;
    }

    // Modifiers

    constexpr void clear()
    {
        for (size_type i = 0; i < m_size; ++i) {
            m_data[i].~T();
        }
        m_size = 0;
    }

    constexpr void push_back(const T& value)
    {
        GEL_ASSERT(m_size < N);
        new (&m_data[m_size++]) T(value);
        // m_data[m_size++] = value;
    }


    constexpr T pop_back()
    {
        GEL_ASSERT(m_size > 0);
        return std::move(m_data[--m_size]);
    }

    template < class... Args >
    constexpr void emplace_back(Args&&... args)
    {
        GEL_ASSERT(m_size < N);
        new (m_data + m_size) T(std::forward<Args>(args)...);
        m_size++;
    }
};
    static_assert(std::ranges::viewable_range<InplaceVector<int>&>);

}

#endif //INPLACE_VECTOR_H
