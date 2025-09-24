/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * @brief Abstract vector class
 * ----------------------------------------------------------------------- */

/** @file ArithVec.h
 * @brief Abstract vector class
 */

#ifndef CGLA_ARITHVEC_H
#define CGLA_ARITHVEC_H

#include <array>
#include <algorithm>
#include <iostream>
#include <GEL/CGLA/CGLA-util.h>

#include <numeric>

namespace CGLA
{

    /** \brief Template representing generic arithmetic vectors.

     The three parameters to the template are

     T - the scalar type (i.e. float, int, double etc.)

     V - the name of the vector type. This template is always (and
     only) used as ancestor of concrete types, and the name of the
     class _inheriting_ _from_ this class is used as the V argument.

     N - The final argument is the dimension N. For instance, N=3 for a
     3D vector.

     This class template contains all functions that are assumed to be
     the same for any arithmetic vector - regardless of dimension or
     the type of scalars used for coordinates.

     The template contains no virtual functions which is important
     since they add overhead.
     */

    template <class T, class V, unsigned int N>
    class ArithVec
    {

    protected:

        /// The actual contents of the vector.
        std::array<T,N> data;

    protected:

        /// @name Constructors
        /// @{

        /// Construct uninitialized vector
        constexpr ArithVec() noexcept = default;

        /// Construct a vector where all coordinates are identical
        constexpr explicit ArithVec(T _a) noexcept
        {
            std::fill_n(data, N, _a);
        }

        /// Construct a 2D vector
        constexpr ArithVec(T _a, T _b) noexcept
        {
            static_assert(N==2);
            data[0] = _a;
            data[1] = _b;
        }

        /// Construct a 3D vector
        constexpr ArithVec(T _a, T _b, T _c) noexcept
        {
            static_assert(N==3);
            data[0] = _a;
            data[1] = _b;
            data[2] = _c;
        }

        /// Construct a 4D vector
        constexpr ArithVec(T _a, T _b, T _c, T _d) noexcept
        {
            static_assert(N==4);
            data[0] = _a;
            data[1] = _b;
            data[2] = _c;
            data[3] = _d;
        }

        /// @}

    public:

        /// For convenience we define a more meaningful name for the scalar type
        using ScalarType = T;

        /// A more meaningful name for vector type
        using VectorType = V;

        /// @name Accessors
        /// @{

        /// Return dimension of vector
        static constexpr unsigned int get_dim() {return N;}

        /// Set all coordinates of a 2D vector.
        constexpr void set(T _a, T _b) noexcept
        {
            static_assert(N==2);
            data[0] = _a;
            data[1] = _b;
        }

        /// Set all coordinates of a 3D vector.
        constexpr void set(T _a, T _b, T _c) noexcept
        {
            static_assert(N==3);
            data[0] = _a;
            data[1] = _b;
            data[2] = _c;
        }

        /// Set all coordinates of a 4D vector.
        constexpr void set(T _a, T _b, T _c, T _d) noexcept
        {
            static_assert(N==4);
            data[0] = _a;
            data[1] = _b;
            data[2] = _c;
            data[3] = _d;
        }

        /// Const index operator
        constexpr const T& operator [] ( unsigned int i ) const
        {
            assert(i<N);
            return data[i];
        }

        /// Non-const index operator
        constexpr T& operator [] ( unsigned int i )
        {
            assert(i<N);
            return data[i];
        }

        /// Const index operator
        constexpr const T& operator () ( unsigned int i ) const
        {
            assert(i<N);
            return data[i];
        }

        /// Non-const index operator
        constexpr T& operator () ( unsigned int i )
        {
            assert(i<N);
            return data[i];
        }

        /// @}
        /// @name Iterators
        /// @{

        using iterator = typename std::array<T, N>::iterator;
        using const_iterator = typename std::array<T, N>::const_iterator;


        /** Get a pointer to first element in data array.
         This function may be useful when interfacing with some other API
         such as OpenGL (TM) */
        constexpr iterator begin() {return data.begin();}
        constexpr iterator end() {return data.end();}
        constexpr T* get() {return data.data();}

        /** Get a const pointer to first element in data array.
         This function may be useful when interfacing with some other API
         such as OpenGL (TM). */
        constexpr const_iterator begin() const {return data.begin();}
        constexpr const_iterator end() const {return data.end();}
        constexpr const T* get() const {return data.data();}

        /// @}
        /// @name Comparison Operators
        /// @{

        /// Equality operator
        constexpr bool operator==(const V& v) const
        {
            return std::equal(begin(),end(), v.begin());
        }

        /// Equality wrt scalar. True if all coords are equal to scalar
        constexpr bool operator==(T k) const
        {
            return std::count(begin(),end(), k)==N;
        }

        /// Inequality operator
        constexpr bool operator!=(const V& v) const
        {
            return !(*this==v);
        }

        /// Inequality wrt scalar. True if any coord not equal to scalar
        constexpr bool operator!=(T k) const
        {
            return !(*this==k);
        }


        /// @}
        /// @name Comparison and geometric significance operators
        /// @{

        /** Compare all coordinates against other vector ( < ) and combine with 'and'.
         Similar to testing whether we are on one side of three planes. */
        constexpr bool  all_l  (const V& v) const
        {
            return std::inner_product(begin(), end(), v.begin(), true,
                                      std::logical_and<>(), std::less<T>());
        }

        /** Compare all coordinates against other vector ( <= ) and combine with 'and'.
         Similar to testing whether we are on one side of three planes. */
        constexpr bool  all_le (const V& v) const
        {
            return std::inner_product(begin(), end(), v.begin(), true,
                                      std::logical_and<>(), std::less_equal<T>());
        }

        /** Compare all coordinates against other vector ( > ) and combine with 'and'.
         Similar to testing whether we are on one side of three planes. */
        constexpr bool  all_g  (const V& v) const
        {
            return std::inner_product(begin(), end(), v.begin(), true,
                                      std::logical_and<>(), std::greater<T>());
        }

        /** Compare all coordinates against other vector  ( >= ) and combine with 'and'.
         Similar to testing whether we are on one side of three planes. */
        constexpr bool  all_ge (const V& v) const
        {
            return std::inner_product(begin(), end(), v.begin(), true,
                                      std::logical_and<>(), std::greater_equal<T>());
        }

        /** Compare all coordinates against other vector ( < ) and combine with 'or'.
         Similar to testing whether we are on one side of three planes. */
        constexpr bool  any_l  (const V& v) const
        {
            return std::inner_product(begin(), end(), v.begin(), true,
                                      std::logical_or<>(), std::less<T>());
        }

        /** Compare all coordinates against other vector ( <= ) and combine with 'or'.
         Similar to testing whether we are on one side of three planes. */
        constexpr bool  any_le (const V& v) const
        {
            return std::inner_product(begin(), end(), v.begin(), true,
                                      std::logical_or<>(), std::less_equal<T>());
        }

        /** Compare all coordinates against other vector ( > ) and combine with 'or'.
         Similar to testing whether we are on one side of three planes. */
        constexpr bool  any_g  (const V& v) const
        {
            return std::inner_product(begin(), end(), v.begin(), true,
                                      std::logical_or<>(), std::greater<T>());
        }

        /** Compare all coordinates against other vector  ( >= ) and combine with 'or'.
         Similar to testing whether we are on one side of three planes. */
        constexpr bool  any_ge (const V& v) const
        {
            return std::inner_product(begin(), end(), v.begin(), true,
                                      std::logical_or<bool>(), std::greater_equal<T>());
        }

        /// Returns true if p(x) is true for any vector element. p :: T -> bool
        template <typename Predicate>
        constexpr bool any(Predicate&& p) const
        {
            return std::any_of(begin(), end(), p);
        }

        /// Returns true if p(x) is true for all vector elements. p :: T -> bool
        template <typename Predicate>
        constexpr bool all(Predicate&& p) const
        {
            return std::all_of(begin(), end(), p);
        }

        /// @}
        /// @name Assignment operators
        /// @{

        /// Assignment multiplication with scalar.
        constexpr V& operator *=(T k)
        {
            for(auto& x : data) {x*=k;}
            return static_cast<V&>(*this);
        }

        /// Assignment division with scalar.
        constexpr V& operator /=(T k)
        {
            for(auto& x : data) {x/=k;}

            return static_cast<V&>(*this);
        }

        /// Assignment addition with scalar. Adds scalar to each coordinate.
        constexpr V& operator +=(T k)
        {
            for(auto& x : data) {x+=k;}

            return static_cast<V&>(*this);
        }

        /// Assignment subtraction with scalar. Subtracts scalar from each coord.
        constexpr V& operator -=(T k)
        {
            for(auto& x : data) {x-=k;}
            return static_cast<V&>(*this);
        }

        /** Assignment multiplication with vector.
         Multiply each coord independently. */
        constexpr V& operator *=(const V& v)
        {
            for (auto i = 0; i < N; ++i) {
                data[i] = data[i] * v[i];
            }
            return static_cast<V&>(*this);
        }

        /// Assigment division with vector. Each coord divided independently.
        constexpr V& operator /=(const V& v)
        {
            for (auto i = 0; i < N; ++i) {
                data[i] = data[i] / v[i];
            }
            return static_cast<V&>(*this);
        }

        /// Assignmment addition with vector.
        constexpr V& operator +=(const V& v)
        {
            for (auto i = 0; i < N; ++i) {
                data[i] = data[i] + v[i];
            }
            return static_cast<V&>(*this);
        }

        /// Assignment subtraction with vector.
        constexpr V& operator -=(const V& v)
        {
            for (auto i = 0; i < N; ++i) {
                data[i] = data[i] - v[i];
            }
            return static_cast<V&>(*this);
        }

        /// @}
        /// @name Unary operators
        /// @{

        /// Negate vector.
        constexpr V operator - () const
        {
            V ret;
            for (auto i = 0; i < N; ++i) {
                ret[i] = -data[i];
            }
            return ret;
        }

        /// @}
        /// @name Binary operations with vectors
        /// @{

        /** Multiply vector with vector. Each coord multiplied independently
         Do not confuse this operation with dot product. */
        constexpr V operator * (const V& v1) const
        {
            V ret;
            for (auto i = 0; i < N; ++i) {
                ret[i] = data[i] * v1[i];
            }
            return ret;
        }

        /// Add two vectors
        constexpr V operator + (const V& v1) const
        {
            V ret;
            for (auto i = 0; i < N; ++i) {
                ret[i] = data[i] + v1[i];
            }
            return ret;
        }

        /// Subtract two vectors.
        constexpr V operator - (const V& v1) const
        {
            V ret;
            for (auto i = 0; i < N; ++i) {
                ret[i] = data[i] - v1[i];
            }
            return ret;
        }

        /// Divide two vectors. Each coord separately
        constexpr V operator / (const V& v1) const
        {
            V ret;
            for (auto i = 0; i < N; ++i) {
                ret[i] = data[i] / v1[i];
            }
            return ret;
        }

        /// @}
        /// @name Binary operations with scalars
        /// @{

        /// Multiply scalar onto vector.
        constexpr V operator * (T k) const
        {
            V ret;
            for (auto i = 0; i < N; ++i) {
                ret[i] = data[i] * k;
            }
            return ret;
        }


        /// Divide vector by scalar.
        constexpr V operator / (T k) const
        {
            V ret;
            for (auto i = 0; i < N; ++i) {
                ret[i] = data[i] / k;
            }
            return ret;
        }


        /// Return the smallest coordinate of the vector
        constexpr T min_coord() const
        {
            return *std::min_element(begin(), end());
        }

        /// Return the largest coordinate of the vector
        constexpr T max_coord() const
        {
            return *std::max_element(begin(), end());
        }

    };

    template <class T, class V, unsigned int N>
    constexpr std::ostream& operator<<(std::ostream&os, const ArithVec<T,V,N>& v)
    {
        os << "[ ";
        for(const T& x : v) os << x << " ";
        os << "]";
        return os;
    }

        /// Get from operator for ArithVec descendants.
        template <class T,class V, unsigned int N>
        constexpr std::istream& operator>>(std::istream&is, ArithVec<T,V,N>& v)
        {
            is >> std::ws;
            if (is.peek() == '[')
                is.ignore();
            is >> std::ws;
            for (int c=0; c<N; ++c)
            {
                is >> v(c) >> std::ws;
            }
            if (is.peek() == ']')
                is.ignore();
            return is;
        }


        /** Dot product for two vectors. The `*' operator is
         reserved for coordinatewise	multiplication of vectors. */
        template <class T,class V, unsigned int N>
        constexpr T dot(const ArithVec<T,V,N>& v0, const ArithVec<T,V,N>& v1)
        {
            T acc = 0;
            for (auto i = 0; i < N; ++i) {
                acc += v0[i] * v1[i];
            }
            return acc;
        }

        /** Compute the sqr length by taking dot product of vector with itself. */
        template <class T,class V, unsigned int N>
        constexpr T sqr_length(const ArithVec<T,V,N>& v)
        {
            return dot(v,v);
        }

        /** Multiply double onto vector. This operator handles the case
         where the vector is on the righ side of the `*'.

         \note It seems to be optimal to put the binary operators inside the
         ArithVec class template, but the operator functions whose
         left operand is _not_ a vector cannot be inside, hence they
         are here.
         We need three operators for scalar * vector although they are
         identical, because, if we use a separate template argument for
         the left operand, it will match any type. If we use just T as
         type for the left operand hoping that other built-in types will
         be automatically converted, we will be disappointed. It seems that
         a float * ArithVec<float,Vec3f,3> function is not found if the
         left operand is really a double.
         */

        template<class T, class V, unsigned int N>
        constexpr V operator * (double k, const ArithVec<T,V,N>& v)
        {
            return v * k;
        }

        /** Multiply float onto vector. See the note in the documentation
         regarding multiplication of a double onto a vector. */
        template<class T, class V, unsigned int N>
        constexpr V operator * (float k, const ArithVec<T,V,N>& v)
        {
            return v * k;
        }

        /** Multiply unsigned int onto vector. See the note in the documentation
         regarding multiplication of a double onto a vector. */
        template<class T, class V, unsigned int N>
        constexpr V operator * (int k, const ArithVec<T,V,N>& v)
        {
            return v * k;
        }

        /** Returns the vector containing for each coordinate the smallest
         value from two vectors. */
        template <class T,class V, unsigned int N>
        constexpr V v_min(const ArithVec<T,V,N>& v0, const ArithVec<T,V,N>& v1)
        {
            V v;
            std::transform(v0.begin(), v0.end(), v1.begin(), v.begin(),
                           [](T a, T b){return (std::min)(a,b);});
            return v;
        }

        /** Returns the vector containing for each coordinate the largest
         value from two vectors. */
        template <class T,class V, unsigned int N>
        constexpr V v_max(const ArithVec<T,V,N>& v0, const ArithVec<T,V,N>& v1)
        {
            V v;
            std::transform(v0.begin(), v0.end(), v1.begin(), v.begin(),
                           [](T a, T b){return (std::max)(a,b);});
            return v;
        }


}

#endif
