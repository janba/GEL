//
//  ls_solve.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 19/03/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef GEL_ls_solve_h
#define GEL_ls_solve_h

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <utility>
#include <vector>

#include <GEL/CGLA/Mat.h>
#include <GEL/CGLA/Vec.h>

namespace CGLA {

    /// Fixed-size vector for least squares past the built-in dimensions 2, 3 and 4.
    template <unsigned int N, std::floating_point Float = double>
    class VecN {
    public:
        using ScalarType = Float;

        static constexpr unsigned int get_dim() { return N; }

        constexpr VecN() noexcept = default;

        constexpr explicit VecN(Float fill) noexcept
        {
            data.fill(fill);
        }

        constexpr Float& operator[](unsigned int i) { return data[i]; }
        constexpr const Float& operator[](unsigned int i) const { return data[i]; }

        constexpr VecN& operator+=(const VecN& other)
        {
            for (unsigned int i = 0; i < N; ++i)
                data[i] += other.data[i];
            return *this;
        }

        constexpr VecN& operator-=(const VecN& other)
        {
            for (unsigned int i = 0; i < N; ++i)
                data[i] -= other.data[i];
            return *this;
        }

        constexpr VecN operator*(Float scale) const
        {
            VecN result;
            for (unsigned int i = 0; i < N; ++i)
                result.data[i] = data[i] * scale;
            return result;
        }

    private:
        std::array<Float, N> data{};
    };

    template <unsigned int N, std::floating_point Float>
    constexpr VecN<N, Float> operator*(Float scale, const VecN<N, Float>& v)
    {
        return v * scale;
    }

    /// Square matrix whose rows are VecN. Exists so lin_solve can run at any dimension.
    template <unsigned int N, std::floating_point Float = double>
    class MatN {
    public:
        using ScalarType = Float;

        static constexpr unsigned int get_v_dim() { return N; }
        static constexpr unsigned int get_h_dim() { return N; }

        constexpr MatN() noexcept = default;

        constexpr explicit MatN(Float fill) noexcept
        {
            for (auto& row : rows)
                row = VecN<N, Float>(fill);
        }

        constexpr VecN<N, Float>& operator[](unsigned int i) { return rows[i]; }
        constexpr const VecN<N, Float>& operator[](unsigned int i) const { return rows[i]; }

    private:
        std::array<VecN<N, Float>, N> rows{};
    };

    template <unsigned int N, std::floating_point Float>
    struct VecT_to_MatT<VecN<N, Float>> {
        using MatT = MatN<N, Float>;
    };

    /// Gaussian elimination with partial pivoting.
    /// V is a built-in vector (Vec2, Vec3, Vec4) or a VecN. The matrix argument
    /// is the square matrix VecT_to_MatT maps that vector to.
    template<class V>
    V lin_solve(const typename VecT_to_MatT<V>::MatT& Ap, const V& b)
    {
        using M = typename VecT_to_MatT<V>::MatT;
        using S = typename V::ScalarType;
        M A = Ap;
        V x = b;
        size_t dim = A.get_v_dim();

        // Forward elimination with partial pivoting
        for (size_t i = 0; i < dim; ++i) {
            // Pivoting
            size_t max_row = i;
            for (size_t k = i + 1; k < dim; ++k) {
                if (std::abs(A[k][i]) > std::abs(A[max_row][i]))
                    max_row = k;
            }
            std::swap(A[i], A[max_row]);
            std::swap(x[i], x[max_row]);

            // Eliminate
            S pivot = A[i][i];
            for (size_t j = i + 1; j < dim; ++j) {
                S factor = A[j][i] / pivot;
                A[j] -= factor * A[i];
                x[j] -= factor * x[i];
            }
        }

        // Back substitution
        V result(dim);
        for (int i = dim - 1; i >= 0; --i) {
            S sum = x[i];
            for (size_t j = i + 1; j < dim; ++j)
                sum -= A[i][j] * result[j];
            result[i] = sum / A[i][i];
        }

        return result;
    }

    /// Least squares for the same vector types as lin_solve. Each row of A is one
    /// sample, and b holds the corresponding scalar observations.
    template<class V>
    V ls_solve(const std::vector<V>& A, const std::vector<typename V::ScalarType>& b)
    {
        using M = typename VecT_to_MatT<V>::MatT;
        using S = typename V::ScalarType;
        V ATb(0);
        M ATA(0);
        for (size_t n = 0; n < A.size(); ++n) {
            for (unsigned int i = 0; i < ATb.get_dim(); ++i)
                ATA[i] += A[n] * A[n][i];
            ATb += b[n] * A[n];
        }
        // One ridge on the normal matrix. Inside the sample loop this grew with the
        // number of rows and biased small stencils.
        for (unsigned int i = 0; i < ATb.get_dim(); ++i)
            ATA[i][i] += S(1e-10);
        return lin_solve(ATA, ATb);
    }

}

#endif
