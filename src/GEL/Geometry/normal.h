#ifndef GEOMETRY_NORMAL_H
#define GEOMETRY_NORMAL_H
#pragma once

#include <GEL/CGLA/Mat.h>
#include <GEL/CGLA/eigensolution.h>
#include <GEL/CGLA/Vec.h>
#include <GEL/CGLA/ls_solve.h>

namespace Geometry
{
    #include <array>
    #include <cmath>
    #include <iostream>

    using namespace CGLA;
    inline void eigen3x3(const Mat3x3d& M, Mat3x3d& Q, std::array<double, 3>& L) {

        double m00 = M[0][0], m01 = M[0][1], m02 = M[0][2];
        double m11 = M[1][1], m12 = M[1][2];
        double m22 = M[2][2];

        double p1 = m01 * m01 + m02 * m02 + m12 * m12;

        if (p1 < 1e-12) // diagonal matrix
        {
            Q[0] = Vec3d(1, 0, 0);
            Q[1] = Vec3d(0, 1, 0);
            Q[2] = Vec3d(0, 0, 1);
            L = { m00,m11,m22 };
            return;
        }

        double q = (m00 + m11 + m22) / 3.0;
        double p2 = (m00 - q) * (m00 - q) + (m11 - q) * (m11 - q) + (m22 - q) * (m22 - q) + 2.0 * p1;
        double p = std::sqrt(p2 / 6.0);

        // normalized matrix B
        Mat3x3d B;
        B[0][0] = (m00 - q) / p; B[0][1] = m01 / p;    B[0][2] = m02 / p;
        B[1][0] = m01 / p;     B[1][1] = (m11 - q) / p; B[1][2] = m12 / p;
        B[2][0] = m02 / p;     B[2][1] = m12 / p;     B[2][2] = (m22 - q) / p;

        double r = 0.5 * (B[0][0] * B[1][1] * B[2][2] + B[0][1] * B[1][2] * B[2][0] + B[0][2] * B[1][0] * B[2][1]
            - B[0][2] * B[1][1] * B[2][0] - B[0][1] * B[1][0] * B[2][2] - B[0][0] * B[1][2] * B[2][1]);
        r = std::clamp(r, -1.0, 1.0);

        double phi = std::acos(r) / 3.0;

        double eig1 = q + 2.0 * p * std::cos(phi);
        double eig3 = q + 2.0 * p * std::cos(phi + 2.0 * M_PI / 3.0);
        double eig2 = 3.0 * q - eig1 - eig3;

        L = { eig1,eig2,eig3 };

        auto compute_eigenvector = [&](double lambda) -> Vec3d {
            Mat3x3d A = M;
            A[0][0] -= lambda; A[1][1] -= lambda; A[2][2] -= lambda;

            Vec3d r0(A[0][0], A[0][1], A[0][2]);
            Vec3d r1(A[1][0], A[1][1], A[1][2]);
            Vec3d r2(A[2][0], A[2][1], A[2][2]);

            Vec3d n0 = cross(r0, r1);
            Vec3d n1 = cross(r0, r2);
            Vec3d n2 = cross(r1, r2);

            double l0 = n0.length(), l1 = n1.length(), l2 = n2.length();
            Vec3d normal = n0;
            if (l1 > l0) normal = n1;
            if (l2 > std::max(l0, l1)) normal = n2;

            return normalize(normal);
            };

        Q[0] = compute_eigenvector(L[0]);
        Q[1] = compute_eigenvector(L[1]);
        Q[2] = compute_eigenvector(L[2]);
        return;
    }

    /// Given an input range of coordinates, estimate a normal vector
    template <std::ranges::input_range Range>
    CGLA::Vec3d estimateNormal(Range&& neighbors, double radius)
    {
        CGLA::Vec3d centroid(0.0f);
        for (const auto& point : neighbors)
            centroid += point;
        centroid /= static_cast<double>(neighbors.size());

        CGLA::Mat3x3d covariance(0.0f);
        for (const auto& point : neighbors) {
            CGLA::Vec3d diff = point - centroid;
            covariance += CGLA::outer_product(diff, diff);// * exp(-4*CGLA::sqr_length(diff)/(radius*radius));
        }

        /*CGLA::Mat3x3d eigenvectors(0), eigenvalues(0);
        CGLA::power_eigensolution<CGLA::Mat3x3d>(covariance, eigenvectors, eigenvalues, 2);
        CGLA::Vec3d norm = CGLA::normalize(CGLA::cross(eigenvectors[0], eigenvectors[1]));*/

        CGLA::Mat3x3d eigenvectors(0);
        std::array<double, 3> eigenvalues;
        eigen3x3(covariance, eigenvectors, eigenvalues);

        int min_idx = 0;
        if (eigenvalues[1] < eigenvalues[min_idx]) min_idx = 1;
        if (eigenvalues[2] < eigenvalues[min_idx]) min_idx = 2;

        CGLA::Vec3d normal = eigenvectors[min_idx];

        if (std::isnan(normal[0]) || std::isnan(normal[1]) || std::isnan(normal[2]))
            std::cout << "Nan normal!" << std::endl;
        if (normal.length() < 1e-8) {
            
            std::cout << "Zero normal!" << std::endl;
        }
        if (std::isnan(normal.length())) {
            std::cout << "Inf normal" << std::endl;
        }

        return normalize(normal);
        // CGLA::Vec3d T,B;
        // CGLA::orthogonal(norm, T, B);

        // CGLA::Mat3x3d M(0.0);
        // CGLA::Vec3d b(0.0);
        // for (auto& point : neighbors) {
        //     CGLA::Vec3d diff = point - centroid;
        //     double u = CGLA::dot(diff, T) / radius;
        //     double v = CGLA::dot(diff, B) / radius;
        //     double z = CGLA::dot(diff, norm) / radius;
        //     M += CGLA::Mat3x3d(
        //         CGLA::Vec3d(1, u, v),
        //         CGLA::Vec3d(u, u * u, u * v),
        //         CGLA::Vec3d(v, v * u, v * v));
        //     b[0] += z;
        //     b[1] += u * z;
        //     b[2] += v * z;
        // }
        // CGLA::Vec3d x = CGLA::lin_solve(M, b);
        // return CGLA::normalize(-x[1] * T - x[2] * B + norm);

    }

    /// MST edge cost for normal orientation. The first term is `1 - |n0·n1|`.
    /// The second is large when the edge leaves the tangent plane. `scale` is
    /// the inserting vertex's farthest-neighbor distance. Below 1e-12 the
    /// second term is omitted.
    inline double orientation_edge_weight(const CGLA::Vec3d& n0, const CGLA::Vec3d& n1, const CGLA::Vec3d& edge, double scale)
    {
        const double n0_length = n0.length();
        const double n1_length = n1.length();
        double angle_cost = 1.0;
        if (n0_length > 1e-12 && n1_length > 1e-12) {
            const double cosine = std::abs(CGLA::dot(n0, n1) / (n0_length * n1_length));
            angle_cost = 1.0 - std::clamp(cosine, 0.0, 1.0);
        }
        if (!(scale > 1e-12))
            return angle_cost;
        const double tunnel_cost = (std::abs(CGLA::dot(edge, n0)) + std::abs(CGLA::dot(edge, n1))) / scale;
        return angle_cost + tunnel_cost;
    }

    /// Largest neighborhood scored by estimateNormalAdaptive.
    inline constexpr int ADAPTIVE_NORMAL_K_MAX = 192;
    /// Prefixes of a distance-sorted neighborhood. Count includes the query point.
    inline constexpr std::array<int, 9> ADAPTIVE_NORMAL_SCHEDULE = {8, 12, 18, 27, 40, 60, 90, 128, 192};
    /// λ2 / λ1 below this is a line. Eigenentropy would prefer it over a plane.
    inline constexpr double ADAPTIVE_NORMAL_MIN_PLANAR_RATIO = 0.15;
    /// Scores within this band of the minimum keep the smaller neighborhood.
    inline constexpr double ADAPTIVE_NORMAL_ENTROPY_TIE = 1e-3;

    struct AdaptiveNormal {
        CGLA::Vec3d normal;
        int k = 0;
    };

    /// Normal from the planar prefix with lowest covariance eigenentropy.
    /// `neighbors` is sorted nearest-first. Below the first schedule size, the
    /// whole range is used. The sign is arbitrary.
    template <std::ranges::random_access_range Range>
        requires std::ranges::sized_range<Range>
    AdaptiveNormal estimateNormalAdaptive(Range&& neighbors)
    {
        const int n = static_cast<int>(std::ranges::size(neighbors));
        if (n <= 0)
            return {CGLA::Vec3d(0, 0, 1), 0};

        int sizes[ADAPTIVE_NORMAL_SCHEDULE.size()];
        int size_count = 0;
        if (n < ADAPTIVE_NORMAL_SCHEDULE.front()) {
            sizes[size_count++] = n;
        } else {
            for (const int k : ADAPTIVE_NORMAL_SCHEDULE) {
                if (k > n)
                    break;
                sizes[size_count++] = k;
            }
        }

        auto eigenentropy = [](double l1, double l2, double l3) {
            l1 = std::max(l1, 0.0);
            l2 = std::max(l2, 0.0);
            l3 = std::max(l3, 0.0);
            const double sum = l1 + l2 + l3;
            if (!(sum > 0.0))
                return 10.0;
            auto term = [sum](double lambda) {
                if (!(lambda > 0.0))
                    return 0.0;
                const double e = lambda / sum;
                return -e * std::log(e);
            };
            return term(l1) + term(l2) + term(l3);
        };

        CGLA::Vec3d mean(0);
        CGLA::Mat3x3d scatter(0.0);
        int seen = 0;
        int next_size = 0;

        bool have_planar = false;
        double best_entropy = 0.0;
        int best_k = sizes[0];
        CGLA::Vec3d best_normal(0, 0, 1);

        double best_planarity = -1.0;
        int fallback_k = sizes[0];
        CGLA::Vec3d fallback_normal(0, 0, 1);

        auto consider = [&](int k) {
            CGLA::Mat3x3d eigenvectors(0);
            std::array<double, 3> eigenvalues{};
            eigen3x3(scatter, eigenvectors, eigenvalues);

            int i0 = 0;
            int i1 = 1;
            int i2 = 2;
            if (eigenvalues[i0] < eigenvalues[i1])
                std::swap(i0, i1);
            if (eigenvalues[i1] < eigenvalues[i2])
                std::swap(i1, i2);
            if (eigenvalues[i0] < eigenvalues[i1])
                std::swap(i0, i1);

            const double l1 = eigenvalues[i0];
            const double l2 = eigenvalues[i1];
            const double l3 = eigenvalues[i2];
            CGLA::Vec3d normal = eigenvectors[i2];
            if (!(normal.length() > 1e-8) || !std::isfinite(normal[0]) || !std::isfinite(normal[1]) ||
                !std::isfinite(normal[2]))
                normal = CGLA::Vec3d(0, 0, 1);
            else
                normal = normalize(normal);

            const double planarity = (l1 > 1e-18) ? (l2 - l3) / l1 : 0.0;
            if (planarity > best_planarity + 1e-9) {
                best_planarity = planarity;
                fallback_k = k;
                fallback_normal = normal;
            }

            // A line has lower eigenentropy than a plane, and its least
            // eigenvector is undefined. Score entropy only on planar scales.
            // The in-plane eigenvalues are averaged first: a longer, thinner
            // neighborhood is not a better plane.
            const bool planar = l1 > 1e-18 && l2 >= ADAPTIVE_NORMAL_MIN_PLANAR_RATIO * l1;
            if (!planar)
                return;
            const double in_plane = 0.5 * (l1 + l2);
            const double entropy = eigenentropy(in_plane, in_plane, l3);
            if (!have_planar || entropy < best_entropy - ADAPTIVE_NORMAL_ENTROPY_TIE) {
                have_planar = true;
                best_entropy = entropy;
                best_k = k;
                best_normal = normal;
            }
        };

        const int last = sizes[size_count - 1];
        for (int i = 0; i < last; ++i) {
            const CGLA::Vec3d point = neighbors[i];
            const double new_count = static_cast<double>(seen + 1);
            const CGLA::Vec3d delta = point - mean;
            mean += delta / new_count;
            scatter += CGLA::outer_product(delta, delta) * (static_cast<double>(seen) / new_count);
            ++seen;
            if (seen == sizes[next_size]) {
                consider(seen);
                ++next_size;
            }
        }

        if (have_planar)
            return {best_normal, best_k};
        return {fallback_normal, fallback_k};
    }

    inline bool solve6x6(double A[6][6], double B[6], double X[6]) {
        // Augmented matrix
        float M[6][7];
        for (int r = 0; r < 6; r++) {
            for (int c = 0; c < 6; c++) M[r][c] = A[r][c];
            M[r][6] = B[r];
        }

        // Gaussian elimination with partial pivoting
        for (int k = 0; k < 6; k++)
        {
            // Pivot
            int pivot = k;
            float maxv = fabs(M[k][k]);
            for (int r = k + 1; r < 6; r++)
                if (fabs(M[r][k]) > maxv) { maxv = fabs(M[r][k]); pivot = r; }

            if (maxv < 1e-12f) return false;

            if (pivot != k)
                for (int c = 0; c < 7; c++) std::swap(M[k][c], M[pivot][c]);

            float div = M[k][k];
            for (int c = k; c < 7; c++) M[k][c] /= div;

            for (int r = 0; r < 6; r++)
            {
                if (r == k) continue;
                float f = M[r][k];
                for (int c = k; c < 7; c++)
                    M[r][c] -= f * M[k][c];
            }
        }

        for (int i = 0; i < 6; i++)
            X[i] = M[i][6];

        return true;
    }

    template <std::ranges::input_range Range>
    Vec3d estimateNormalJet(Range&& neighbors, double radius) {
        CGLA::Vec3d centroid(0.0f);
        for (const auto& point : neighbors)
            centroid += point;
        centroid /= static_cast<double>(neighbors.size());

        Vec3d nPCA = estimateNormal(neighbors, radius); 

        // Build local frame
        Vec3d tmp = fabs(nPCA[0]) < 0.9f ? Vec3d(1, 0, 0) : Vec3d(0, 1, 0);
        Vec3d u = cross(nPCA, tmp); 
        u.normalize();
        Vec3d v = cross(nPCA, u); 
        v.normalize();

        // Prepare arrays
        double A[6][6] = { 0 };
        double B[6] = { 0 };

        for (const auto& point : neighbors)
        {
            Vec3d d = point - centroid;

            double uu = dot(d, u);
            double vv = dot(d, v);
            double ww = dot(d, nPCA);

            double row[6] = { 1, uu, vv, uu * uu, uu * vv, vv * vv };

            for (int r = 0; r < 6; r++)
            {
                for (int c = 0; c < 6; c++) A[r][c] += row[r] * row[c];
                B[r] += row[r] * ww;
            }
        }

        double X[6];
        if (!solve6x6(A, B, X))
            return nPCA;

        // Get local normal
        double a1 = X[1], a2 = X[2];
        Vec3d n_local(-a1, -a2, 1.);
        if (sqr_length(n_local) < 1e-10f)
            return nPCA;
        n_local.normalize();

        Vec3d output = Vec3d(n_local[0] * u +
            n_local[1] * v +
            n_local[2] * nPCA);
        output.normalize();
        // Back to world coords
        return output;
    }

}

#endif
