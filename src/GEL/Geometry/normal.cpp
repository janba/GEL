#include <GEL/Geometry/normal.h>

namespace Geometry {

    Vec3d smallestEigenVector(const Mat3x3d& matrix) {
        Mat3x3d inverseMatrix = invert(matrix);  // Invert the matrix
        Mat3x3d eigenvectors(0);
        Mat3x3d eigenvalues(0);

        // Call GEL's power_eigensolution on the inverse matrix
        power_eigensolution<Mat3x3d>(inverseMatrix, eigenvectors, eigenvalues, 1);

        return eigenvectors[0];
    }

	Vec3d estimateNormal(const std::vector<Vec3d>& neighbors) {
        Vec3d centroid(0.0f);
        for (const auto& point : neighbors) {
            centroid += point;
        }
        centroid /= static_cast<double>(neighbors.size());

        Mat3x3d covariance(0.0f);
        for (const auto& point : neighbors) {
            Vec3d diff = point - centroid;
            covariance += outer_product(diff, diff);
        }
        covariance /= static_cast<double>(neighbors.size());

        Vec3d normal(0.0f);

        normal = smallestEigenVector(covariance);

        if (normal.length() < 1e-8)
            std::cout << "Zero normal!" << std::endl;

        return normalize(normal); // Normalize to ensure unit length
	}
}