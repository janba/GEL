#include <GEL/Geometry/normal.h>
#include <GEL/CGLA/eigensolution.h>

namespace Geometry
{
using namespace CGLA;

Vec3d smallestEigenVector(const Mat3x3d& matrix)
{
    Mat3x3d inverseMatrix = invert(matrix); // Invert the matrix
    Mat3x3d eigenvectors(0);
    Mat3x3d eigenvalues(0);

    // Call GEL's power_eigensolution on the inverse matrix
    power_eigensolution<Mat3x3d>(inverseMatrix, eigenvectors, eigenvalues, 1);

    return eigenvectors[0];
}
}
