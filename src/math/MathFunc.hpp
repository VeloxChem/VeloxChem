#ifndef MathFunc_hpp
#define MathFunc_hpp

#include <cmath>
#include <cstdint>
#include <vector>

namespace mathfunc {  // mathfunc namespace

/**
 Counts number of significant elements in vector mask.

 @param mask the vector mask.
 @return the number of significant elements in vector mask.
 */
auto countSignificantElements(const std::vector<int64_t>& mask) -> int64_t;

/**
 Sets elements of vector to zero.

 @param values the vector of values.
 */
auto zero(std::vector<double>& values) -> void;

/**
 Computes Chebtshev quadrature of second kind in [-1,1] interval.

 @param coordinates the vector of quadature coordinates.
 @param weights the vector of quadrature weights.
 @param nPoints the number of points in quadrature.
 */
auto quadChebyshevOfKindTwo(double* coordinates, double* weights, const int64_t nPoints) -> void;

/**
 Computes distance between two 3D vectors.

 @param aCoordX the Cartesian X coordinate of first vector.
 @param aCoordY the Cartesian Y coordinate of first vector.
 @param aCoordZ the Cartesian Z coordinate of first vector.
 @param bCoordX the Cartesian X coordinate of second vector.
 @param bCoordY the Cartesian Y coordinate of second vector.
 @param bCoordZ the Cartesian Z coordinate of second vector.
 @return the distance between vectors.
 */
inline auto
distance(const double aCoordX, const double aCoordY, const double aCoordZ, const double bCoordX, const double bCoordY, const double bCoordZ) -> double
{
    auto rx = aCoordX - bCoordX;

    auto ry = aCoordY - bCoordY;

    auto rz = aCoordZ - bCoordZ;

    return std::sqrt(rx * rx + ry * ry + rz * rz);
}

/**
 Determines batch size.

 @param nElements the size of data vector.
 @param rank the rank of thread or process.
 @param nodes the number of threads or processes.
 @return the size of data batch.
 */
auto batch_size(const int64_t nElements, const int64_t rank, const int64_t nodes) -> int64_t;

/**
 Determines offset of batch.

 @param nElements the number of elements in data vector.
 @param rank the rank of thread or process.
 @param nodes the number of threads or processes.
 @return the offset of batch.
 */
auto batch_offset(const int64_t nElements, const int64_t rank, const int64_t nodes) -> int64_t;

}  // namespace mathfunc

#endif /* MathFunc_hpp */
