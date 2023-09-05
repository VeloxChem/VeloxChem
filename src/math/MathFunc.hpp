#ifndef MathFunc_hpp
#define MathFunc_hpp

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
auto quadChebyshevOfKindTwo(double* coordinates, double* weights, const int32_t nPoints) -> void;

}  // namespace mathfunc

#endif /* MathFunc_hpp */
