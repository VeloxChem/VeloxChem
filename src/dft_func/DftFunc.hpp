#ifndef DftFunc_hpp
#define DftFunc_hpp

#include <cstddef>
#include <vector>

#include "SubMatrix.hpp"

namespace gtoval {  // gtoval namespace

/// @brief Distributes basis function values into given submatrix.
/// @param matrix The submatrix to distribute basis function values.
/// @param values The vector of basis function values to distribute.
/// @param irow The index of row to distribute basis function values.
auto distribute(CSubMatrix* matrix, const std::vector<double>& values, const size_t irow) -> void;

// @brief Distributes basis function values into given submatrix.
/// @param matrix The submatrix to distribute basis function values.
/// @param values The vector of basis function values to distribute.
/// @param factor The scaling factor of values.
/// @param irow The index of row to distribute basis function values.
auto distribute(CSubMatrix* matrix, const std::vector<double>& values, const double factor, const size_t irow) -> void;

}  // namespace gtoval

#endif /* DftFunc_hpp */
