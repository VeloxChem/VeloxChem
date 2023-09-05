#ifndef DftFunc_hpp
#define DftFunc_hpp

#include <cstdint>
#include <vector>

#include "SubMatrix.hpp"

namespace dft {  // dft namespace

/**
 Distributes GTO values into given submatrix.

 @param submatrix the submatrix to distribute GTO values.
 @param values the vector of GTO values to distribute.
 @param irow the index of row to distribute GTO values.
*/
auto distribute(      CSubMatrix*          matrix,
                const std::vector<double>& values,
                const int64_t              irow) -> void;

/**
 Distributes GTO values into given submatrix.

 @param submatrix the submatrix to distribute GTO values.
 @param values the vector of GTO values to distribute.
 @param factor the scaling factor of values.
 @param irow the index of row to distribute GTO values.
*/
auto distribute(      CSubMatrix*          matrix,
                const std::vector<double>& values,
                const double               factor, 
                const int64_t              irow) -> void;

}  // namespace dft

#endif /* DftFunc_hpp */
