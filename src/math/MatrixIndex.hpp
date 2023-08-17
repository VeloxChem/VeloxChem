#ifndef MatrixIndex_hpp
#define MatrixIndex_hpp

#include <cstdint>

namespace mathfunc {  // mathfunc namespace

/**
 Gets upper triangular matrix index (C++ indexing scheme).

 @param i the index of row in matrix.
 @param j the index of collumn in matrix..
 @return the upper triangular matrix index.
 */
inline auto
uplo_index(const int64_t i, const int64_t j) -> int64_t
{
    return i + j * (j + 1) / 2;
}

}  // namespace mathfunc

#endif /* MatrixIndex_hpp */
