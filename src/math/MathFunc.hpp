#ifndef MathFunc_hpp
#define MathFunc_hpp

#include <array>
#include <cmath>
#include <vector>

#include "CustomConstrains.hpp"

namespace mathfunc {  // mathfunc

/// @brief Compares two floating point values.
/// @tparam T The floating point type.
/// @param flhs The left-hand-side floating point value to compare.
/// @param frhs The rigth-hand-side floating point value to compare.
/// @param rtol The relative tollerance.
/// @param atol The absolute tollerance.
/// @return True if floating point values are equal, False otherwise.
template <FloatingPoint T>
inline auto
equal(const T flhs, const T frhs, const T rtol, const T atol) -> bool
{
    return std::fabs(flhs - frhs) < std::max(atol, rtol * std::max(std::fabs(flhs), std::fabs(frhs)));
}

}  // namespace mathfunc

#endif /* MathFunc_hpp */
