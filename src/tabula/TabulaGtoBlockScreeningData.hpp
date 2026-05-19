//
//  Tabula — custom-recursion molecular-integral machinery.
//  Per-CGTO screening bounds of a basis-function block.
//

#ifndef TabulaGtoBlockScreeningData_hpp
#define TabulaGtoBlockScreeningData_hpp

#include <vector>

#include "GtoBlock.hpp"

namespace tabula {  // tabula namespace

/// @brief Per-contracted-GTO screening bounds derived from one CGtoBlock.
///
/// These are the operator-agnostic magnitude/geometry inputs to the integral
/// screening estimates (see docs/screening-predicates.md): the largest
/// contraction-coefficient magnitude and the exponent range of each
/// contracted GTO in the block. One entry per CGTO, in block-local order.
struct GtoBlockScreeningData
{
    /// @brief Per CGTO, the largest contraction-coefficient magnitude
    /// (`max |norm|` over its primitives).
    std::vector<double> maxCoefficient;

    /// @brief Per CGTO, the smallest primitive exponent.
    std::vector<double> minExponent;

    /// @brief Per CGTO, the largest primitive exponent.
    std::vector<double> maxExponent;
};

/// @brief Computes the per-CGTO screening bounds of a basis-function block.
///
/// VeloxChem's CGtoBlock does not guarantee sorted exponents, so the exponent
/// range is taken by scanning the block's primitive-major arrays.
/// @param block The basis-function block.
/// @return The per-CGTO screening bounds.
auto make_screening_data(const CGtoBlock& block) -> GtoBlockScreeningData;

}  // namespace tabula

#endif /* TabulaGtoBlockScreeningData_hpp */
