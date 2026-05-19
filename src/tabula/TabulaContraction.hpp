//
//  Tabula — custom-recursion molecular-integral machinery.
//  Primitive-pair contraction.
//

#ifndef TabulaContraction_hpp
#define TabulaContraction_hpp

#include <cstddef>
#include <vector>

namespace tabula {  // tabula namespace

/// @brief Contracts a primitive-pair buffer to a contracted-pair buffer —
/// step (b) of the late-contraction recursion, general to any integral.
///
/// The input is a row-major buffer of `rows` rows; each row holds
/// `cdim · nppairs` values in primitive-pair-major order (`ijoff = pp·cdim +
/// ij`), the row stride padded to a multiple of 8. For every row the routine
/// sums, for each contracted pair `ij`, its `nppairs` primitive-pair entries.
///
/// The combined primitive contraction weight is assumed already folded into
/// the input (as in `compute_overlap_seed`), so this is a plain sum. With a
/// single primitive pair per contracted pair the contraction is the identity
/// and the input buffer is returned directly.
///
/// @param primitive The primitive-pair buffer.
/// @param rows The number of rows.
/// @param cdim The number of contracted pairs.
/// @param nppairs The number of primitive pairs per contracted pair.
/// @return The contracted-pair buffer — `rows` rows of `cdim` values, the row
/// stride padded to a multiple of 8.
auto contract_primitive_pairs(const std::vector<double> &primitive,
                              const std::size_t          rows,
                              const std::size_t          cdim,
                              const std::size_t          nppairs) -> std::vector<double>;

}  // namespace tabula

#endif /* TabulaContraction_hpp */
