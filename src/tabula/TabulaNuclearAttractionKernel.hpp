//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center nuclear-attraction kernel — `(a | N | c) = ∫ χ_a (1/|r−N|) χ_c`,
//  summed over a set of point charges `N`.
//

#ifndef TabulaNuclearAttractionKernel_hpp
#define TabulaNuclearAttractionKernel_hpp

#include "TabulaChargeSet.hpp"
#include "TabulaKernelBlockData.hpp"

namespace tabula {  // tabula namespace

/// @brief The nuclear-attraction kernel — `Σ_N Z_N (a|1/|r−N||c)` (VeloxChem's
/// positive convention; no physical minus).
///
/// Unlike overlap / kinetic / Coulomb this is an *external-center* integral:
/// its custom recursion runs over two vectors (`R_AC` and `Q = R_AN −
/// (γ/ζ)R_AC`), so a single straight-line kernel explodes (g|g ≈ 140k terms).
/// Instead it is **table-driven, early-contraction**: the Q-into-V transform
/// (done in codegen) makes the M matrix a pure-`R_AC` polynomial, so a generic
/// engine interprets the `(l_a, l_c)` table —
///   1. the charge sum folds into a contracted expanded-V grid `V[ev]`,
///   2. the pure-`R_AC` M·V runs per contracted pair → spherical.
/// `F_m` (the Boys function) is scalar, so the seed phase is not vectorized
/// over the ket tile. Covers `l = 0 … 4`.
///
/// @param l_a The bra angular momentum.
/// @param l_c The ket angular momentum.
/// @param bra The bra basis-function block data.
/// @param bra_begin The first bra contracted-GTO index of the task range.
/// @param bra_end The one-past-last bra contracted-GTO index of the range.
/// @param ket The ket basis-function block data.
/// @param charges The point charges to sum over.
/// @param spherical The output spherical block — `(2l_a+1)·(2l_c+1)` rows of
/// `cdim = (bra_end−bra_begin)·ket.ncgtos` values, row stride padded to a
/// multiple of 8, component-major. The contracted pair `(i, j)` is column
/// `(i−bra_begin)·ket.ncgtos + j`.
auto nuclear_attraction_kernel(const int              l_a,
                               const int              l_c,
                               const KernelBlockData &bra,
                               const int              bra_begin,
                               const int              bra_end,
                               const KernelBlockData &ket,
                               const ChargeSet       &charges,
                               double                *spherical) -> void;

}  // namespace tabula

#endif /* TabulaNuclearAttractionKernel_hpp */
