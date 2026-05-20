//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center charge-dipole kernel — `Σ_N d_N · (a | (r−N)/|r−N|³ | c)`, the
//  field of the bra·ket density at each site contracted with the site's dipole.
//

#ifndef TabulaChargeDipoleKernel_hpp
#define TabulaChargeDipoleKernel_hpp

#include "TabulaDipoleSet.hpp"
#include "TabulaKernelBlockData.hpp"

namespace tabula {  // tabula namespace

/// @brief The charge-dipole kernel — `Σ_N d_N · (a|(r−N)/|r−N|³|c)`.
///
/// An external-center operator like nuclear attraction: each per-axis field
/// integral is the nuclear-attraction recursion with a `∇_Q` shift along that
/// axis. The engine evaluates the three field axes (each its own table) with
/// the per-site weight `d_{N,axis}`, summed into the spherical block (the
/// `−∇_{R_N}` sign folded in). Reuses the nuclear three-phase shape, the coef
/// dictionary, and transpose symmetry. Covers `l = 0 … 4`.
///
/// @param l_a The bra angular momentum.
/// @param l_c The ket angular momentum.
/// @param bra The bra basis-function block data.
/// @param bra_begin The first bra contracted-GTO index of the task range.
/// @param bra_end The one-past-last bra contracted-GTO index of the range.
/// @param ket The ket basis-function block data.
/// @param dipoles The point dipoles to sum over.
/// @param spherical The output spherical block — `(2l_a+1)·(2l_c+1)` rows of
/// `cdim = (bra_end−bra_begin)·ket.ncgtos` values, row stride padded to a
/// multiple of 8, component-major. The contracted pair `(i, j)` is column
/// `(i−bra_begin)·ket.ncgtos + j`.
auto charge_dipole_kernel(const int              l_a,
                          const int              l_c,
                          const KernelBlockData &bra,
                          const int              bra_begin,
                          const int              bra_end,
                          const KernelBlockData &ket,
                          const DipoleSet       &dipoles,
                          double                *spherical) -> void;

}  // namespace tabula

#endif /* TabulaChargeDipoleKernel_hpp */
