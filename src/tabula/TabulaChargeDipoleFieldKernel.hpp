//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center charge-dipole field kernel — the transpose of the matrix kernel:
//  the same per-axis field integral `(a|(r−N)/|r−N|³|c)`, but contracted on the
//  AO side with a density block to give the electric field at each external
//  point `E_i(R_N) = Σ_ac (a|(r_i−R_N)/|r−R_N|³|c) D_ac`, instead of contracted
//  on the point side with dipole moments to give an AO matrix.
//

#ifndef TabulaChargeDipoleFieldKernel_hpp
#define TabulaChargeDipoleFieldKernel_hpp

#include "TabulaKernelBlockData.hpp"

namespace tabula {  // tabula namespace

/// @brief Accumulates one (bra-block, ket-block) pair's contribution to the
/// per-point electric field `E_i(R_N) = Σ_ac (a|(r_i−R_N)/|r−R_N|³|c) D_ac`.
///
/// The transpose of `charge_dipole_kernel`: phase 1 builds the per-point
/// expanded-V (kept per point, the three axes sharing Boys / Q-powers /
/// geometry), then phase 3 folds the density block into the M-rows to a
/// per-`ev` weight and accumulates `E_i[N] += Σ_ev w_i[ev]·V_i[ev,N]`.
///
/// @param l_a The bra angular momentum.
/// @param l_c The ket angular momentum.
/// @param bra The bra basis-function block data.
/// @param bra_begin The first bra contracted-GTO index of the task range.
/// @param bra_end The one-past-last bra contracted-GTO index of the range.
/// @param ket The ket basis-function block data (the full ket block — dense).
/// @param point_x The x coordinate of each external point.
/// @param point_y The y coordinate of each external point.
/// @param point_z The z coordinate of each external point.
/// @param n_points The number of external points.
/// @param density_block The gathered density sub-block, component-major:
/// `density_block[out_row·cdim + (i−bra_begin)·ket.ncgtos + jj]`, with
/// `out_row = (bra-component)·(2l_c+1) + (ket-component)` and
/// `cdim = (bra_end−bra_begin)·ket.ncgtos`.
/// @param weight The block-pair multiplicity (1 on the diagonal, 2 off it —
/// the field and the density are both symmetric).
/// @param field The field accumulator, axis-major: `field[axis·n_points + N]`.
/// It is accumulated into (the caller zeroes and reduces it).
auto charge_dipole_field_kernel(int                    l_a,
                                int                    l_c,
                                const KernelBlockData &bra,
                                int                    bra_begin,
                                int                    bra_end,
                                const KernelBlockData &ket,
                                const double          *point_x,
                                const double          *point_y,
                                const double          *point_z,
                                int                    n_points,
                                const double          *density_block,
                                double                 weight,
                                double                *field) -> void;

}  // namespace tabula

#endif /* TabulaChargeDipoleFieldKernel_hpp */
