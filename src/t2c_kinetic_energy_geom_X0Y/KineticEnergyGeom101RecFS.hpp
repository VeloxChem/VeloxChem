#ifndef KineticEnergyGeom101RecFS_hpp
#define KineticEnergyGeom101RecFS_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "OverlapPrimRecSS.hpp"
#include "KineticEnergyPrimRecSS.hpp"
#include "OverlapPrimRecSP.hpp"
#include "KineticEnergyPrimRecSP.hpp"
#include "OverlapPrimRecPS.hpp"
#include "KineticEnergyPrimRecPS.hpp"
#include "OverlapPrimRecPP.hpp"
#include "KineticEnergyPrimRecPP.hpp"
#include "OverlapPrimRecDS.hpp"
#include "KineticEnergyPrimRecDS.hpp"
#include "OverlapPrimRecDP.hpp"
#include "KineticEnergyPrimRecDP.hpp"
#include "OverlapPrimRecFS.hpp"
#include "KineticEnergyPrimRecFS.hpp"
#include "OverlapPrimRecFP.hpp"
#include "KineticEnergyPrimRecFP.hpp"
#include "OverlapPrimRecGP.hpp"
#include "KineticEnergyPrimRecGP.hpp"
#include "GeometricalDerivatives1X1ForFS.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace kinrec { // kinrec namespace

/// @brief Computes (d^(1)/dA^(1)F|T|d^(1)/dB^(1)S)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_kinetic_energy_geom_11_fs(T& distributor,
                               const CGtoBlock& bra_gto_block,
                               const CGtoBlock& ket_gto_block,
                               const std::pair<size_t, size_t>& bra_indices,
                               const std::pair<size_t, size_t>& ket_indices,
                               const bool bra_eq_ket) -> void
{
    // intialize GTOs data on bra side

    const auto bra_gto_coords = bra_gto_block.coordinates();

    const auto bra_gto_exps = bra_gto_block.exponents();

    const auto bra_gto_norms = bra_gto_block.normalization_factors();

    const auto bra_gto_indices = bra_gto_block.orbital_indices();

    const auto bra_ncgtos = bra_gto_block.number_of_basis_functions();

    const auto bra_npgtos = bra_gto_block.number_of_primitives();

    // intialize GTOs data on ket side

    const auto ket_gto_coords = ket_gto_block.coordinates();

    const auto ket_gto_exps = ket_gto_block.exponents();

    const auto ket_gto_norms = ket_gto_block.normalization_factors();

    const auto ket_gto_indices = ket_gto_block.orbital_indices();

    const auto ket_npgtos = ket_gto_block.number_of_primitives();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> factors(14, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(340, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(90, 1);

    CSimdArray<double> sbuffer(63, 1);

    // set up ket partitioning

    const auto ket_dim = ket_indices.second - ket_indices.first;

    const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());

    for (size_t i = 0; i < ket_blocks; i++)
    {
        auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);

        factors.load(ket_gto_exps, ket_range, 0, ket_npgtos);

        factors.load(ket_gto_norms, ket_range, 1, ket_npgtos);

        factors.replicate_points(ket_gto_coords, ket_range, 2, ket_npgtos);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        sbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        pbuffer.set_active_width(ket_width);

        // loop over contracted basis functions on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            cbuffer.zero();

            sbuffer.zero();

            const auto r_a = bra_gto_coords[j];

            t2cfunc::comp_distances_ab(factors, 5, 2, r_a);

            for (size_t k = 0; k < bra_npgtos; k++)
            {
                const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];

                const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];

                t2cfunc::comp_distances_pa(factors, 8, 5, a_exp);

                t2cfunc::comp_distances_pb(factors, 11, 5, a_exp);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                kinrec::comp_prim_kinetic_energy_ss(pbuffer, 1, 0, factors, a_exp);

                ovlrec::comp_prim_overlap_sp(pbuffer, 2, 0, factors, 11);

                kinrec::comp_prim_kinetic_energy_sp(pbuffer, 5, 1, 2, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_ps(pbuffer, 8, 0, factors, 8);

                kinrec::comp_prim_kinetic_energy_ps(pbuffer, 11, 1, 8, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pp(pbuffer, 14, 0, 2, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pp(pbuffer, 23, 1, 5, 14, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_ds(pbuffer, 32, 0, 8, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_ds(pbuffer, 38, 0, 1, 11, 32, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dp(pbuffer, 44, 2, 8, 14, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_dp(pbuffer, 62, 2, 5, 11, 23, 44, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fs(pbuffer, 80, 8, 32, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_fs(pbuffer, 90, 8, 11, 38, 80, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fp(pbuffer, 100, 14, 32, 44, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_fp(pbuffer, 130, 14, 23, 38, 62, 100, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gp(pbuffer, 160, 44, 80, 100, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_gp(pbuffer, 205, 44, 62, 90, 130, 160, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_11_fs(pbuffer, 250, 62, 205, 1, factors, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 250, ket_width, ket_npgtos);
            }

            t2cfunc::transform<3, 0>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 0, j, ket_range, bra_eq_ket);
        }
    }
}

} // kinrec namespace

#endif /* KineticEnergyGeom101RecFS_hpp */
