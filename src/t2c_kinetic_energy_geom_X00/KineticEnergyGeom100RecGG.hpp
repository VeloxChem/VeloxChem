#ifndef KineticEnergyGeom100RecGG_hpp
#define KineticEnergyGeom100RecGG_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "OverlapPrimRecSS.hpp"
#include "KineticEnergyPrimRecSS.hpp"
#include "OverlapPrimRecSP.hpp"
#include "KineticEnergyPrimRecSP.hpp"
#include "OverlapPrimRecSD.hpp"
#include "KineticEnergyPrimRecSD.hpp"
#include "OverlapPrimRecSF.hpp"
#include "KineticEnergyPrimRecSF.hpp"
#include "OverlapPrimRecSG.hpp"
#include "KineticEnergyPrimRecSG.hpp"
#include "OverlapPrimRecPS.hpp"
#include "KineticEnergyPrimRecPS.hpp"
#include "OverlapPrimRecPP.hpp"
#include "KineticEnergyPrimRecPP.hpp"
#include "OverlapPrimRecPD.hpp"
#include "KineticEnergyPrimRecPD.hpp"
#include "OverlapPrimRecPF.hpp"
#include "KineticEnergyPrimRecPF.hpp"
#include "OverlapPrimRecPG.hpp"
#include "KineticEnergyPrimRecPG.hpp"
#include "OverlapPrimRecDP.hpp"
#include "KineticEnergyPrimRecDP.hpp"
#include "OverlapPrimRecDD.hpp"
#include "KineticEnergyPrimRecDD.hpp"
#include "OverlapPrimRecDF.hpp"
#include "KineticEnergyPrimRecDF.hpp"
#include "OverlapPrimRecDG.hpp"
#include "KineticEnergyPrimRecDG.hpp"
#include "OverlapPrimRecFD.hpp"
#include "KineticEnergyPrimRecFD.hpp"
#include "OverlapPrimRecFF.hpp"
#include "KineticEnergyPrimRecFF.hpp"
#include "OverlapPrimRecFG.hpp"
#include "KineticEnergyPrimRecFG.hpp"
#include "OverlapPrimRecGF.hpp"
#include "KineticEnergyPrimRecGF.hpp"
#include "OverlapPrimRecGG.hpp"
#include "KineticEnergyPrimRecGG.hpp"
#include "OverlapPrimRecHG.hpp"
#include "KineticEnergyPrimRecHG.hpp"
#include "GeometricalDerivatives1X0ForGY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace kinrec { // kinrec namespace

/// @brief Computes (d^(1)/dA^(1)G|T|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_kinetic_energy_geom_10_gg(T& distributor,
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

    CSimdArray<double> pbuffer(3363, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(675, 1);

    CSimdArray<double> sbuffer(243, 1);

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

                ovlrec::comp_prim_overlap_sd(pbuffer, 8, 0, 2, factors, 11, a_exp);

                kinrec::comp_prim_kinetic_energy_sd(pbuffer, 14, 0, 1, 5, 8, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_sf(pbuffer, 20, 2, 8, factors, 11, a_exp);

                kinrec::comp_prim_kinetic_energy_sf(pbuffer, 30, 2, 5, 14, 20, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_sg(pbuffer, 40, 8, 20, factors, 11, a_exp);

                kinrec::comp_prim_kinetic_energy_sg(pbuffer, 55, 8, 14, 30, 40, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_ps(pbuffer, 70, 0, factors, 8);

                kinrec::comp_prim_kinetic_energy_ps(pbuffer, 73, 1, 70, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pp(pbuffer, 76, 0, 2, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pp(pbuffer, 85, 1, 5, 76, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pd(pbuffer, 94, 2, 8, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pd(pbuffer, 112, 5, 14, 94, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pf(pbuffer, 130, 8, 20, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pf(pbuffer, 160, 14, 30, 130, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pg(pbuffer, 190, 20, 40, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pg(pbuffer, 235, 30, 55, 190, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dp(pbuffer, 280, 2, 70, 76, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_dp(pbuffer, 298, 2, 5, 73, 85, 280, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dd(pbuffer, 316, 8, 76, 94, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_dd(pbuffer, 352, 8, 14, 85, 112, 316, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_df(pbuffer, 388, 20, 94, 130, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_df(pbuffer, 448, 20, 30, 112, 160, 388, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dg(pbuffer, 508, 40, 130, 190, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_dg(pbuffer, 598, 40, 55, 160, 235, 508, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fd(pbuffer, 688, 94, 280, 316, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_fd(pbuffer, 748, 94, 112, 298, 352, 688, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_ff(pbuffer, 808, 130, 316, 388, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_ff(pbuffer, 908, 130, 160, 352, 448, 808, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fg(pbuffer, 1008, 190, 388, 508, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_fg(pbuffer, 1158, 190, 235, 448, 598, 1008, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gf(pbuffer, 1308, 388, 688, 808, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_gf(pbuffer, 1458, 388, 448, 748, 908, 1308, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gg(pbuffer, 1608, 508, 808, 1008, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_gg(pbuffer, 1833, 508, 598, 908, 1158, 1608, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_hg(pbuffer, 2058, 1008, 1308, 1608, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_hg(pbuffer, 2373, 1008, 1158, 1458, 1833, 2058, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_gx(pbuffer, 2688, 1158, 2373, 1, 15, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 2688, ket_width, ket_npgtos);
            }

            t2cfunc::transform<4, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // kinrec namespace

#endif /* KineticEnergyGeom100RecGG_hpp */
