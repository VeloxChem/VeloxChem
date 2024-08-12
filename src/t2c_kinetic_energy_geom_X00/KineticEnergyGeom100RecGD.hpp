#ifndef KineticEnergyGeom100RecGD_hpp
#define KineticEnergyGeom100RecGD_hpp

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
#include "OverlapPrimRecPS.hpp"
#include "KineticEnergyPrimRecPS.hpp"
#include "OverlapPrimRecPP.hpp"
#include "KineticEnergyPrimRecPP.hpp"
#include "OverlapPrimRecPD.hpp"
#include "KineticEnergyPrimRecPD.hpp"
#include "OverlapPrimRecDS.hpp"
#include "KineticEnergyPrimRecDS.hpp"
#include "OverlapPrimRecDP.hpp"
#include "KineticEnergyPrimRecDP.hpp"
#include "OverlapPrimRecDD.hpp"
#include "KineticEnergyPrimRecDD.hpp"
#include "OverlapPrimRecFS.hpp"
#include "KineticEnergyPrimRecFS.hpp"
#include "OverlapPrimRecFP.hpp"
#include "KineticEnergyPrimRecFP.hpp"
#include "OverlapPrimRecFD.hpp"
#include "KineticEnergyPrimRecFD.hpp"
#include "OverlapPrimRecGP.hpp"
#include "KineticEnergyPrimRecGP.hpp"
#include "OverlapPrimRecGD.hpp"
#include "KineticEnergyPrimRecGD.hpp"
#include "OverlapPrimRecHD.hpp"
#include "KineticEnergyPrimRecHD.hpp"
#include "GeometricalDerivatives1X0ForGY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace kinrec { // kinrec namespace

/// @brief Computes (d^(1)/dA^(1)G|T|D)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_kinetic_energy_geom_10_gd(T& distributor,
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

    CSimdArray<double> pbuffer(1192, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(270, 1);

    CSimdArray<double> sbuffer(135, 1);

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

                ovlrec::comp_prim_overlap_ps(pbuffer, 20, 0, factors, 8);

                kinrec::comp_prim_kinetic_energy_ps(pbuffer, 23, 1, 20, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pp(pbuffer, 26, 0, 2, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pp(pbuffer, 35, 1, 5, 26, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pd(pbuffer, 44, 2, 8, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_pd(pbuffer, 62, 5, 14, 44, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_ds(pbuffer, 80, 0, 20, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_ds(pbuffer, 86, 0, 1, 23, 80, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dp(pbuffer, 92, 2, 20, 26, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_dp(pbuffer, 110, 2, 5, 23, 35, 92, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dd(pbuffer, 128, 8, 26, 44, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_dd(pbuffer, 164, 8, 14, 35, 62, 128, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fs(pbuffer, 200, 20, 80, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_fs(pbuffer, 210, 20, 23, 86, 200, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fp(pbuffer, 220, 26, 80, 92, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_fp(pbuffer, 250, 26, 35, 86, 110, 220, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fd(pbuffer, 280, 44, 92, 128, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_fd(pbuffer, 340, 44, 62, 110, 164, 280, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gp(pbuffer, 400, 92, 200, 220, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_gp(pbuffer, 445, 92, 110, 210, 250, 400, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gd(pbuffer, 490, 128, 220, 280, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_gd(pbuffer, 580, 128, 164, 250, 340, 490, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_hd(pbuffer, 670, 280, 400, 490, factors, 8, a_exp);

                kinrec::comp_prim_kinetic_energy_hd(pbuffer, 796, 280, 340, 445, 580, 670, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_gx(pbuffer, 922, 340, 796, 1, 6, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 922, ket_width, ket_npgtos);
            }

            t2cfunc::transform<4, 2>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 2, j, ket_range, bra_eq_ket);
        }
    }
}

} // kinrec namespace

#endif /* KineticEnergyGeom100RecGD_hpp */
