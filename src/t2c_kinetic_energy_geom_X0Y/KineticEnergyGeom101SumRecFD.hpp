#ifndef KineticEnergyGeom101SumRecFD_hpp
#define KineticEnergyGeom101SumRecFD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "GeometricalDerivatives1X1ForFD.hpp"
#include "GtoBlock.hpp"
#include "KineticEnergyPrimRecDD.hpp"
#include "KineticEnergyPrimRecDF.hpp"
#include "KineticEnergyPrimRecDP.hpp"
#include "KineticEnergyPrimRecDS.hpp"
#include "KineticEnergyPrimRecFD.hpp"
#include "KineticEnergyPrimRecFF.hpp"
#include "KineticEnergyPrimRecFP.hpp"
#include "KineticEnergyPrimRecFS.hpp"
#include "KineticEnergyPrimRecGF.hpp"
#include "KineticEnergyPrimRecGP.hpp"
#include "KineticEnergyPrimRecPD.hpp"
#include "KineticEnergyPrimRecPF.hpp"
#include "KineticEnergyPrimRecPP.hpp"
#include "KineticEnergyPrimRecPS.hpp"
#include "KineticEnergyPrimRecSD.hpp"
#include "KineticEnergyPrimRecSF.hpp"
#include "KineticEnergyPrimRecSP.hpp"
#include "KineticEnergyPrimRecSS.hpp"
#include "OverlapPrimRecDD.hpp"
#include "OverlapPrimRecDF.hpp"
#include "OverlapPrimRecDP.hpp"
#include "OverlapPrimRecDS.hpp"
#include "OverlapPrimRecFD.hpp"
#include "OverlapPrimRecFF.hpp"
#include "OverlapPrimRecFP.hpp"
#include "OverlapPrimRecFS.hpp"
#include "OverlapPrimRecGF.hpp"
#include "OverlapPrimRecGP.hpp"
#include "OverlapPrimRecPD.hpp"
#include "OverlapPrimRecPF.hpp"
#include "OverlapPrimRecPP.hpp"
#include "OverlapPrimRecPS.hpp"
#include "OverlapPrimRecSD.hpp"
#include "OverlapPrimRecSF.hpp"
#include "OverlapPrimRecSP.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace kinrec {  // kinrec namespace

/// @brief Computes (d^(1)/dA^(1)F|T|d^(1)/dB^(1)D)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_kinetic_energy_geom_11_fd(T&                               distributor,
                                   const CGtoBlock&                 bra_gto_block,
                                   const CGtoBlock&                 ket_gto_block,
                                   const std::pair<size_t, size_t>& bra_indices,
                                   const std::pair<size_t, size_t>& ket_indices,
                                   const bool                       bra_eq_ket) -> void
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

    CSimdArray<double> pbuffer(1730, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(540, 1);

    CSimdArray<double> sbuffer(315, 1);

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

                for (size_t l = 0; l < coords.size(); l++)
                {
                    ovlrec::comp_prim_overlap_sp(pbuffer, 2, 0, factors, 11);

                    kinrec::comp_prim_kinetic_energy_sp(pbuffer, 5, 1, 2, factors, 11, a_exp);

                    ovlrec::comp_prim_overlap_sd(pbuffer, 8, 0, 2, factors, 11, a_exp);

                    kinrec::comp_prim_kinetic_energy_sd(pbuffer, 14, 0, 1, 5, 8, factors, 11, a_exp);

                    ovlrec::comp_prim_overlap_sf(pbuffer, 20, 2, 8, factors, 11, a_exp);

                    kinrec::comp_prim_kinetic_energy_sf(pbuffer, 30, 2, 5, 14, 20, factors, 11, a_exp);

                    ovlrec::comp_prim_overlap_ps(pbuffer, 40, 0, factors, 8);

                    kinrec::comp_prim_kinetic_energy_ps(pbuffer, 43, 1, 40, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_pp(pbuffer, 46, 0, 2, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_pp(pbuffer, 55, 1, 5, 46, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_pd(pbuffer, 64, 2, 8, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_pd(pbuffer, 82, 5, 14, 64, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_pf(pbuffer, 100, 8, 20, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_pf(pbuffer, 130, 14, 30, 100, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_ds(pbuffer, 160, 0, 40, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_ds(pbuffer, 166, 0, 1, 43, 160, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_dp(pbuffer, 172, 2, 40, 46, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_dp(pbuffer, 190, 2, 5, 43, 55, 172, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_dd(pbuffer, 208, 8, 46, 64, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_dd(pbuffer, 244, 8, 14, 55, 82, 208, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_df(pbuffer, 280, 20, 64, 100, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_df(pbuffer, 340, 20, 30, 82, 130, 280, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_fs(pbuffer, 400, 40, 160, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_fs(pbuffer, 410, 40, 43, 166, 400, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_fp(pbuffer, 420, 46, 160, 172, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_fp(pbuffer, 450, 46, 55, 166, 190, 420, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_fd(pbuffer, 480, 64, 172, 208, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_fd(pbuffer, 540, 64, 82, 190, 244, 480, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_ff(pbuffer, 600, 100, 208, 280, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_ff(pbuffer, 700, 100, 130, 244, 340, 600, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_gp(pbuffer, 800, 172, 400, 420, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_gp(pbuffer, 845, 172, 190, 410, 450, 800, factors, 8, a_exp);

                    ovlrec::comp_prim_overlap_gf(pbuffer, 890, 280, 480, 600, factors, 8, a_exp);

                    kinrec::comp_prim_kinetic_energy_gf(pbuffer, 1040, 280, 340, 540, 700, 890, factors, 8, a_exp);

                    t2cgeom::comp_prim_op_geom_11_fd(pbuffer, 1190, 190, 340, 845, 1040, 1, factors, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 0, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<3, 2>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 2, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace kinrec

#endif /* KineticEnergyGeom101SumRecFD_hpp */
