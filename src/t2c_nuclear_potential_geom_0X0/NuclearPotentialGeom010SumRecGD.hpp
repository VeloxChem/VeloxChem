#ifndef NuclearPotentialGeom010SumRecGD_hpp
#define NuclearPotentialGeom010SumRecGD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "GtoBlock.hpp"
#include "NuclearPotentialGeom010PrimRecDD.hpp"
#include "NuclearPotentialGeom010PrimRecDP.hpp"
#include "NuclearPotentialGeom010PrimRecDS.hpp"
#include "NuclearPotentialGeom010PrimRecFD.hpp"
#include "NuclearPotentialGeom010PrimRecFP.hpp"
#include "NuclearPotentialGeom010PrimRecGD.hpp"
#include "NuclearPotentialGeom010PrimRecPD.hpp"
#include "NuclearPotentialGeom010PrimRecPP.hpp"
#include "NuclearPotentialGeom010PrimRecPS.hpp"
#include "NuclearPotentialGeom010PrimRecSD.hpp"
#include "NuclearPotentialGeom010PrimRecSP.hpp"
#include "NuclearPotentialGeom010PrimRecSS.hpp"
#include "NuclearPotentialPrimRecDD.hpp"
#include "NuclearPotentialPrimRecDP.hpp"
#include "NuclearPotentialPrimRecFD.hpp"
#include "NuclearPotentialPrimRecPD.hpp"
#include "NuclearPotentialPrimRecPP.hpp"
#include "NuclearPotentialPrimRecPS.hpp"
#include "NuclearPotentialPrimRecSD.hpp"
#include "NuclearPotentialPrimRecSP.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes (G|AG(1)|D)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_010_gd(T&                               distributor,
                                       const CGtoBlock&                 bra_gto_block,
                                       const CGtoBlock&                 ket_gto_block,
                                       const std::pair<size_t, size_t>& bra_indices,
                                       const std::pair<size_t, size_t>& ket_indices,
                                       const bool                       bra_eq_ket) -> void
{
    // intialize external coordinate(s)

    const auto coords = distributor.coordinates();

    // intialize external dipoles data

    const auto dipoles = distributor.data();

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

    CSimdArray<double> factors(20, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(2180, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(270, 1);

    CSimdArray<double> sbuffer(135, 1);

    // setup Boys function data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

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

        bf_data.set_active_width(ket_width);

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

                t2cfunc::comp_coordinates_p(factors, 8, 2, r_a, a_exp);

                t2cfunc::comp_distances_pa_from_p(factors, 11, 8, r_a);

                t2cfunc::comp_distances_pb_from_p(factors, 14, 8, 2);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                for (size_t l = 0; l < coords.size(); l++)
                {
                    t2cfunc::comp_distances_pc(factors, 17, 8, coords[l]);

                    t2cfunc::comp_boys_args(bf_data, 8, factors, 17, a_exp);

                    bf_table.compute(bf_data, 0, 8);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 1, 0, bf_data, 1, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 2, 0, bf_data, 2, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 3, 0, bf_data, 3, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 4, 0, bf_data, 4, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 5, 0, bf_data, 5, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 6, 0, bf_data, 6, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 7, 0, bf_data, 7, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 8, 1, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 11, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 14, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 17, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 20, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 23, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 26, 7, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 29, 1, 2, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 32, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 35, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 38, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 41, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 44, 1, 8, 11, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 53, 2, 11, 14, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 62, 3, 14, 17, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 71, 4, 17, 20, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 80, 5, 20, 23, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 89, 6, 23, 26, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 98, 1, 2, 29, 32, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 104, 2, 3, 32, 35, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 110, 3, 4, 35, 38, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 116, 4, 5, 38, 41, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(
                        pbuffer, 122, 8, 11, 29, 44, 53, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(
                        pbuffer, 140, 11, 14, 32, 53, 62, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(
                        pbuffer, 158, 14, 17, 35, 62, 71, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(
                        pbuffer, 176, 17, 20, 38, 71, 80, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(
                        pbuffer, 194, 20, 23, 41, 80, 89, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 212, 1, 2, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 215, 2, 3, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 218, 3, 4, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 221, 1, 8, 11, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 230, 2, 11, 14, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 239, 3, 14, 17, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 248, 4, 17, 20, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 257, 1, 2, 29, 32, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 266, 2, 3, 32, 35, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 275, 3, 4, 35, 38, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(
                        pbuffer, 284, 8, 11, 29, 44, 53, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(
                        pbuffer, 311, 11, 14, 32, 53, 62, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(
                        pbuffer, 338, 14, 17, 35, 62, 71, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(
                        pbuffer, 365, 17, 20, 38, 71, 80, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 392, 29, 32, 98, 104, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 410, 32, 35, 104, 110, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 428, 35, 38, 110, 116, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(
                        pbuffer, 446, 44, 53, 98, 122, 140, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(
                        pbuffer, 500, 53, 62, 104, 140, 158, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(
                        pbuffer, 554, 62, 71, 110, 158, 176, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(
                        pbuffer, 608, 71, 80, 116, 176, 194, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(
                        pbuffer, 662, 8, 11, 212, 221, 230, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(
                        pbuffer, 680, 11, 14, 215, 230, 239, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(
                        pbuffer, 698, 14, 17, 218, 239, 248, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(
                        pbuffer, 716, 29, 32, 212, 215, 257, 266, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(
                        pbuffer, 734, 32, 35, 215, 218, 266, 275, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(
                        pbuffer, 752, 44, 53, 221, 230, 257, 284, 311, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(
                        pbuffer, 806, 53, 62, 230, 239, 266, 311, 338, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(
                        pbuffer, 860, 62, 71, 239, 248, 275, 338, 365, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(
                        pbuffer, 914, 98, 104, 257, 266, 392, 410, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(
                        pbuffer, 950, 104, 110, 266, 275, 410, 428, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(
                        pbuffer, 986, 122, 140, 284, 311, 392, 446, 500, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(
                        pbuffer, 1094, 140, 158, 311, 338, 410, 500, 554, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(
                        pbuffer, 1202, 158, 176, 338, 365, 428, 554, 608, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fp(
                        pbuffer, 1310, 284, 311, 662, 680, 716, 752, 806, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fp(
                        pbuffer, 1400, 311, 338, 680, 698, 734, 806, 860, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(
                        pbuffer, 1490, 392, 410, 716, 734, 914, 950, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fd(
                        pbuffer, 1550, 446, 500, 752, 806, 914, 986, 1094, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fd(
                        pbuffer, 1730, 500, 554, 806, 860, 950, 1094, 1202, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_gd(
                        pbuffer, 1910, 986, 1094, 1310, 1400, 1490, 1550, 1730, factors, 11, 17, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 1910, dipoles, 3, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 2>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 2, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace npotrec

#endif /* NuclearPotentialGeom010SumRecGD_hpp */
