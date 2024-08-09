#ifndef NuclearPotentialGeom101SumRecGF_hpp
#define NuclearPotentialGeom101SumRecGF_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "OverlapPrimRecSS.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "NuclearPotentialPrimRecSP.hpp"
#include "NuclearPotentialPrimRecSD.hpp"
#include "NuclearPotentialPrimRecSF.hpp"
#include "NuclearPotentialPrimRecSG.hpp"
#include "NuclearPotentialPrimRecPS.hpp"
#include "NuclearPotentialPrimRecPP.hpp"
#include "NuclearPotentialPrimRecPD.hpp"
#include "NuclearPotentialPrimRecPF.hpp"
#include "NuclearPotentialPrimRecPG.hpp"
#include "NuclearPotentialPrimRecDS.hpp"
#include "NuclearPotentialPrimRecDP.hpp"
#include "NuclearPotentialPrimRecDD.hpp"
#include "NuclearPotentialPrimRecDF.hpp"
#include "NuclearPotentialPrimRecDG.hpp"
#include "NuclearPotentialPrimRecFS.hpp"
#include "NuclearPotentialPrimRecFP.hpp"
#include "NuclearPotentialPrimRecFD.hpp"
#include "NuclearPotentialPrimRecFF.hpp"
#include "NuclearPotentialPrimRecFG.hpp"
#include "NuclearPotentialPrimRecGP.hpp"
#include "NuclearPotentialPrimRecGD.hpp"
#include "NuclearPotentialPrimRecGF.hpp"
#include "NuclearPotentialPrimRecGG.hpp"
#include "NuclearPotentialPrimRecHD.hpp"
#include "NuclearPotentialPrimRecHG.hpp"
#include "GeometricalDerivatives1X1ForGF.hpp"

#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes (d^(1)/dA^(1)G|A|d^(1)/dB^(1)F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_11_gf(T& distributor,
                                      const CGtoBlock& bra_gto_block,
                                      const CGtoBlock& ket_gto_block,
                                      const std::pair<size_t, size_t>& bra_indices,
                                      const std::pair<size_t, size_t>& ket_indices,
                                      const bool bra_eq_ket) -> void
{
    // intialize external coordinate(s)

    const auto coords = distributor.coordinates();

    // intialize external charge(s)

    const auto charges = distributor.data();

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

    CSimdArray<double> pbuffer(5472, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(1350, 1);

    CSimdArray<double> sbuffer(567, 1);

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

                t2cfunc::comp_distances_pa_from_p(factors, 11 , 8, r_a);

                t2cfunc::comp_distances_pb_from_p(factors, 14 , 8, 2);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                for (size_t l = 0; l < coords.size(); l++)
                {
                    t2cfunc::comp_distances_pc(factors, 17, 8, coords[l]);

                    t2cfunc::comp_boys_args(bf_data, 8, factors, 17, a_exp);

                    bf_table.compute(bf_data, 0, 8);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 1, 0, bf_data, 0, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 2, 0, bf_data, 1, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 3, 0, bf_data, 2, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 4, 0, bf_data, 3, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 5, 0, bf_data, 4, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 6, 0, bf_data, 5, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 7, 0, bf_data, 6, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 8, 0, bf_data, 7, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 9, 0, bf_data, 8, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 10, 0, bf_data, 9, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 11, 1, 2, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 14, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 17, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 20, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 23, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 26, 6, 7, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 29, 7, 8, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 32, 8, 9, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 35, 9, 10, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 38, 1, 2, 11, 14, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 44, 2, 3, 14, 17, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 50, 3, 4, 17, 20, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 56, 4, 5, 20, 23, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 62, 5, 6, 23, 26, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 68, 6, 7, 26, 29, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 74, 7, 8, 29, 32, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 80, 8, 9, 32, 35, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 86, 11, 14, 38, 44, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 96, 14, 17, 44, 50, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 106, 17, 20, 50, 56, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 116, 20, 23, 56, 62, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 126, 23, 26, 62, 68, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 136, 26, 29, 68, 74, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 146, 29, 32, 74, 80, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 156, 38, 44, 86, 96, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 171, 44, 50, 96, 106, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 186, 50, 56, 106, 116, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 201, 56, 62, 116, 126, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 216, 62, 68, 126, 136, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 231, 68, 74, 136, 146, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 246, 1, 2, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 249, 2, 3, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 252, 3, 4, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 255, 4, 5, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 258, 5, 6, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 261, 1, 2, 11, 14, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 270, 2, 3, 14, 17, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 279, 3, 4, 17, 20, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 288, 4, 5, 20, 23, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 297, 5, 6, 23, 26, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 306, 11, 14, 38, 44, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 324, 14, 17, 44, 50, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 342, 17, 20, 50, 56, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 360, 20, 23, 56, 62, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 378, 23, 26, 62, 68, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 396, 38, 44, 86, 96, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 426, 44, 50, 96, 106, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 456, 50, 56, 106, 116, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 486, 56, 62, 116, 126, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 516, 62, 68, 126, 136, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 546, 86, 96, 156, 171, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 591, 96, 106, 171, 186, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 636, 106, 116, 186, 201, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 681, 116, 126, 201, 216, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 726, 126, 136, 216, 231, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 771, 1, 2, 246, 249, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 777, 2, 3, 249, 252, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 783, 3, 4, 252, 255, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 789, 4, 5, 255, 258, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 795, 11, 14, 246, 249, 261, 270, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 813, 14, 17, 249, 252, 270, 279, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 831, 17, 20, 252, 255, 279, 288, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 849, 20, 23, 255, 258, 288, 297, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 867, 38, 44, 261, 270, 306, 324, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 903, 44, 50, 270, 279, 324, 342, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 939, 50, 56, 279, 288, 342, 360, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 975, 56, 62, 288, 297, 360, 378, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 1011, 86, 96, 306, 324, 396, 426, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 1071, 96, 106, 324, 342, 426, 456, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 1131, 106, 116, 342, 360, 456, 486, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 1191, 116, 126, 360, 378, 486, 516, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 1251, 156, 171, 396, 426, 546, 591, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 1341, 171, 186, 426, 456, 591, 636, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 1431, 186, 201, 456, 486, 636, 681, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 1521, 201, 216, 486, 516, 681, 726, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 1611, 246, 249, 771, 777, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 1621, 249, 252, 777, 783, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 1631, 252, 255, 783, 789, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 1641, 261, 270, 771, 777, 795, 813, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 1671, 270, 279, 777, 783, 813, 831, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 1701, 279, 288, 783, 789, 831, 849, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 1731, 306, 324, 795, 813, 867, 903, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 1791, 324, 342, 813, 831, 903, 939, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 1851, 342, 360, 831, 849, 939, 975, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 1911, 396, 426, 867, 903, 1011, 1071, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 2011, 426, 456, 903, 939, 1071, 1131, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 2111, 456, 486, 939, 975, 1131, 1191, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fg(pbuffer, 2211, 546, 591, 1011, 1071, 1251, 1341, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fg(pbuffer, 2361, 591, 636, 1071, 1131, 1341, 1431, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fg(pbuffer, 2511, 636, 681, 1131, 1191, 1431, 1521, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gp(pbuffer, 2661, 795, 813, 1611, 1621, 1641, 1671, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gp(pbuffer, 2706, 813, 831, 1621, 1631, 1671, 1701, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gd(pbuffer, 2751, 867, 903, 1641, 1671, 1731, 1791, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gd(pbuffer, 2841, 903, 939, 1671, 1701, 1791, 1851, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gf(pbuffer, 2931, 1011, 1071, 1731, 1791, 1911, 2011, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gf(pbuffer, 3081, 1071, 1131, 1791, 1851, 2011, 2111, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gg(pbuffer, 3231, 1251, 1341, 1911, 2011, 2211, 2361, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gg(pbuffer, 3456, 1341, 1431, 2011, 2111, 2361, 2511, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hd(pbuffer, 3681, 1731, 1791, 2661, 2706, 2751, 2841, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hg(pbuffer, 3807, 2211, 2361, 2931, 3081, 3231, 3456, factors, 11, 17, a_exp);

                    t2cgeom::comp_prim_op_geom_11_gf(pbuffer, 4122, 1731, 2211, 3681, 3807, 1, factors, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 0, charges[l], ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 3, j, ket_range, bra_eq_ket);
        }
    }
}

} // npotrec namespace

#endif /* NuclearPotentialGeom101SumRecGF_hpp */
