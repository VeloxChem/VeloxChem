#ifndef NuclearPotentialGeom110SumRecGD_hpp
#define NuclearPotentialGeom110SumRecGD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "GeometricalDerivatives1X0ForGY.hpp"
#include "GtoBlock.hpp"
#include "NuclearPotentialGeom010PrimRecDD.hpp"
#include "NuclearPotentialGeom010PrimRecDP.hpp"
#include "NuclearPotentialGeom010PrimRecDS.hpp"
#include "NuclearPotentialGeom010PrimRecFD.hpp"
#include "NuclearPotentialGeom010PrimRecFP.hpp"
#include "NuclearPotentialGeom010PrimRecFS.hpp"
#include "NuclearPotentialGeom010PrimRecGD.hpp"
#include "NuclearPotentialGeom010PrimRecGP.hpp"
#include "NuclearPotentialGeom010PrimRecHD.hpp"
#include "NuclearPotentialGeom010PrimRecPD.hpp"
#include "NuclearPotentialGeom010PrimRecPP.hpp"
#include "NuclearPotentialGeom010PrimRecPS.hpp"
#include "NuclearPotentialGeom010PrimRecSD.hpp"
#include "NuclearPotentialGeom010PrimRecSP.hpp"
#include "NuclearPotentialGeom010PrimRecSS.hpp"
#include "NuclearPotentialPrimRecDD.hpp"
#include "NuclearPotentialPrimRecDP.hpp"
#include "NuclearPotentialPrimRecDS.hpp"
#include "NuclearPotentialPrimRecFD.hpp"
#include "NuclearPotentialPrimRecFP.hpp"
#include "NuclearPotentialPrimRecGD.hpp"
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

/// @brief Computes (d^(1)/dA^(1)G|AG(1)|D)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_110_gd(T&                               distributor,
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

    CSimdArray<double> pbuffer(4890, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(810, 1);

    CSimdArray<double> sbuffer(405, 1);

    // setup Boys function data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 8, 0, bf_data, 8, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 9, 1, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 12, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 15, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 18, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 21, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 24, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 27, 7, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 30, 8, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 33, 1, 2, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 36, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 39, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 42, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 45, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 48, 6, 7, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 51, 1, 9, 12, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 60, 2, 12, 15, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 69, 3, 15, 18, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 78, 4, 18, 21, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 87, 5, 21, 24, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 96, 6, 24, 27, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 105, 7, 27, 30, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 114, 1, 2, 33, 36, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 120, 2, 3, 36, 39, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 126, 3, 4, 39, 42, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 132, 4, 5, 42, 45, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 138, 5, 6, 45, 48, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 144, 9, 12, 33, 51, 60, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 162, 12, 15, 36, 60, 69, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 180, 15, 18, 39, 69, 78, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 198, 18, 21, 42, 78, 87, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 216, 21, 24, 45, 87, 96, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 234, 24, 27, 48, 96, 105, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 252, 1, 2, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 255, 2, 3, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 258, 3, 4, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 261, 4, 5, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 264, 1, 9, 12, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 273, 2, 12, 15, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 282, 3, 15, 18, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 291, 4, 18, 21, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 300, 5, 21, 24, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 309, 1, 2, 33, 36, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 318, 2, 3, 36, 39, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 327, 3, 4, 39, 42, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 336, 4, 5, 42, 45, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 345, 9, 12, 33, 51, 60, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 372, 12, 15, 36, 60, 69, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 399, 15, 18, 39, 69, 78, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 426, 18, 21, 42, 78, 87, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 453, 21, 24, 45, 87, 96, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 480, 33, 36, 114, 120, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 498, 36, 39, 120, 126, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 516, 39, 42, 126, 132, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 534, 42, 45, 132, 138, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 552, 51, 60, 114, 144, 162, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 606, 60, 69, 120, 162, 180, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 660, 69, 78, 126, 180, 198, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 714, 78, 87, 132, 198, 216, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 768, 87, 96, 138, 216, 234, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 822, 1, 2, 252, 255, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 828, 2, 3, 255, 258, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 834, 3, 4, 258, 261, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 840, 9, 12, 252, 264, 273, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 858, 12, 15, 255, 273, 282, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 876, 15, 18, 258, 282, 291, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 894, 18, 21, 261, 291, 300, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 912, 33, 36, 252, 255, 309, 318, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 930, 36, 39, 255, 258, 318, 327, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 948, 39, 42, 258, 261, 327, 336, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(pbuffer, 966, 51, 60, 264, 273, 309, 345, 372, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(pbuffer, 1020, 60, 69, 273, 282, 318, 372, 399, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(pbuffer, 1074, 69, 78, 282, 291, 327, 399, 426, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(pbuffer, 1128, 78, 87, 291, 300, 336, 426, 453, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 1182, 114, 120, 309, 318, 480, 498, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 1218, 120, 126, 318, 327, 498, 516, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 1254, 126, 132, 327, 336, 516, 534, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(pbuffer, 1290, 144, 162, 345, 372, 480, 552, 606, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(pbuffer, 1398, 162, 180, 372, 399, 498, 606, 660, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(pbuffer, 1506, 180, 198, 399, 426, 516, 660, 714, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(pbuffer, 1614, 198, 216, 426, 453, 534, 714, 768, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fs(pbuffer, 1722, 264, 273, 822, 840, 858, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fs(pbuffer, 1752, 273, 282, 828, 858, 876, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fs(pbuffer, 1782, 282, 291, 834, 876, 894, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 1812, 309, 318, 822, 828, 912, 930, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 1842, 318, 327, 828, 834, 930, 948, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fp(pbuffer, 1872, 345, 372, 840, 858, 912, 966, 1020, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fp(pbuffer, 1962, 372, 399, 858, 876, 930, 1020, 1074, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fp(pbuffer, 2052, 399, 426, 876, 894, 948, 1074, 1128, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 2142, 480, 498, 912, 930, 1182, 1218, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 2202, 498, 516, 930, 948, 1218, 1254, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fd(pbuffer, 2262, 552, 606, 966, 1020, 1182, 1290, 1398, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fd(pbuffer, 2442, 606, 660, 1020, 1074, 1218, 1398, 1506, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fd(pbuffer, 2622, 660, 714, 1074, 1128, 1254, 1506, 1614, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_gp(pbuffer, 2802, 966, 1020, 1722, 1752, 1812, 1872, 1962, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_gp(pbuffer, 2937, 1020, 1074, 1752, 1782, 1842, 1962, 2052, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gd(pbuffer, 3072, 1182, 1218, 1812, 1842, 2142, 2202, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_gd(pbuffer, 3162, 1290, 1398, 1872, 1962, 2142, 2262, 2442, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_gd(pbuffer, 3432, 1398, 1506, 1962, 2052, 2202, 2442, 2622, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_hd(pbuffer, 3702, 2262, 2442, 2802, 2937, 3072, 3162, 3432, factors, 11, 17, a_exp);

                    t2cgeom::comp_prim_op_geom_10_gx(pbuffer, 4080, 2262, 3702, 3, 6, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 4080, dipoles, 3, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 2>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 2, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace npotrec

#endif /* NuclearPotentialGeom110SumRecGD_hpp */
