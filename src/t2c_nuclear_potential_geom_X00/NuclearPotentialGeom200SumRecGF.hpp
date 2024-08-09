#ifndef NuclearPotentialGeom200SumRecGF_hpp
#define NuclearPotentialGeom200SumRecGF_hpp

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
#include "NuclearPotentialPrimRecPS.hpp"
#include "NuclearPotentialPrimRecPP.hpp"
#include "NuclearPotentialPrimRecPD.hpp"
#include "NuclearPotentialPrimRecPF.hpp"
#include "NuclearPotentialPrimRecDS.hpp"
#include "NuclearPotentialPrimRecDP.hpp"
#include "NuclearPotentialPrimRecDD.hpp"
#include "NuclearPotentialPrimRecDF.hpp"
#include "NuclearPotentialPrimRecFS.hpp"
#include "NuclearPotentialPrimRecFP.hpp"
#include "NuclearPotentialPrimRecFD.hpp"
#include "NuclearPotentialPrimRecFF.hpp"
#include "NuclearPotentialPrimRecGP.hpp"
#include "NuclearPotentialPrimRecGD.hpp"
#include "NuclearPotentialPrimRecGF.hpp"
#include "NuclearPotentialPrimRecHD.hpp"
#include "NuclearPotentialPrimRecHF.hpp"
#include "NuclearPotentialPrimRecIF.hpp"
#include "GeometricalDerivatives2X0ForGY.hpp"

#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes (d^(2)/dA^(2)G|A|F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_20_gf(T& distributor,
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

    CSimdArray<double> pbuffer(4623, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(900, 1);

    CSimdArray<double> sbuffer(378, 1);

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

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 156, 1, 2, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 159, 2, 3, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 162, 3, 4, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 165, 4, 5, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 168, 5, 6, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 171, 6, 7, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 174, 1, 2, 11, 14, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 183, 2, 3, 14, 17, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 192, 3, 4, 17, 20, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 201, 4, 5, 20, 23, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 210, 5, 6, 23, 26, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 219, 6, 7, 26, 29, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 228, 11, 14, 38, 44, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 246, 14, 17, 44, 50, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 264, 17, 20, 50, 56, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 282, 20, 23, 56, 62, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 300, 23, 26, 62, 68, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 318, 26, 29, 68, 74, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 336, 38, 44, 86, 96, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 366, 44, 50, 96, 106, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 396, 50, 56, 106, 116, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 426, 56, 62, 116, 126, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 456, 62, 68, 126, 136, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 486, 68, 74, 136, 146, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 516, 1, 2, 156, 159, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 522, 2, 3, 159, 162, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 528, 3, 4, 162, 165, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 534, 4, 5, 165, 168, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 540, 5, 6, 168, 171, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 546, 11, 14, 156, 159, 174, 183, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 564, 14, 17, 159, 162, 183, 192, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 582, 17, 20, 162, 165, 192, 201, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 600, 20, 23, 165, 168, 201, 210, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 618, 23, 26, 168, 171, 210, 219, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 636, 38, 44, 174, 183, 228, 246, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 672, 44, 50, 183, 192, 246, 264, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 708, 50, 56, 192, 201, 264, 282, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 744, 56, 62, 201, 210, 282, 300, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 780, 62, 68, 210, 219, 300, 318, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 816, 86, 96, 228, 246, 336, 366, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 876, 96, 106, 246, 264, 366, 396, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 936, 106, 116, 264, 282, 396, 426, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 996, 116, 126, 282, 300, 426, 456, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 1056, 126, 136, 300, 318, 456, 486, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 1116, 156, 159, 516, 522, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 1126, 159, 162, 522, 528, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 1136, 162, 165, 528, 534, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 1146, 165, 168, 534, 540, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 1156, 174, 183, 516, 522, 546, 564, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 1186, 183, 192, 522, 528, 564, 582, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 1216, 192, 201, 528, 534, 582, 600, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 1246, 201, 210, 534, 540, 600, 618, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 1276, 228, 246, 546, 564, 636, 672, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 1336, 246, 264, 564, 582, 672, 708, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 1396, 264, 282, 582, 600, 708, 744, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 1456, 282, 300, 600, 618, 744, 780, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 1516, 336, 366, 636, 672, 816, 876, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 1616, 366, 396, 672, 708, 876, 936, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 1716, 396, 426, 708, 744, 936, 996, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 1816, 426, 456, 744, 780, 996, 1056, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gp(pbuffer, 1916, 546, 564, 1116, 1126, 1156, 1186, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gp(pbuffer, 1961, 564, 582, 1126, 1136, 1186, 1216, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gp(pbuffer, 2006, 582, 600, 1136, 1146, 1216, 1246, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gd(pbuffer, 2051, 636, 672, 1156, 1186, 1276, 1336, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gd(pbuffer, 2141, 672, 708, 1186, 1216, 1336, 1396, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gd(pbuffer, 2231, 708, 744, 1216, 1246, 1396, 1456, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gf(pbuffer, 2321, 816, 876, 1276, 1336, 1516, 1616, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gf(pbuffer, 2471, 876, 936, 1336, 1396, 1616, 1716, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gf(pbuffer, 2621, 936, 996, 1396, 1456, 1716, 1816, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hd(pbuffer, 2771, 1276, 1336, 1916, 1961, 2051, 2141, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hd(pbuffer, 2897, 1336, 1396, 1961, 2006, 2141, 2231, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hf(pbuffer, 3023, 1516, 1616, 2051, 2141, 2321, 2471, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hf(pbuffer, 3233, 1616, 1716, 2141, 2231, 2471, 2621, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_if(pbuffer, 3443, 2321, 2471, 2771, 2897, 3023, 3233, factors, 11, 17, a_exp);

                    t2cgeom::comp_prim_op_geom_20_gx(pbuffer, 3723, 816, 2321, 3443, 1, 10, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 0, charges[l], ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 3, j, ket_range, bra_eq_ket);
        }
    }
}

} // npotrec namespace

#endif /* NuclearPotentialGeom200SumRecGF_hpp */
