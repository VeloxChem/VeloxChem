#ifndef NuclearPotentialGeom020SumRecGF_hpp
#define NuclearPotentialGeom020SumRecGF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "GtoBlock.hpp"
#include "NuclearPotentialGeom010PrimRecDD.hpp"
#include "NuclearPotentialGeom010PrimRecDF.hpp"
#include "NuclearPotentialGeom010PrimRecFF.hpp"
#include "NuclearPotentialGeom010PrimRecPD.hpp"
#include "NuclearPotentialGeom010PrimRecPF.hpp"
#include "NuclearPotentialGeom010PrimRecPP.hpp"
#include "NuclearPotentialGeom010PrimRecSD.hpp"
#include "NuclearPotentialGeom010PrimRecSF.hpp"
#include "NuclearPotentialGeom010PrimRecSP.hpp"
#include "NuclearPotentialGeom010PrimRecSS.hpp"
#include "NuclearPotentialGeom020PrimRecDD.hpp"
#include "NuclearPotentialGeom020PrimRecDF.hpp"
#include "NuclearPotentialGeom020PrimRecDP.hpp"
#include "NuclearPotentialGeom020PrimRecFD.hpp"
#include "NuclearPotentialGeom020PrimRecFF.hpp"
#include "NuclearPotentialGeom020PrimRecGF.hpp"
#include "NuclearPotentialGeom020PrimRecPD.hpp"
#include "NuclearPotentialGeom020PrimRecPF.hpp"
#include "NuclearPotentialGeom020PrimRecPP.hpp"
#include "NuclearPotentialGeom020PrimRecPS.hpp"
#include "NuclearPotentialGeom020PrimRecSD.hpp"
#include "NuclearPotentialGeom020PrimRecSF.hpp"
#include "NuclearPotentialGeom020PrimRecSP.hpp"
#include "NuclearPotentialGeom020PrimRecSS.hpp"
#include "NuclearPotentialPrimRecDF.hpp"
#include "NuclearPotentialPrimRecPD.hpp"
#include "NuclearPotentialPrimRecPF.hpp"
#include "NuclearPotentialPrimRecSD.hpp"
#include "NuclearPotentialPrimRecSF.hpp"
#include "NuclearPotentialPrimRecSP.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes (G|AG(2)|F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_020_gf(T&                               distributor,
                                       const CGtoBlock&                 bra_gto_block,
                                       const CGtoBlock&                 ket_gto_block,
                                       const std::pair<size_t, size_t>& bra_indices,
                                       const std::pair<size_t, size_t>& ket_indices,
                                       const bool                       bra_eq_ket) -> void
{
    // intialize external coordinate(s)

    const auto coords = distributor.coordinates();

    // intialize external quadrupoles data

    const auto quadrupoles = distributor.data();

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

    CSimdArray<double> pbuffer(8911, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(900, 1);

    CSimdArray<double> sbuffer(378, 1);

    // setup Boys function data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_data(11, ket_npgtos);

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

                    t2cfunc::comp_boys_args(bf_data, 10, factors, 17, a_exp);

                    bf_table.compute(bf_data, 0, 10);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 1, 0, bf_data, 1, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 2, 0, bf_data, 2, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 3, 0, bf_data, 3, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 4, 0, bf_data, 4, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 5, 0, bf_data, 5, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 6, 0, bf_data, 6, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 7, 0, bf_data, 7, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 8, 0, bf_data, 8, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 9, 0, bf_data, 9, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 10, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 13, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 16, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 19, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 22, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 25, 7, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 28, 8, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 31, 1, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 37, 2, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 43, 3, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 49, 4, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 55, 5, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 61, 6, 7, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 67, 7, 8, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 73, 8, 9, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 79, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 82, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 85, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 88, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 91, 6, 7, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 94, 2, 10, 13, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 103, 3, 13, 16, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 112, 4, 16, 19, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 121, 5, 19, 22, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 130, 6, 22, 25, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 139, 7, 25, 28, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 148, 10, 31, 37, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 166, 13, 37, 43, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 184, 16, 43, 49, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 202, 19, 49, 55, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 220, 22, 55, 61, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 238, 25, 61, 67, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 256, 28, 67, 73, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 274, 2, 3, 79, 82, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 280, 3, 4, 82, 85, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 286, 4, 5, 85, 88, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 292, 5, 6, 88, 91, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 298, 10, 13, 79, 94, 103, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 316, 13, 16, 82, 103, 112, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 334, 16, 19, 85, 112, 121, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 352, 19, 22, 88, 121, 130, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 370, 22, 25, 91, 130, 139, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 388, 31, 37, 94, 148, 166, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 424, 37, 43, 103, 166, 184, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 460, 43, 49, 112, 184, 202, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 496, 49, 55, 121, 202, 220, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 532, 55, 61, 130, 220, 238, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 568, 61, 67, 139, 238, 256, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 604, 79, 82, 274, 280, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 614, 82, 85, 280, 286, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 624, 85, 88, 286, 292, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 634, 94, 103, 274, 298, 316, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 664, 103, 112, 280, 316, 334, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 694, 112, 121, 286, 334, 352, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 724, 121, 130, 292, 352, 370, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 754, 148, 166, 298, 388, 424, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 814, 166, 184, 316, 424, 460, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 874, 184, 202, 334, 460, 496, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 934, 202, 220, 352, 496, 532, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 994, 220, 238, 370, 532, 568, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 1054, 10, 31, 37, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 1072, 13, 37, 43, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 1090, 16, 43, 49, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 1108, 19, 49, 55, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 1126, 10, 13, 79, 94, 103, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 1153, 13, 16, 82, 103, 112, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 1180, 16, 19, 85, 112, 121, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 1207, 31, 37, 94, 148, 166, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 1261, 37, 43, 103, 166, 184, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 1315, 43, 49, 112, 184, 202, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 1369, 49, 55, 121, 202, 220, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 1423, 79, 82, 274, 280, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 1441, 82, 85, 280, 286, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 1459, 94, 103, 274, 298, 316, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 1513, 103, 112, 280, 316, 334, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 1567, 112, 121, 286, 334, 352, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 1621, 148, 166, 298, 388, 424, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 1729, 166, 184, 316, 424, 460, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 1837, 184, 202, 334, 460, 496, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 1945, 202, 220, 352, 496, 532, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 2053, 274, 280, 604, 614, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 2083, 280, 286, 614, 624, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 2113, 298, 316, 604, 634, 664, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 2203, 316, 334, 614, 664, 694, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 2293, 334, 352, 624, 694, 724, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 2383, 388, 424, 634, 754, 814, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 2563, 424, 460, 664, 814, 874, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 2743, 460, 496, 694, 874, 934, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 2923, 496, 532, 724, 934, 994, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dp(pbuffer, 3103, 148, 166, 1054, 1072, 1126, 1207, 1261, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dp(pbuffer, 3211, 166, 184, 1072, 1090, 1153, 1261, 1315, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dp(pbuffer, 3319, 184, 202, 1090, 1108, 1180, 1315, 1369, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(pbuffer, 3427, 298, 316, 1126, 1153, 1423, 1459, 1513, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(pbuffer, 3535, 316, 334, 1153, 1180, 1441, 1513, 1567, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dd(pbuffer, 3643, 388, 424, 1207, 1261, 1459, 1621, 1729, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dd(pbuffer, 3859, 424, 460, 1261, 1315, 1513, 1729, 1837, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dd(pbuffer, 4075, 460, 496, 1315, 1369, 1567, 1837, 1945, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 4291, 604, 614, 1423, 1441, 2053, 2083, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_df(pbuffer, 4351, 634, 664, 1459, 1513, 2053, 2113, 2203, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_df(pbuffer, 4531, 664, 694, 1513, 1567, 2083, 2203, 2293, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_df(pbuffer, 4711, 754, 814, 1621, 1729, 2113, 2383, 2563, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_df(pbuffer, 5071, 814, 874, 1729, 1837, 2203, 2563, 2743, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_df(pbuffer, 5431, 874, 934, 1837, 1945, 2293, 2743, 2923, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fd(pbuffer, 5791, 1621, 1729, 3103, 3211, 3427, 3643, 3859, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fd(pbuffer, 6151, 1729, 1837, 3211, 3319, 3535, 3859, 4075, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ff(pbuffer, 6511, 2113, 2203, 3427, 3535, 4291, 4351, 4531, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ff(pbuffer, 6811, 2383, 2563, 3643, 3859, 4351, 4711, 5071, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ff(pbuffer, 7411, 2563, 2743, 3859, 4075, 4531, 5071, 5431, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_gf(pbuffer, 8011, 4711, 5071, 5791, 6151, 6511, 6811, 7411, factors, 11, 17, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 8011, quadrupoles, 6, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 3, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace npotrec

#endif /* NuclearPotentialGeom020SumRecGF_hpp */
