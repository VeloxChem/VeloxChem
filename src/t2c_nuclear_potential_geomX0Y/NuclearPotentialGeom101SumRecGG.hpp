#ifndef NuclearPotentialGeom101SumRecGG_hpp
#define NuclearPotentialGeom101SumRecGG_hpp

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
#include "NuclearPotentialPrimRecSH.hpp"
#include "NuclearPotentialPrimRecPS.hpp"
#include "NuclearPotentialPrimRecPP.hpp"
#include "NuclearPotentialPrimRecPD.hpp"
#include "NuclearPotentialPrimRecPF.hpp"
#include "NuclearPotentialPrimRecPG.hpp"
#include "NuclearPotentialPrimRecPH.hpp"
#include "NuclearPotentialPrimRecDS.hpp"
#include "NuclearPotentialPrimRecDP.hpp"
#include "NuclearPotentialPrimRecDD.hpp"
#include "NuclearPotentialPrimRecDF.hpp"
#include "NuclearPotentialPrimRecDG.hpp"
#include "NuclearPotentialPrimRecDH.hpp"
#include "NuclearPotentialPrimRecFP.hpp"
#include "NuclearPotentialPrimRecFD.hpp"
#include "NuclearPotentialPrimRecFF.hpp"
#include "NuclearPotentialPrimRecFG.hpp"
#include "NuclearPotentialPrimRecFH.hpp"
#include "NuclearPotentialPrimRecGD.hpp"
#include "NuclearPotentialPrimRecGF.hpp"
#include "NuclearPotentialPrimRecGG.hpp"
#include "NuclearPotentialPrimRecGH.hpp"
#include "NuclearPotentialPrimRecHF.hpp"
#include "NuclearPotentialPrimRecHH.hpp"
#include "GeometricalDerivatives1X1ForGG.hpp"

#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes (d^(1)/dA^(1)G|A|d^(1)/dB^(1)G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_11_gg(T& distributor,
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

    CSimdArray<double> pbuffer(8477, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(2025, 1);

    CSimdArray<double> sbuffer(729, 1);

    // setup Boys function data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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

                    t2cfunc::comp_boys_args(bf_data, 9, factors, 17, a_exp);

                    bf_table.compute(bf_data, 0, 9);

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

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 11, 0, bf_data, 10, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 12, 1, 2, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 15, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 18, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 21, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 24, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 27, 6, 7, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 30, 7, 8, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 33, 8, 9, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 36, 9, 10, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 39, 10, 11, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 42, 1, 2, 12, 15, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 48, 2, 3, 15, 18, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 54, 3, 4, 18, 21, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 60, 4, 5, 21, 24, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 66, 5, 6, 24, 27, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 72, 6, 7, 27, 30, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 78, 7, 8, 30, 33, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 84, 8, 9, 33, 36, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 90, 9, 10, 36, 39, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 96, 12, 15, 42, 48, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 106, 15, 18, 48, 54, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 116, 18, 21, 54, 60, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 126, 21, 24, 60, 66, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 136, 24, 27, 66, 72, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 146, 27, 30, 72, 78, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 156, 30, 33, 78, 84, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 166, 33, 36, 84, 90, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 176, 42, 48, 96, 106, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 191, 48, 54, 106, 116, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 206, 54, 60, 116, 126, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 221, 60, 66, 126, 136, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 236, 66, 72, 136, 146, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 251, 72, 78, 146, 156, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 266, 78, 84, 156, 166, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sh(pbuffer, 281, 96, 106, 176, 191, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sh(pbuffer, 302, 106, 116, 191, 206, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sh(pbuffer, 323, 116, 126, 206, 221, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sh(pbuffer, 344, 126, 136, 221, 236, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sh(pbuffer, 365, 136, 146, 236, 251, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sh(pbuffer, 386, 146, 156, 251, 266, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 407, 1, 2, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 410, 2, 3, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 413, 3, 4, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 416, 4, 5, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 419, 5, 6, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 422, 1, 2, 12, 15, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 431, 2, 3, 15, 18, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 440, 3, 4, 18, 21, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 449, 4, 5, 21, 24, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 458, 5, 6, 24, 27, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 467, 12, 15, 42, 48, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 485, 15, 18, 48, 54, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 503, 18, 21, 54, 60, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 521, 21, 24, 60, 66, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 539, 24, 27, 66, 72, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 557, 42, 48, 96, 106, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 587, 48, 54, 106, 116, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 617, 54, 60, 116, 126, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 647, 60, 66, 126, 136, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 677, 66, 72, 136, 146, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 707, 96, 106, 176, 191, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 752, 106, 116, 191, 206, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 797, 116, 126, 206, 221, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 842, 126, 136, 221, 236, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 887, 136, 146, 236, 251, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ph(pbuffer, 932, 176, 191, 281, 302, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ph(pbuffer, 995, 191, 206, 302, 323, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ph(pbuffer, 1058, 206, 221, 323, 344, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ph(pbuffer, 1121, 221, 236, 344, 365, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ph(pbuffer, 1184, 236, 251, 365, 386, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 1247, 1, 2, 407, 410, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 1253, 2, 3, 410, 413, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 1259, 3, 4, 413, 416, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 1265, 4, 5, 416, 419, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 1271, 12, 15, 407, 410, 422, 431, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 1289, 15, 18, 410, 413, 431, 440, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 1307, 18, 21, 413, 416, 440, 449, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 1325, 21, 24, 416, 419, 449, 458, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 1343, 42, 48, 422, 431, 467, 485, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 1379, 48, 54, 431, 440, 485, 503, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 1415, 54, 60, 440, 449, 503, 521, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 1451, 60, 66, 449, 458, 521, 539, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 1487, 96, 106, 467, 485, 557, 587, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 1547, 106, 116, 485, 503, 587, 617, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 1607, 116, 126, 503, 521, 617, 647, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 1667, 126, 136, 521, 539, 647, 677, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 1727, 176, 191, 557, 587, 707, 752, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 1817, 191, 206, 587, 617, 752, 797, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 1907, 206, 221, 617, 647, 797, 842, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 1997, 221, 236, 647, 677, 842, 887, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dh(pbuffer, 2087, 281, 302, 707, 752, 932, 995, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dh(pbuffer, 2213, 302, 323, 752, 797, 995, 1058, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dh(pbuffer, 2339, 323, 344, 797, 842, 1058, 1121, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dh(pbuffer, 2465, 344, 365, 842, 887, 1121, 1184, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 2591, 422, 431, 1247, 1253, 1271, 1289, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 2621, 431, 440, 1253, 1259, 1289, 1307, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 2651, 440, 449, 1259, 1265, 1307, 1325, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 2681, 467, 485, 1271, 1289, 1343, 1379, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 2741, 485, 503, 1289, 1307, 1379, 1415, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 2801, 503, 521, 1307, 1325, 1415, 1451, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 2861, 557, 587, 1343, 1379, 1487, 1547, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 2961, 587, 617, 1379, 1415, 1547, 1607, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 3061, 617, 647, 1415, 1451, 1607, 1667, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fg(pbuffer, 3161, 707, 752, 1487, 1547, 1727, 1817, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fg(pbuffer, 3311, 752, 797, 1547, 1607, 1817, 1907, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fg(pbuffer, 3461, 797, 842, 1607, 1667, 1907, 1997, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fh(pbuffer, 3611, 932, 995, 1727, 1817, 2087, 2213, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fh(pbuffer, 3821, 995, 1058, 1817, 1907, 2213, 2339, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fh(pbuffer, 4031, 1058, 1121, 1907, 1997, 2339, 2465, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gd(pbuffer, 4241, 1343, 1379, 2591, 2621, 2681, 2741, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gd(pbuffer, 4331, 1379, 1415, 2621, 2651, 2741, 2801, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gf(pbuffer, 4421, 1487, 1547, 2681, 2741, 2861, 2961, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gf(pbuffer, 4571, 1547, 1607, 2741, 2801, 2961, 3061, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gg(pbuffer, 4721, 1727, 1817, 2861, 2961, 3161, 3311, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gg(pbuffer, 4946, 1817, 1907, 2961, 3061, 3311, 3461, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gh(pbuffer, 5171, 2087, 2213, 3161, 3311, 3611, 3821, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gh(pbuffer, 5486, 2213, 2339, 3311, 3461, 3821, 4031, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hf(pbuffer, 5801, 2861, 2961, 4241, 4331, 4421, 4571, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hh(pbuffer, 6011, 3611, 3821, 4721, 4946, 5171, 5486, factors, 11, 17, a_exp);

                    t2cgeom::comp_prim_op_geom_11_gg(pbuffer, 6452, 2861, 3611, 5801, 6011, 1, factors, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 6452, charges[l], ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // npotrec namespace

#endif /* NuclearPotentialGeom101SumRecGG_hpp */
