#ifndef NuclearPotentialGeom010SumErfRecGG_hpp
#define NuclearPotentialGeom010SumErfRecGG_hpp

#include <cstddef>
#include <array>
#include <vector>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "OverlapPrimRecSS.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "NuclearPotentialGeom010PrimRecSS.hpp"
#include "NuclearPotentialPrimRecSP.hpp"
#include "NuclearPotentialGeom010PrimRecSP.hpp"
#include "NuclearPotentialPrimRecSD.hpp"
#include "NuclearPotentialGeom010PrimRecSD.hpp"
#include "NuclearPotentialPrimRecSF.hpp"
#include "NuclearPotentialGeom010PrimRecSF.hpp"
#include "NuclearPotentialPrimRecSG.hpp"
#include "NuclearPotentialGeom010PrimRecSG.hpp"
#include "NuclearPotentialGeom010PrimRecPP.hpp"
#include "NuclearPotentialPrimRecPD.hpp"
#include "NuclearPotentialGeom010PrimRecPD.hpp"
#include "NuclearPotentialPrimRecPF.hpp"
#include "NuclearPotentialGeom010PrimRecPF.hpp"
#include "NuclearPotentialPrimRecPG.hpp"
#include "NuclearPotentialGeom010PrimRecPG.hpp"
#include "NuclearPotentialGeom010PrimRecDD.hpp"
#include "NuclearPotentialPrimRecDF.hpp"
#include "NuclearPotentialGeom010PrimRecDF.hpp"
#include "NuclearPotentialPrimRecDG.hpp"
#include "NuclearPotentialGeom010PrimRecDG.hpp"
#include "NuclearPotentialGeom010PrimRecFF.hpp"
#include "NuclearPotentialPrimRecFG.hpp"
#include "NuclearPotentialGeom010PrimRecFG.hpp"
#include "NuclearPotentialGeom010PrimRecGG.hpp"
#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes (G|Erf(AG(1))|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param omegas The vector of range-separation factors.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_erf_nuclear_potential_geom_010_gg(T& distributor,
                                           const std::vector<double>& omegas,
                                           const CGtoBlock& bra_gto_block,
                                           const CGtoBlock& ket_gto_block,
                                           const std::pair<size_t, size_t>& bra_indices,
                                           const std::pair<size_t, size_t>& ket_indices,
                                           const bool bra_eq_ket) -> void
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

    CSimdArray<double> pbuffer(6609, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(675, 1);

    CSimdArray<double> sbuffer(243, 1);

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

                t2cfunc::comp_distances_pa_from_p(factors, 11 , 8, r_a);

                t2cfunc::comp_distances_pb_from_p(factors, 14 , 8, 2);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                for (size_t l = 0; l < coords.size(); l++)
                {
                    t2cfunc::comp_distances_pc(factors, 17, 8, coords[l]);

                    t2cfunc::comp_boys_args(bf_data, 10, factors, 17, a_exp, omegas[l]);

                    bf_table.compute(bf_data, 0, 10, factors, a_exp, omegas[l]);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 1, 0, bf_data, 1, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 2, 0, bf_data, 2, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 3, 0, bf_data, 3, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 4, 0, bf_data, 4, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 5, 0, bf_data, 5, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 6, 0, bf_data, 6, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 7, 0, bf_data, 7, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 8, 0, bf_data, 8, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 9, 0, bf_data, 9, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 10, 1, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 13, 2, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 16, 3, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 19, 4, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 22, 5, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 25, 6, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 28, 7, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 31, 8, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 34, 9, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 37, 1, 2, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 40, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 43, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 46, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 49, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 52, 6, 7, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 55, 7, 8, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 58, 1, 10, 13, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 67, 2, 13, 16, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 76, 3, 16, 19, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 85, 4, 19, 22, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 94, 5, 22, 25, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 103, 6, 25, 28, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 112, 7, 28, 31, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 121, 8, 31, 34, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 130, 1, 2, 37, 40, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 136, 2, 3, 40, 43, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 142, 3, 4, 43, 46, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 148, 4, 5, 46, 49, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 154, 5, 6, 49, 52, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 160, 6, 7, 52, 55, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 166, 10, 13, 37, 58, 67, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 184, 13, 16, 40, 67, 76, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 202, 16, 19, 43, 76, 85, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 220, 19, 22, 46, 85, 94, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 238, 22, 25, 49, 94, 103, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 256, 25, 28, 52, 103, 112, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 274, 28, 31, 55, 112, 121, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 292, 37, 40, 130, 136, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 302, 40, 43, 136, 142, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 312, 43, 46, 142, 148, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 322, 46, 49, 148, 154, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 332, 49, 52, 154, 160, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 342, 58, 67, 130, 166, 184, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 372, 67, 76, 136, 184, 202, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 402, 76, 85, 142, 202, 220, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 432, 85, 94, 148, 220, 238, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 462, 94, 103, 154, 238, 256, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 492, 103, 112, 160, 256, 274, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 522, 130, 136, 292, 302, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 537, 136, 142, 302, 312, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 552, 142, 148, 312, 322, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 567, 148, 154, 322, 332, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 582, 166, 184, 292, 342, 372, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 627, 184, 202, 302, 372, 402, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 672, 202, 220, 312, 402, 432, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 717, 220, 238, 322, 432, 462, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 762, 238, 256, 332, 462, 492, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 807, 10, 13, 37, 58, 67, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 834, 13, 16, 40, 67, 76, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 861, 16, 19, 43, 76, 85, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 888, 19, 22, 46, 85, 94, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 915, 37, 40, 130, 136, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 933, 40, 43, 136, 142, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 951, 43, 46, 142, 148, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 969, 58, 67, 130, 166, 184, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 1023, 67, 76, 136, 184, 202, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 1077, 76, 85, 142, 202, 220, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 1131, 85, 94, 148, 220, 238, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 1185, 130, 136, 292, 302, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 1215, 136, 142, 302, 312, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 1245, 142, 148, 312, 322, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 1275, 166, 184, 292, 342, 372, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 1365, 184, 202, 302, 372, 402, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 1455, 202, 220, 312, 402, 432, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 1545, 220, 238, 322, 432, 462, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 1635, 292, 302, 522, 537, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 1680, 302, 312, 537, 552, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 1725, 312, 322, 552, 567, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 1770, 342, 372, 522, 582, 627, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 1905, 372, 402, 537, 627, 672, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 2040, 402, 432, 552, 672, 717, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 2175, 432, 462, 567, 717, 762, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(pbuffer, 2310, 166, 184, 807, 834, 915, 969, 1023, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(pbuffer, 2418, 184, 202, 834, 861, 933, 1023, 1077, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dd(pbuffer, 2526, 202, 220, 861, 888, 951, 1077, 1131, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 2634, 292, 302, 915, 933, 1185, 1215, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 2694, 302, 312, 933, 951, 1215, 1245, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_df(pbuffer, 2754, 342, 372, 969, 1023, 1185, 1275, 1365, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_df(pbuffer, 2934, 372, 402, 1023, 1077, 1215, 1365, 1455, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_df(pbuffer, 3114, 402, 432, 1077, 1131, 1245, 1455, 1545, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 3294, 522, 537, 1185, 1215, 1635, 1680, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 3384, 537, 552, 1215, 1245, 1680, 1725, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dg(pbuffer, 3474, 582, 627, 1275, 1365, 1635, 1770, 1905, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dg(pbuffer, 3744, 627, 672, 1365, 1455, 1680, 1905, 2040, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dg(pbuffer, 4014, 672, 717, 1455, 1545, 1725, 2040, 2175, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ff(pbuffer, 4284, 1275, 1365, 2310, 2418, 2634, 2754, 2934, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ff(pbuffer, 4584, 1365, 1455, 2418, 2526, 2694, 2934, 3114, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fg(pbuffer, 4884, 1635, 1680, 2634, 2694, 3294, 3384, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fg(pbuffer, 5034, 1770, 1905, 2754, 2934, 3294, 3474, 3744, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fg(pbuffer, 5484, 1905, 2040, 2934, 3114, 3384, 3744, 4014, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_gg(pbuffer, 5934, 3474, 3744, 4284, 4584, 4884, 5034, 5484, factors, 11, 17, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 5934, dipoles, 3, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // npotrec namespace

#endif /* NuclearPotentialGeom010SumErfRecGG_hpp */
