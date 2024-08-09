#ifndef NuclearPotentialGeom020SumRecFG_hpp
#define NuclearPotentialGeom020SumRecFG_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "OverlapPrimRecSS.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "NuclearPotentialGeom010PrimRecSS.hpp"
#include "NuclearPotentialGeom020PrimRecSS.hpp"
#include "NuclearPotentialPrimRecSP.hpp"
#include "NuclearPotentialGeom010PrimRecSP.hpp"
#include "NuclearPotentialGeom020PrimRecSP.hpp"
#include "NuclearPotentialPrimRecSD.hpp"
#include "NuclearPotentialGeom010PrimRecSD.hpp"
#include "NuclearPotentialGeom020PrimRecSD.hpp"
#include "NuclearPotentialPrimRecSF.hpp"
#include "NuclearPotentialGeom010PrimRecSF.hpp"
#include "NuclearPotentialGeom020PrimRecSF.hpp"
#include "NuclearPotentialPrimRecSG.hpp"
#include "NuclearPotentialGeom010PrimRecSG.hpp"
#include "NuclearPotentialGeom020PrimRecSG.hpp"
#include "NuclearPotentialGeom020PrimRecPD.hpp"
#include "NuclearPotentialGeom010PrimRecPF.hpp"
#include "NuclearPotentialGeom020PrimRecPF.hpp"
#include "NuclearPotentialPrimRecPG.hpp"
#include "NuclearPotentialGeom010PrimRecPG.hpp"
#include "NuclearPotentialGeom020PrimRecPG.hpp"
#include "NuclearPotentialGeom020PrimRecDF.hpp"
#include "NuclearPotentialGeom010PrimRecDG.hpp"
#include "NuclearPotentialGeom020PrimRecDG.hpp"
#include "NuclearPotentialGeom020PrimRecFG.hpp"
#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes (F|AG(2)|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_020_fg(T& distributor,
                                       const CGtoBlock& bra_gto_block,
                                       const CGtoBlock& ket_gto_block,
                                       const std::pair<size_t, size_t>& bra_indices,
                                       const std::pair<size_t, size_t>& ket_indices,
                                       const bool bra_eq_ket) -> void
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

    CSimdArray<double> pbuffer(6718, ket_npgtos);

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

                t2cfunc::comp_distances_pa_from_p(factors, 11 , 8, r_a);

                t2cfunc::comp_distances_pb_from_p(factors, 14 , 8, 2);

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

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 1054, 274, 280, 604, 614, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 1069, 280, 286, 614, 624, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 1084, 298, 316, 604, 634, 664, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 1129, 316, 334, 614, 664, 694, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 1174, 334, 352, 624, 694, 724, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 1219, 388, 424, 634, 754, 814, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 1309, 424, 460, 664, 814, 874, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 1399, 460, 496, 694, 874, 934, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 1489, 496, 532, 724, 934, 994, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 1579, 148, 166, 298, 388, 424, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 1687, 166, 184, 316, 424, 460, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pd(pbuffer, 1795, 184, 202, 334, 460, 496, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 1903, 298, 316, 604, 634, 664, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 1993, 316, 334, 614, 664, 694, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 2083, 388, 424, 634, 754, 814, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 2263, 424, 460, 664, 814, 874, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 2443, 460, 496, 694, 874, 934, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 2623, 604, 614, 1054, 1069, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 2668, 634, 664, 1054, 1084, 1129, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 2803, 664, 694, 1069, 1129, 1174, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pg(pbuffer, 2938, 754, 814, 1084, 1219, 1309, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pg(pbuffer, 3208, 814, 874, 1129, 1309, 1399, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pg(pbuffer, 3478, 874, 934, 1174, 1399, 1489, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_df(pbuffer, 3748, 754, 814, 1579, 1687, 1903, 2083, 2263, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_df(pbuffer, 4108, 814, 874, 1687, 1795, 1993, 2263, 2443, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dg(pbuffer, 4468, 1084, 1129, 1903, 1993, 2623, 2668, 2803, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dg(pbuffer, 4738, 1219, 1309, 2083, 2263, 2668, 2938, 3208, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dg(pbuffer, 5278, 1309, 1399, 2263, 2443, 2803, 3208, 3478, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fg(pbuffer, 5818, 2938, 3208, 3748, 4108, 4468, 4738, 5278, factors, 11, 17, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 5818, quadrupoles, 6, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<3, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // npotrec namespace

#endif /* NuclearPotentialGeom020SumRecFG_hpp */
