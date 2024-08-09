#ifndef NuclearPotentialGeom020SumRecDG_hpp
#define NuclearPotentialGeom020SumRecDG_hpp

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
#include "NuclearPotentialGeom020PrimRecPF.hpp"
#include "NuclearPotentialGeom010PrimRecPG.hpp"
#include "NuclearPotentialGeom020PrimRecPG.hpp"
#include "NuclearPotentialGeom020PrimRecDG.hpp"
#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes (D|AG(2)|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_020_dg(T& distributor,
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

    CSimdArray<double> pbuffer(2804, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(540, 1);

    CSimdArray<double> sbuffer(270, 1);

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

                t2cfunc::comp_distances_pa_from_p(factors, 11 , 8, r_a);

                t2cfunc::comp_distances_pb_from_p(factors, 14 , 8, 2);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                for (size_t l = 0; l < coords.size(); l++)
                {
                    t2cfunc::comp_distances_pc(factors, 17, 8, coords[l]);

                    t2cfunc::comp_boys_args(bf_data, 9, factors, 17, a_exp);

                    bf_table.compute(bf_data, 0, 9);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 1, 0, bf_data, 1, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 2, 0, bf_data, 2, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 3, 0, bf_data, 3, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 4, 0, bf_data, 4, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 5, 0, bf_data, 5, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 6, 0, bf_data, 6, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 7, 0, bf_data, 7, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 8, 0, bf_data, 8, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 9, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 12, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 15, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 18, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 21, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 24, 7, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 27, 1, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 33, 2, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 39, 3, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 45, 4, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 51, 5, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 57, 6, 7, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 63, 7, 8, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 69, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 72, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 75, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 78, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 81, 2, 9, 12, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 90, 3, 12, 15, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 99, 4, 15, 18, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 108, 5, 18, 21, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 117, 6, 21, 24, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 126, 9, 27, 33, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 144, 12, 33, 39, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 162, 15, 39, 45, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 180, 18, 45, 51, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 198, 21, 51, 57, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 216, 24, 57, 63, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 234, 2, 3, 69, 72, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 240, 3, 4, 72, 75, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 246, 4, 5, 75, 78, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 252, 9, 12, 69, 81, 90, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 270, 12, 15, 72, 90, 99, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 288, 15, 18, 75, 99, 108, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 306, 18, 21, 78, 108, 117, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 324, 27, 33, 81, 126, 144, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 360, 33, 39, 90, 144, 162, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 396, 39, 45, 99, 162, 180, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 432, 45, 51, 108, 180, 198, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sd(pbuffer, 468, 51, 57, 117, 198, 216, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 504, 69, 72, 234, 240, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 514, 72, 75, 240, 246, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 524, 81, 90, 234, 252, 270, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 554, 90, 99, 240, 270, 288, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 584, 99, 108, 246, 288, 306, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 614, 126, 144, 252, 324, 360, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 674, 144, 162, 270, 360, 396, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 734, 162, 180, 288, 396, 432, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sf(pbuffer, 794, 180, 198, 306, 432, 468, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 854, 234, 240, 504, 514, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 869, 252, 270, 504, 524, 554, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 914, 270, 288, 514, 554, 584, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 959, 324, 360, 524, 614, 674, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 1049, 360, 396, 554, 674, 734, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_sg(pbuffer, 1139, 396, 432, 584, 734, 794, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 1229, 324, 360, 524, 614, 674, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pf(pbuffer, 1409, 360, 396, 554, 674, 734, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 1589, 524, 554, 854, 869, 914, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pg(pbuffer, 1724, 614, 674, 869, 959, 1049, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pg(pbuffer, 1994, 674, 734, 914, 1049, 1139, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dg(pbuffer, 2264, 959, 1049, 1229, 1409, 1589, 1724, 1994, factors, 11, 17, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 2264, quadrupoles, 6, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<2, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 2, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // npotrec namespace

#endif /* NuclearPotentialGeom020SumRecDG_hpp */
