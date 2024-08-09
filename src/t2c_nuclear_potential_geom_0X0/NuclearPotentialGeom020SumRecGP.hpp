#ifndef NuclearPotentialGeom020SumRecGP_hpp
#define NuclearPotentialGeom020SumRecGP_hpp

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
#include "NuclearPotentialPrimRecPS.hpp"
#include "NuclearPotentialGeom010PrimRecPS.hpp"
#include "NuclearPotentialGeom020PrimRecPS.hpp"
#include "NuclearPotentialPrimRecPP.hpp"
#include "NuclearPotentialGeom010PrimRecPP.hpp"
#include "NuclearPotentialGeom020PrimRecPP.hpp"
#include "NuclearPotentialGeom010PrimRecDS.hpp"
#include "NuclearPotentialGeom020PrimRecDS.hpp"
#include "NuclearPotentialPrimRecDP.hpp"
#include "NuclearPotentialGeom010PrimRecDP.hpp"
#include "NuclearPotentialGeom020PrimRecDP.hpp"
#include "NuclearPotentialGeom020PrimRecFS.hpp"
#include "NuclearPotentialGeom010PrimRecFP.hpp"
#include "NuclearPotentialGeom020PrimRecFP.hpp"
#include "NuclearPotentialGeom020PrimRecGP.hpp"
#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes (G|AG(2)|P)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_020_gp(T& distributor,
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

    CSimdArray<double> pbuffer(2048, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(270, 1);

    CSimdArray<double> sbuffer(162, 1);

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

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 1, 0, bf_data, 1, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 2, 0, bf_data, 2, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 3, 0, bf_data, 3, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 4, 0, bf_data, 4, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 5, 0, bf_data, 5, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 6, 0, bf_data, 6, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 7, 0, bf_data, 7, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 8, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 11, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 14, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 17, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 20, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 23, 1, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 29, 2, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 35, 3, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 41, 4, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 47, 5, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 53, 6, 7, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 59, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 62, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 65, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 68, 2, 8, 11, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 77, 3, 11, 14, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 86, 4, 14, 17, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 95, 5, 17, 20, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 104, 8, 23, 29, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 122, 11, 29, 35, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 140, 14, 35, 41, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 158, 17, 41, 47, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_sp(pbuffer, 176, 20, 47, 53, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 194, 2, 3, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 197, 3, 4, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 200, 2, 8, 11, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 209, 3, 11, 14, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 218, 4, 14, 17, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 227, 8, 23, 29, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 245, 11, 29, 35, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 263, 14, 35, 41, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 281, 17, 41, 47, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 299, 2, 3, 59, 62, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 308, 3, 4, 62, 65, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 317, 8, 11, 59, 68, 77, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 344, 11, 14, 62, 77, 86, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 371, 14, 17, 65, 86, 95, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 398, 23, 29, 68, 104, 122, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 452, 29, 35, 77, 122, 140, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 506, 35, 41, 86, 140, 158, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_pp(pbuffer, 560, 41, 47, 95, 158, 176, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 614, 8, 11, 194, 200, 209, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 632, 11, 14, 197, 209, 218, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ds(pbuffer, 650, 23, 29, 200, 227, 245, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ds(pbuffer, 686, 29, 35, 209, 245, 263, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ds(pbuffer, 722, 35, 41, 218, 263, 281, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 758, 59, 62, 194, 197, 299, 308, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(pbuffer, 776, 68, 77, 200, 209, 299, 317, 344, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(pbuffer, 830, 77, 86, 209, 218, 308, 344, 371, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dp(pbuffer, 884, 104, 122, 227, 245, 317, 398, 452, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dp(pbuffer, 992, 122, 140, 245, 263, 344, 452, 506, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_dp(pbuffer, 1100, 140, 158, 263, 281, 371, 506, 560, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fs(pbuffer, 1208, 227, 245, 614, 650, 686, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fs(pbuffer, 1268, 245, 263, 632, 686, 722, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fp(pbuffer, 1328, 317, 344, 614, 632, 758, 776, 830, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fp(pbuffer, 1418, 398, 452, 650, 686, 776, 884, 992, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fp(pbuffer, 1598, 452, 506, 686, 722, 830, 992, 1100, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_gp(pbuffer, 1778, 884, 992, 1208, 1268, 1328, 1418, 1598, factors, 11, 17, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 1778, quadrupoles, 6, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 1>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 1, j, ket_range, bra_eq_ket);
        }
    }
}

} // npotrec namespace

#endif /* NuclearPotentialGeom020SumRecGP_hpp */
