#ifndef NuclearPotentialGeom010SumErfRecGP_hpp
#define NuclearPotentialGeom010SumErfRecGP_hpp

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
#include "NuclearPotentialPrimRecPS.hpp"
#include "NuclearPotentialGeom010PrimRecPS.hpp"
#include "NuclearPotentialPrimRecPP.hpp"
#include "NuclearPotentialGeom010PrimRecPP.hpp"
#include "NuclearPotentialPrimRecDS.hpp"
#include "NuclearPotentialGeom010PrimRecDS.hpp"
#include "NuclearPotentialPrimRecDP.hpp"
#include "NuclearPotentialGeom010PrimRecDP.hpp"
#include "NuclearPotentialGeom010PrimRecFS.hpp"
#include "NuclearPotentialPrimRecFP.hpp"
#include "NuclearPotentialGeom010PrimRecFP.hpp"
#include "NuclearPotentialGeom010PrimRecGP.hpp"
#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes (G|Erf(AG(1))|P)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param omegas The vector of range-separation factors.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_erf_nuclear_potential_geom_010_gp(T& distributor,
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

    CSimdArray<double> pbuffer(931, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(135, 1);

    CSimdArray<double> sbuffer(81, 1);

    // setup Boys function data

    const CBoysFunc<6> bf_table;

    CSimdArray<double> bf_data(8, ket_npgtos);

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

                    t2cfunc::comp_boys_args(bf_data, 7, factors, 17, a_exp, omegas[l]);

                    bf_table.compute(bf_data, 0, 7, factors, a_exp, omegas[l]);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 1, 0, bf_data, 1, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 2, 0, bf_data, 2, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 3, 0, bf_data, 3, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 4, 0, bf_data, 4, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 5, 0, bf_data, 5, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 6, 0, bf_data, 6, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 7, 1, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 10, 2, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 13, 3, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 16, 4, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 19, 5, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 22, 6, factors, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 25, 1, 2, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 28, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 31, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 34, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 37, 1, 7, 10, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 46, 2, 10, 13, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 55, 3, 13, 16, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 64, 4, 16, 19, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 73, 5, 19, 22, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 82, 1, 2, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 85, 2, 3, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 88, 3, 4, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 91, 1, 7, 10, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 100, 2, 10, 13, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 109, 3, 13, 16, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 118, 4, 16, 19, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 127, 1, 2, 25, 28, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 136, 2, 3, 28, 31, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 145, 3, 4, 31, 34, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 154, 7, 10, 25, 37, 46, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 181, 10, 13, 28, 46, 55, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 208, 13, 16, 31, 55, 64, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pp(pbuffer, 235, 16, 19, 34, 64, 73, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 262, 1, 2, 82, 85, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 268, 2, 3, 85, 88, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 274, 7, 10, 82, 91, 100, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 292, 10, 13, 85, 100, 109, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 310, 13, 16, 88, 109, 118, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 328, 25, 28, 82, 85, 127, 136, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 346, 28, 31, 85, 88, 136, 145, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(pbuffer, 364, 37, 46, 91, 100, 127, 154, 181, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(pbuffer, 418, 46, 55, 100, 109, 136, 181, 208, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dp(pbuffer, 472, 55, 64, 109, 118, 145, 208, 235, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fs(pbuffer, 526, 91, 100, 262, 274, 292, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fs(pbuffer, 556, 100, 109, 268, 292, 310, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 586, 127, 136, 262, 268, 328, 346, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fp(pbuffer, 616, 154, 181, 274, 292, 328, 364, 418, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fp(pbuffer, 706, 181, 208, 292, 310, 346, 418, 472, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_gp(pbuffer, 796, 364, 418, 526, 556, 586, 616, 706, factors, 11, 17, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 796, dipoles, 3, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 1>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 1, j, ket_range, bra_eq_ket);
        }
    }
}

} // npotrec namespace

#endif /* NuclearPotentialGeom010SumErfRecGP_hpp */
