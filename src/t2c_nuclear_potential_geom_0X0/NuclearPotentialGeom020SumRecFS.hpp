#ifndef NuclearPotentialGeom020SumRecFS_hpp
#define NuclearPotentialGeom020SumRecFS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "GtoBlock.hpp"
#include "NuclearPotentialGeom010PrimRecDS.hpp"
#include "NuclearPotentialGeom010PrimRecPS.hpp"
#include "NuclearPotentialGeom010PrimRecSS.hpp"
#include "NuclearPotentialGeom020PrimRecDS.hpp"
#include "NuclearPotentialGeom020PrimRecFS.hpp"
#include "NuclearPotentialGeom020PrimRecPS.hpp"
#include "NuclearPotentialGeom020PrimRecSS.hpp"
#include "NuclearPotentialPrimRecPS.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes (F|AG(2)|S)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_020_fs(T&                               distributor,
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

    CSimdArray<double> factors(17, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(264, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(60, 1);

    CSimdArray<double> sbuffer(42, 1);

    // setup Boys function data

    const CBoysFunc<5> bf_table;

    CSimdArray<double> bf_data(7, ket_npgtos);

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

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                for (size_t l = 0; l < coords.size(); l++)
                {
                    t2cfunc::comp_distances_pc(factors, 14, 8, coords[l]);

                    t2cfunc::comp_boys_args(bf_data, 6, factors, 14, a_exp);

                    bf_table.compute(bf_data, 0, 6);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 1, 0, bf_data, 1, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 2, 0, bf_data, 2, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 3, 0, bf_data, 3, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 4, 0, bf_data, 4, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 5, 0, bf_data, 5, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 6, 2, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 9, 3, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 12, 4, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 15, 1, 2, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 21, 2, 3, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 27, 3, 4, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, 33, 4, 5, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 39, 2, 3, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 42, 2, 6, 9, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 51, 3, 9, 12, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 60, 6, 15, 21, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 78, 9, 21, 27, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_geom_020_ps(pbuffer, 96, 12, 27, 33, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 114, 6, 9, 39, 42, 51, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ds(pbuffer, 132, 15, 21, 42, 60, 78, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_ds(pbuffer, 168, 21, 27, 51, 78, 96, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_020_fs(pbuffer, 204, 60, 78, 114, 132, 168, factors, 11, 14, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 204, quadrupoles, 6, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<3, 0>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 0, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace npotrec

#endif /* NuclearPotentialGeom020SumRecFS_hpp */
