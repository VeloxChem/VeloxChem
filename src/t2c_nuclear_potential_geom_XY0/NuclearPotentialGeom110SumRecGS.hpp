#ifndef NuclearPotentialGeom110SumRecGS_hpp
#define NuclearPotentialGeom110SumRecGS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "GeometricalDerivatives1X0ForGY.hpp"
#include "GtoBlock.hpp"
#include "NuclearPotentialGeom010PrimRecDS.hpp"
#include "NuclearPotentialGeom010PrimRecFS.hpp"
#include "NuclearPotentialGeom010PrimRecGS.hpp"
#include "NuclearPotentialGeom010PrimRecHS.hpp"
#include "NuclearPotentialGeom010PrimRecPS.hpp"
#include "NuclearPotentialGeom010PrimRecSS.hpp"
#include "NuclearPotentialPrimRecDS.hpp"
#include "NuclearPotentialPrimRecFS.hpp"
#include "NuclearPotentialPrimRecGS.hpp"
#include "NuclearPotentialPrimRecPS.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes (d^(1)/dA^(1)G|AG(1)|S)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_110_gs(T&                               distributor,
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

    CSimdArray<double> factors(17, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(585, ket_npgtos);

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

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 6, 0, bf_data, 6, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 7, 1, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 10, 2, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 13, 3, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 16, 4, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 19, 5, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 22, 6, factors, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 25, 1, 2, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 28, 2, 3, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 31, 3, 4, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 34, 4, 5, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 37, 1, 7, 10, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 46, 2, 10, 13, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 55, 3, 13, 16, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 64, 4, 16, 19, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_geom_010_ps(pbuffer, 73, 5, 19, 22, factors, 11, 14);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 82, 1, 2, 25, 28, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 88, 2, 3, 28, 31, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 94, 3, 4, 31, 34, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 100, 7, 10, 25, 37, 46, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 118, 10, 13, 28, 46, 55, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 136, 13, 16, 31, 55, 64, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ds(pbuffer, 154, 16, 19, 34, 64, 73, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 172, 25, 28, 82, 88, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 182, 28, 31, 88, 94, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fs(pbuffer, 192, 37, 46, 82, 100, 118, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fs(pbuffer, 222, 46, 55, 88, 118, 136, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fs(pbuffer, 252, 55, 64, 94, 136, 154, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_gs(pbuffer, 282, 82, 88, 172, 182, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_gs(pbuffer, 297, 100, 118, 172, 192, 222, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_gs(pbuffer, 342, 118, 136, 182, 222, 252, factors, 11, 14, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_hs(pbuffer, 387, 192, 222, 282, 297, 342, factors, 11, 14, a_exp);

                    t2cgeom::comp_prim_op_geom_10_gx(pbuffer, 450, 192, 387, 3, 1, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 450, dipoles, 3, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 0>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 0, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace npotrec

#endif /* NuclearPotentialGeom110SumRecGS_hpp */
