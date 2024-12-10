#ifndef NuclearPotentialGeom200SumRecGP_hpp
#define NuclearPotentialGeom200SumRecGP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "GeometricalDerivatives2X0ForGY.hpp"
#include "GtoBlock.hpp"
#include "NuclearPotentialPrimRecDP.hpp"
#include "NuclearPotentialPrimRecDS.hpp"
#include "NuclearPotentialPrimRecFP.hpp"
#include "NuclearPotentialPrimRecFS.hpp"
#include "NuclearPotentialPrimRecGP.hpp"
#include "NuclearPotentialPrimRecGS.hpp"
#include "NuclearPotentialPrimRecHP.hpp"
#include "NuclearPotentialPrimRecHS.hpp"
#include "NuclearPotentialPrimRecIP.hpp"
#include "NuclearPotentialPrimRecPP.hpp"
#include "NuclearPotentialPrimRecPS.hpp"
#include "NuclearPotentialPrimRecSP.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes (d^(2)/dA^(2)G|A|P)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_20_gp(T&                               distributor,
                                      const CGtoBlock&                 bra_gto_block,
                                      const CGtoBlock&                 ket_gto_block,
                                      const std::pair<size_t, size_t>& bra_indices,
                                      const std::pair<size_t, size_t>& ket_indices,
                                      const bool                       bra_eq_ket) -> void
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

    CSimdArray<double> pbuffer(1084, ket_npgtos);

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

                t2cfunc::comp_distances_pa_from_p(factors, 11, 8, r_a);

                t2cfunc::comp_distances_pb_from_p(factors, 14, 8, 2);

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

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 9, 1, 2, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 12, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 15, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 18, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 21, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 24, 6, 7, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 27, 7, 8, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 30, 1, 2, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 33, 2, 3, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 36, 3, 4, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 39, 4, 5, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 42, 5, 6, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 45, 6, 7, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 48, 1, 2, 9, 12, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 57, 2, 3, 12, 15, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 66, 3, 4, 15, 18, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 75, 4, 5, 18, 21, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 84, 5, 6, 21, 24, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 93, 6, 7, 24, 27, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 102, 1, 2, 30, 33, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 108, 2, 3, 33, 36, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 114, 3, 4, 36, 39, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 120, 4, 5, 39, 42, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 126, 5, 6, 42, 45, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 132, 9, 12, 30, 33, 48, 57, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 150, 12, 15, 33, 36, 57, 66, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 168, 15, 18, 36, 39, 66, 75, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 186, 18, 21, 39, 42, 75, 84, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 204, 21, 24, 42, 45, 84, 93, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 222, 30, 33, 102, 108, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 232, 33, 36, 108, 114, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 242, 36, 39, 114, 120, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fs(pbuffer, 252, 39, 42, 120, 126, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 262, 48, 57, 102, 108, 132, 150, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 292, 57, 66, 108, 114, 150, 168, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 322, 66, 75, 114, 120, 168, 186, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 352, 75, 84, 120, 126, 186, 204, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gs(pbuffer, 382, 102, 108, 222, 232, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gs(pbuffer, 397, 108, 114, 232, 242, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gs(pbuffer, 412, 114, 120, 242, 252, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gp(pbuffer, 427, 132, 150, 222, 232, 262, 292, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gp(pbuffer, 472, 150, 168, 232, 242, 292, 322, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gp(pbuffer, 517, 168, 186, 242, 252, 322, 352, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hs(pbuffer, 562, 222, 232, 382, 397, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hs(pbuffer, 583, 232, 242, 397, 412, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hp(pbuffer, 604, 262, 292, 382, 397, 427, 472, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_hp(pbuffer, 667, 292, 322, 397, 412, 472, 517, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ip(pbuffer, 730, 427, 472, 562, 583, 604, 667, factors, 11, 17, a_exp);

                    t2cgeom::comp_prim_op_geom_20_gx(pbuffer, 814, 132, 427, 730, 1, 3, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 814, charges[l], ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 1>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 1, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace npotrec

#endif /* NuclearPotentialGeom200SumRecGP_hpp */
