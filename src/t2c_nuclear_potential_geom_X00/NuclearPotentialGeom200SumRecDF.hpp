#ifndef NuclearPotentialGeom200SumRecDF_hpp
#define NuclearPotentialGeom200SumRecDF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "GeometricalDerivatives2X0ForDY.hpp"
#include "GtoBlock.hpp"
#include "NuclearPotentialPrimRecDD.hpp"
#include "NuclearPotentialPrimRecDF.hpp"
#include "NuclearPotentialPrimRecDP.hpp"
#include "NuclearPotentialPrimRecFD.hpp"
#include "NuclearPotentialPrimRecFF.hpp"
#include "NuclearPotentialPrimRecGF.hpp"
#include "NuclearPotentialPrimRecPD.hpp"
#include "NuclearPotentialPrimRecPF.hpp"
#include "NuclearPotentialPrimRecPP.hpp"
#include "NuclearPotentialPrimRecPS.hpp"
#include "NuclearPotentialPrimRecSD.hpp"
#include "NuclearPotentialPrimRecSF.hpp"
#include "NuclearPotentialPrimRecSP.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes (d^(2)/dA^(2)D|A|F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_20_df(T&                               distributor,
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

    CSimdArray<double> pbuffer(1528, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(360, 1);

    CSimdArray<double> sbuffer(210, 1);

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

                t2cfunc::comp_distances_pb_from_p(factors, 14, 8, 2);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                for (size_t l = 0; l < coords.size(); l++)
                {
                    t2cfunc::comp_distances_pc(factors, 17, 8, coords[l]);

                    t2cfunc::comp_boys_args(bf_data, 6, factors, 17, a_exp);

                    bf_table.compute(bf_data, 0, 6);

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

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 30, 1, 2, 9, 12, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 36, 2, 3, 12, 15, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 42, 3, 4, 15, 18, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 48, 4, 5, 18, 21, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 54, 5, 6, 21, 24, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 60, 6, 7, 24, 27, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 66, 9, 12, 30, 36, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 76, 12, 15, 36, 42, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 86, 15, 18, 42, 48, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 96, 18, 21, 48, 54, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 106, 21, 24, 54, 60, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 116, 1, 2, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 119, 2, 3, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 122, 3, 4, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 125, 4, 5, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 128, 1, 2, 9, 12, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 137, 2, 3, 12, 15, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 146, 3, 4, 15, 18, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 155, 4, 5, 18, 21, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 164, 9, 12, 30, 36, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 182, 12, 15, 36, 42, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 200, 15, 18, 42, 48, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 218, 18, 21, 48, 54, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 236, 30, 36, 66, 76, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 266, 36, 42, 76, 86, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 296, 42, 48, 86, 96, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 326, 48, 54, 96, 106, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 356, 9, 12, 116, 119, 128, 137, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 374, 12, 15, 119, 122, 137, 146, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 392, 15, 18, 122, 125, 146, 155, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 410, 30, 36, 128, 137, 164, 182, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 446, 36, 42, 137, 146, 182, 200, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 482, 42, 48, 146, 155, 200, 218, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 518, 66, 76, 164, 182, 236, 266, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 578, 76, 86, 182, 200, 266, 296, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 638, 86, 96, 200, 218, 296, 326, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 698, 164, 182, 356, 374, 410, 446, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 758, 182, 200, 374, 392, 446, 482, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 818, 236, 266, 410, 446, 518, 578, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 918, 266, 296, 446, 482, 578, 638, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gf(pbuffer, 1018, 518, 578, 698, 758, 818, 918, factors, 11, 17, a_exp);

                    t2cgeom::comp_prim_op_geom_20_dx(pbuffer, 1168, 66, 518, 1018, 1, 10, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 1168, charges[l], ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<2, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 2, 3, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace npotrec

#endif /* NuclearPotentialGeom200SumRecDF_hpp */
