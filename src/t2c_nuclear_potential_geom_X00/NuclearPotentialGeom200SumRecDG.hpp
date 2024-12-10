#ifndef NuclearPotentialGeom200SumRecDG_hpp
#define NuclearPotentialGeom200SumRecDG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "GeometricalDerivatives2X0ForDY.hpp"
#include "GtoBlock.hpp"
#include "NuclearPotentialPrimRecDD.hpp"
#include "NuclearPotentialPrimRecDF.hpp"
#include "NuclearPotentialPrimRecDG.hpp"
#include "NuclearPotentialPrimRecFF.hpp"
#include "NuclearPotentialPrimRecFG.hpp"
#include "NuclearPotentialPrimRecGG.hpp"
#include "NuclearPotentialPrimRecPD.hpp"
#include "NuclearPotentialPrimRecPF.hpp"
#include "NuclearPotentialPrimRecPG.hpp"
#include "NuclearPotentialPrimRecPP.hpp"
#include "NuclearPotentialPrimRecSD.hpp"
#include "NuclearPotentialPrimRecSF.hpp"
#include "NuclearPotentialPrimRecSG.hpp"
#include "NuclearPotentialPrimRecSP.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"

namespace npotrec {  // npotrec namespace

/// @brief Computes (d^(2)/dA^(2)D|A|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_20_dg(T&                               distributor,
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

    CSimdArray<double> pbuffer(2442, ket_npgtos);

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

                t2cfunc::comp_distances_pa_from_p(factors, 11, 8, r_a);

                t2cfunc::comp_distances_pb_from_p(factors, 14, 8, 2);

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

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 10, 1, 2, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 13, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 16, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 19, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 22, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 25, 6, 7, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 28, 7, 8, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 31, 8, 9, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 34, 1, 2, 10, 13, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 40, 2, 3, 13, 16, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 46, 3, 4, 16, 19, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 52, 4, 5, 19, 22, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 58, 5, 6, 22, 25, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 64, 6, 7, 25, 28, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 70, 7, 8, 28, 31, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 76, 10, 13, 34, 40, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 86, 13, 16, 40, 46, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 96, 16, 19, 46, 52, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 106, 19, 22, 52, 58, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 116, 22, 25, 58, 64, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 126, 25, 28, 64, 70, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 136, 34, 40, 76, 86, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 151, 40, 46, 86, 96, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 166, 46, 52, 96, 106, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 181, 52, 58, 106, 116, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 196, 58, 64, 116, 126, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 211, 1, 2, 10, 13, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 220, 2, 3, 13, 16, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 229, 3, 4, 16, 19, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 238, 4, 5, 19, 22, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 247, 10, 13, 34, 40, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 265, 13, 16, 40, 46, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 283, 16, 19, 46, 52, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 301, 19, 22, 52, 58, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 319, 34, 40, 76, 86, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 349, 40, 46, 86, 96, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 379, 46, 52, 96, 106, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 409, 52, 58, 106, 116, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 439, 76, 86, 136, 151, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 484, 86, 96, 151, 166, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 529, 96, 106, 166, 181, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 574, 106, 116, 181, 196, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 619, 34, 40, 211, 220, 247, 265, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 655, 40, 46, 220, 229, 265, 283, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 691, 46, 52, 229, 238, 283, 301, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 727, 76, 86, 247, 265, 319, 349, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 787, 86, 96, 265, 283, 349, 379, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 847, 96, 106, 283, 301, 379, 409, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 907, 136, 151, 319, 349, 439, 484, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 997, 151, 166, 349, 379, 484, 529, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 1087, 166, 181, 379, 409, 529, 574, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 1177, 319, 349, 619, 655, 727, 787, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 1277, 349, 379, 655, 691, 787, 847, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fg(pbuffer, 1377, 439, 484, 727, 787, 907, 997, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fg(pbuffer, 1527, 484, 529, 787, 847, 997, 1087, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gg(pbuffer, 1677, 907, 997, 1177, 1277, 1377, 1527, factors, 11, 17, a_exp);

                    t2cgeom::comp_prim_op_geom_20_dx(pbuffer, 1902, 136, 907, 1677, 1, 15, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 1902, charges[l], ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<2, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 2, 4, j, ket_range, bra_eq_ket);
        }
    }
}

}  // namespace npotrec

#endif /* NuclearPotentialGeom200SumRecDG_hpp */
