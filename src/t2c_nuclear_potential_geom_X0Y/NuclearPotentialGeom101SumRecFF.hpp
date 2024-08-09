#ifndef NuclearPotentialGeom101SumRecFF_hpp
#define NuclearPotentialGeom101SumRecFF_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "OverlapPrimRecSS.hpp"
#include "NuclearPotentialPrimRecSS.hpp"
#include "NuclearPotentialPrimRecSP.hpp"
#include "NuclearPotentialPrimRecSD.hpp"
#include "NuclearPotentialPrimRecSF.hpp"
#include "NuclearPotentialPrimRecSG.hpp"
#include "NuclearPotentialPrimRecPS.hpp"
#include "NuclearPotentialPrimRecPP.hpp"
#include "NuclearPotentialPrimRecPD.hpp"
#include "NuclearPotentialPrimRecPF.hpp"
#include "NuclearPotentialPrimRecPG.hpp"
#include "NuclearPotentialPrimRecDS.hpp"
#include "NuclearPotentialPrimRecDP.hpp"
#include "NuclearPotentialPrimRecDD.hpp"
#include "NuclearPotentialPrimRecDF.hpp"
#include "NuclearPotentialPrimRecDG.hpp"
#include "NuclearPotentialPrimRecFP.hpp"
#include "NuclearPotentialPrimRecFD.hpp"
#include "NuclearPotentialPrimRecFF.hpp"
#include "NuclearPotentialPrimRecFG.hpp"
#include "NuclearPotentialPrimRecGD.hpp"
#include "NuclearPotentialPrimRecGG.hpp"
#include "GeometricalDerivatives1X1ForFF.hpp"

#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes (d^(1)/dA^(1)F|A|d^(1)/dB^(1)F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_nuclear_potential_geom_11_ff(T& distributor,
                                      const CGtoBlock& bra_gto_block,
                                      const CGtoBlock& ket_gto_block,
                                      const std::pair<size_t, size_t>& bra_indices,
                                      const std::pair<size_t, size_t>& ket_indices,
                                      const bool bra_eq_ket) -> void
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

    CSimdArray<double> pbuffer(3156, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(900, 1);

    CSimdArray<double> sbuffer(441, 1);

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

                    t2cfunc::comp_boys_args(bf_data, 7, factors, 17, a_exp);

                    bf_table.compute(bf_data, 0, 7);

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

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 211, 1, 2, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 214, 2, 3, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 217, 3, 4, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_ps(pbuffer, 220, 4, 5, factors, 11, 17);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 223, 1, 2, 10, 13, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 232, 2, 3, 13, 16, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 241, 3, 4, 16, 19, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pp(pbuffer, 250, 4, 5, 19, 22, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 259, 10, 13, 34, 40, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 277, 13, 16, 40, 46, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 295, 16, 19, 46, 52, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pd(pbuffer, 313, 19, 22, 52, 58, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 331, 34, 40, 76, 86, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 361, 40, 46, 86, 96, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 391, 46, 52, 96, 106, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 421, 52, 58, 106, 116, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 451, 76, 86, 136, 151, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 496, 86, 96, 151, 166, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 541, 96, 106, 166, 181, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 586, 106, 116, 181, 196, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 631, 1, 2, 211, 214, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 637, 2, 3, 214, 217, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ds(pbuffer, 643, 3, 4, 217, 220, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 649, 10, 13, 211, 214, 223, 232, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 667, 13, 16, 214, 217, 232, 241, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dp(pbuffer, 685, 16, 19, 217, 220, 241, 250, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 703, 34, 40, 223, 232, 259, 277, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 739, 40, 46, 232, 241, 277, 295, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dd(pbuffer, 775, 46, 52, 241, 250, 295, 313, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 811, 76, 86, 259, 277, 331, 361, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 871, 86, 96, 277, 295, 361, 391, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_df(pbuffer, 931, 96, 106, 295, 313, 391, 421, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 991, 136, 151, 331, 361, 451, 496, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 1081, 151, 166, 361, 391, 496, 541, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 1171, 166, 181, 391, 421, 541, 586, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 1261, 223, 232, 631, 637, 649, 667, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fp(pbuffer, 1291, 232, 241, 637, 643, 667, 685, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 1321, 259, 277, 649, 667, 703, 739, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fd(pbuffer, 1381, 277, 295, 667, 685, 739, 775, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 1441, 331, 361, 703, 739, 811, 871, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_ff(pbuffer, 1541, 361, 391, 739, 775, 871, 931, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fg(pbuffer, 1641, 451, 496, 811, 871, 991, 1081, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_fg(pbuffer, 1791, 496, 541, 871, 931, 1081, 1171, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gd(pbuffer, 1941, 703, 739, 1261, 1291, 1321, 1381, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_gg(pbuffer, 2031, 991, 1081, 1441, 1541, 1641, 1791, factors, 11, 17, a_exp);

                    t2cgeom::comp_prim_op_geom_11_ff(pbuffer, 2256, 703, 991, 1941, 2031, 1, factors, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 0, charges[l], ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<3, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 3, j, ket_range, bra_eq_ket);
        }
    }
}

} // npotrec namespace

#endif /* NuclearPotentialGeom101SumRecFF_hpp */
