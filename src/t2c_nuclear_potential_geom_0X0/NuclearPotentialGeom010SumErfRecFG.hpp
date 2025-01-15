#ifndef NuclearPotentialGeom010SumErfRecFG_hpp
#define NuclearPotentialGeom010SumErfRecFG_hpp

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
#include "NuclearPotentialPrimRecSD.hpp"
#include "NuclearPotentialGeom010PrimRecSD.hpp"
#include "NuclearPotentialPrimRecSF.hpp"
#include "NuclearPotentialGeom010PrimRecSF.hpp"
#include "NuclearPotentialPrimRecSG.hpp"
#include "NuclearPotentialGeom010PrimRecSG.hpp"
#include "NuclearPotentialGeom010PrimRecPD.hpp"
#include "NuclearPotentialPrimRecPF.hpp"
#include "NuclearPotentialGeom010PrimRecPF.hpp"
#include "NuclearPotentialPrimRecPG.hpp"
#include "NuclearPotentialGeom010PrimRecPG.hpp"
#include "NuclearPotentialGeom010PrimRecDF.hpp"
#include "NuclearPotentialPrimRecDG.hpp"
#include "NuclearPotentialGeom010PrimRecDG.hpp"
#include "NuclearPotentialGeom010PrimRecFG.hpp"
#include "BoysFunc.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes (F|Erf(AG(1))|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param omegas The vector of range-separation factors.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_erf_nuclear_potential_geom_010_fg(T& distributor,
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

    CSimdArray<double> pbuffer(3094, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(450, 1);

    CSimdArray<double> sbuffer(189, 1);

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

                    t2cfunc::comp_boys_args(bf_data, 9, factors, 17, a_exp, omegas[l]);

                    bf_table.compute(bf_data, 0, 9, factors, a_exp, omegas[l]);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 1, 0, bf_data, 1, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 2, 0, bf_data, 2, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 3, 0, bf_data, 3, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 4, 0, bf_data, 4, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 5, 0, bf_data, 5, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 6, 0, bf_data, 6, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 7, 0, bf_data, 7, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_ss(pbuffer, 8, 0, bf_data, 8, factors, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 9, 1, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 12, 2, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 15, 3, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 18, 4, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 21, 5, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 24, 6, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 27, 7, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, 30, 8, factors, 17, omegas[l], a_exp);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 33, 1, 2, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 36, 2, 3, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 39, 3, 4, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 42, 4, 5, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 45, 5, 6, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sp(pbuffer, 48, 6, 7, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 51, 1, 9, 12, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 60, 2, 12, 15, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 69, 3, 15, 18, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 78, 4, 18, 21, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 87, 5, 21, 24, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 96, 6, 24, 27, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_geom_010_sp(pbuffer, 105, 7, 27, 30, factors, 14, 17);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 114, 1, 2, 33, 36, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 120, 2, 3, 36, 39, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 126, 3, 4, 39, 42, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 132, 4, 5, 42, 45, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sd(pbuffer, 138, 5, 6, 45, 48, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 144, 9, 12, 33, 51, 60, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 162, 12, 15, 36, 60, 69, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 180, 15, 18, 39, 69, 78, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 198, 18, 21, 42, 78, 87, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 216, 21, 24, 45, 87, 96, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sd(pbuffer, 234, 24, 27, 48, 96, 105, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 252, 33, 36, 114, 120, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 262, 36, 39, 120, 126, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 272, 39, 42, 126, 132, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sf(pbuffer, 282, 42, 45, 132, 138, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 292, 51, 60, 114, 144, 162, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 322, 60, 69, 120, 162, 180, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 352, 69, 78, 126, 180, 198, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 382, 78, 87, 132, 198, 216, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sf(pbuffer, 412, 87, 96, 138, 216, 234, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 442, 114, 120, 252, 262, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 457, 120, 126, 262, 272, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_sg(pbuffer, 472, 126, 132, 272, 282, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 487, 144, 162, 252, 292, 322, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 532, 162, 180, 262, 322, 352, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 577, 180, 198, 272, 352, 382, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_sg(pbuffer, 622, 198, 216, 282, 382, 412, factors, 14, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 667, 51, 60, 114, 144, 162, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 721, 60, 69, 120, 162, 180, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pd(pbuffer, 775, 69, 78, 126, 180, 198, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 829, 114, 120, 252, 262, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pf(pbuffer, 859, 120, 126, 262, 272, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 889, 144, 162, 252, 292, 322, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 979, 162, 180, 262, 322, 352, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pf(pbuffer, 1069, 180, 198, 272, 352, 382, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 1159, 252, 262, 442, 457, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_pg(pbuffer, 1204, 262, 272, 457, 472, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 1249, 292, 322, 442, 487, 532, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 1384, 322, 352, 457, 532, 577, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_pg(pbuffer, 1519, 352, 382, 472, 577, 622, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_df(pbuffer, 1654, 292, 322, 667, 721, 829, 889, 979, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_df(pbuffer, 1834, 322, 352, 721, 775, 859, 979, 1069, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_dg(pbuffer, 2014, 442, 457, 829, 859, 1159, 1204, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dg(pbuffer, 2104, 487, 532, 889, 979, 1159, 1249, 1384, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_dg(pbuffer, 2374, 532, 577, 979, 1069, 1204, 1384, 1519, factors, 11, 17, a_exp);

                    npotrec::comp_prim_nuclear_potential_geom_010_fg(pbuffer, 2644, 1249, 1384, 1654, 1834, 2014, 2104, 2374, factors, 11, 17, a_exp);

                    t2cfunc::reduce(cbuffer, pbuffer, 2644, dipoles, 3, l, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<3, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // npotrec namespace

#endif /* NuclearPotentialGeom010SumErfRecFG_hpp */
