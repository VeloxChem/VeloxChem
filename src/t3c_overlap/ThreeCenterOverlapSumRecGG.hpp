#ifndef ThreeCenterOverlapSumRecGG_hpp
#define ThreeCenterOverlapSumRecGG_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "OverlapPrimRecSS.hpp"
#include "ThreeCenterOverlapPrimRecSS.hpp"
#include "ThreeCenterOverlapPrimRecSP.hpp"
#include "ThreeCenterOverlapPrimRecSD.hpp"
#include "ThreeCenterOverlapPrimRecSF.hpp"
#include "ThreeCenterOverlapPrimRecSG.hpp"
#include "ThreeCenterOverlapPrimRecPP.hpp"
#include "ThreeCenterOverlapPrimRecPD.hpp"
#include "ThreeCenterOverlapPrimRecPF.hpp"
#include "ThreeCenterOverlapPrimRecPG.hpp"
#include "ThreeCenterOverlapPrimRecDD.hpp"
#include "ThreeCenterOverlapPrimRecDF.hpp"
#include "ThreeCenterOverlapPrimRecDG.hpp"
#include "ThreeCenterOverlapPrimRecFF.hpp"
#include "ThreeCenterOverlapPrimRecFG.hpp"
#include "ThreeCenterOverlapPrimRecGG.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t3ovlrec { // t3ovlrec namespace

/// @brief Computes (G|G(r)|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_sum_overlap_gg(T& distributor,
                    const CGtoBlock& bra_gto_block,
                    const CGtoBlock& ket_gto_block,
                    const std::pair<size_t, size_t>& bra_indices,
                    const std::pair<size_t, size_t>& ket_indices,
                    const bool bra_eq_ket) -> void
{
    // intialize external coordinate(s)

    const auto coords = distributor.coordinates();

    // intialize external Gaussian(s)

    const auto exgtos = distributor.data();

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

    CSimdArray<double> pbuffer(799, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(225, 1);

    CSimdArray<double> sbuffer(81, 1);

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

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                const size_t npoints = coords.size();

                for (size_t l = 0; l < npoints; l++)
                {
                    t2cfunc::comp_distances_pc(factors, 17, 8, coords[l]);

                    t2cfunc::comp_distances_ga(factors, 11, 8, r_a, coords[l], a_exp, exgtos[l]);

                    t2cfunc::comp_distances_gb(factors, 14, 8, 2, coords[l], a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_ss(pbuffer, 1, 0, factors, 17, a_exp, exgtos[l], exgtos[npoints + l]);

                    t3ovlrec::comp_prim_overlap_sp(pbuffer, 2, 1, factors, 14);

                    t3ovlrec::comp_prim_overlap_sd(pbuffer, 5, 1, 2, factors, 14, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_sf(pbuffer, 11, 2, 5, factors, 14, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_sg(pbuffer, 21, 5, 11, factors, 14, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_pp(pbuffer, 36, 1, 2, factors, 11, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_pd(pbuffer, 45, 2, 5, factors, 11, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_pf(pbuffer, 63, 5, 11, factors, 11, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_pg(pbuffer, 93, 11, 21, factors, 11, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_dd(pbuffer, 138, 5, 36, 45, factors, 11, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_df(pbuffer, 174, 11, 45, 63, factors, 11, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_dg(pbuffer, 234, 21, 63, 93, factors, 11, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_ff(pbuffer, 324, 63, 138, 174, factors, 11, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_fg(pbuffer, 424, 93, 174, 234, factors, 11, a_exp, exgtos[l]);

                    t3ovlrec::comp_prim_overlap_gg(pbuffer, 574, 234, 324, 424, factors, 11, a_exp, exgtos[l]);

                    t2cfunc::reduce(cbuffer, pbuffer, 574, ket_width, ket_npgtos);
                }
            }

            t2cfunc::transform<4, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // t3ovlrec namespace

#endif /* ThreeCenterOverlapSumRecGG_hpp */
