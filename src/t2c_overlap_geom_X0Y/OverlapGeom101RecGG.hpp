#ifndef OverlapGeom101RecGG_hpp
#define OverlapGeom101RecGG_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "OverlapPrimRecSS.hpp"
#include "OverlapPrimRecSP.hpp"
#include "OverlapPrimRecSD.hpp"
#include "OverlapPrimRecSF.hpp"
#include "OverlapPrimRecSG.hpp"
#include "OverlapPrimRecSH.hpp"
#include "OverlapPrimRecPS.hpp"
#include "OverlapPrimRecPP.hpp"
#include "OverlapPrimRecPD.hpp"
#include "OverlapPrimRecPF.hpp"
#include "OverlapPrimRecPG.hpp"
#include "OverlapPrimRecPH.hpp"
#include "OverlapPrimRecDS.hpp"
#include "OverlapPrimRecDP.hpp"
#include "OverlapPrimRecDD.hpp"
#include "OverlapPrimRecDF.hpp"
#include "OverlapPrimRecDG.hpp"
#include "OverlapPrimRecDH.hpp"
#include "OverlapPrimRecFP.hpp"
#include "OverlapPrimRecFD.hpp"
#include "OverlapPrimRecFF.hpp"
#include "OverlapPrimRecFG.hpp"
#include "OverlapPrimRecFH.hpp"
#include "OverlapPrimRecGD.hpp"
#include "OverlapPrimRecGF.hpp"
#include "OverlapPrimRecGG.hpp"
#include "OverlapPrimRecGH.hpp"
#include "OverlapPrimRecHF.hpp"
#include "OverlapPrimRecHH.hpp"
#include "GeometricalDerivatives1X1ForGG.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace ovlrec { // ovlrec namespace

/// @brief Computes (d^(1)/dA^(1)G|d^(1)/dB^(1)G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_overlap_geom_11_gg(T& distributor,
                        const CGtoBlock& bra_gto_block,
                        const CGtoBlock& ket_gto_block,
                        const std::pair<size_t, size_t>& bra_indices,
                        const std::pair<size_t, size_t>& ket_indices,
                        const bool bra_eq_ket) -> void
{
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

    CSimdArray<double> factors(14, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(4566, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(2025, 1);

    CSimdArray<double> sbuffer(729, 1);

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

                t2cfunc::comp_distances_pa(factors, 8, 5, a_exp);

                t2cfunc::comp_distances_pb(factors, 11, 5, a_exp);

                ovlrec::comp_prim_overlap_ss(pbuffer, 0, factors, a_exp, a_norm);

                ovlrec::comp_prim_overlap_sp(pbuffer, 1, 0, factors, 11);

                ovlrec::comp_prim_overlap_sd(pbuffer, 4, 0, 1, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_sf(pbuffer, 10, 1, 4, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_sg(pbuffer, 20, 4, 10, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_sh(pbuffer, 35, 10, 20, factors, 11, a_exp);

                ovlrec::comp_prim_overlap_ps(pbuffer, 56, 0, factors, 8);

                ovlrec::comp_prim_overlap_pp(pbuffer, 59, 0, 1, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pd(pbuffer, 68, 1, 4, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pf(pbuffer, 86, 4, 10, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pg(pbuffer, 116, 10, 20, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_ph(pbuffer, 161, 20, 35, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_ds(pbuffer, 224, 0, 56, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dp(pbuffer, 230, 1, 56, 59, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dd(pbuffer, 248, 4, 59, 68, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_df(pbuffer, 284, 10, 68, 86, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dg(pbuffer, 344, 20, 86, 116, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dh(pbuffer, 434, 35, 116, 161, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fp(pbuffer, 560, 59, 224, 230, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fd(pbuffer, 590, 68, 230, 248, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_ff(pbuffer, 650, 86, 248, 284, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fg(pbuffer, 750, 116, 284, 344, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fh(pbuffer, 900, 161, 344, 434, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gd(pbuffer, 1110, 248, 560, 590, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gf(pbuffer, 1200, 284, 590, 650, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gg(pbuffer, 1350, 344, 650, 750, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gh(pbuffer, 1575, 434, 750, 900, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_hf(pbuffer, 1890, 650, 1110, 1200, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_hh(pbuffer, 2100, 900, 1350, 1575, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_11_gg(pbuffer, 2541, 650, 900, 1890, 2100, 1, factors, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 0, ket_width, ket_npgtos);
            }

            t2cfunc::transform<4, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // ovlrec namespace

#endif /* OverlapGeom101RecGG_hpp */
