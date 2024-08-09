#ifndef OverlapGeom100RecGG_hpp
#define OverlapGeom100RecGG_hpp

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
#include "OverlapPrimRecPS.hpp"
#include "OverlapPrimRecPP.hpp"
#include "OverlapPrimRecPD.hpp"
#include "OverlapPrimRecPF.hpp"
#include "OverlapPrimRecPG.hpp"
#include "OverlapPrimRecDP.hpp"
#include "OverlapPrimRecDD.hpp"
#include "OverlapPrimRecDF.hpp"
#include "OverlapPrimRecDG.hpp"
#include "OverlapPrimRecFD.hpp"
#include "OverlapPrimRecFF.hpp"
#include "OverlapPrimRecFG.hpp"
#include "OverlapPrimRecGF.hpp"
#include "OverlapPrimRecGG.hpp"
#include "OverlapPrimRecHG.hpp"
#include "GeometricalDerivatives1X0ForGY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace ovlrec { // ovlrec namespace

/// @brief Computes (d^(1)/dA^(1)G|G)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_overlap_geom_10_gg(T& distributor,
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

    CSimdArray<double> pbuffer(2019, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(675, 1);

    CSimdArray<double> sbuffer(243, 1);

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

                ovlrec::comp_prim_overlap_ps(pbuffer, 35, 0, factors, 8);

                ovlrec::comp_prim_overlap_pp(pbuffer, 38, 0, 1, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pd(pbuffer, 47, 1, 4, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pf(pbuffer, 65, 4, 10, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_pg(pbuffer, 95, 10, 20, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dp(pbuffer, 140, 1, 35, 38, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dd(pbuffer, 158, 4, 38, 47, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_df(pbuffer, 194, 10, 47, 65, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_dg(pbuffer, 254, 20, 65, 95, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fd(pbuffer, 344, 47, 140, 158, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_ff(pbuffer, 404, 65, 158, 194, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_fg(pbuffer, 504, 95, 194, 254, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gf(pbuffer, 654, 194, 344, 404, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_gg(pbuffer, 804, 254, 404, 504, factors, 8, a_exp);

                ovlrec::comp_prim_overlap_hg(pbuffer, 1029, 504, 654, 804, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_gx(pbuffer, 1344, 504, 1029, 1, 15, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 0, ket_width, ket_npgtos);
            }

            t2cfunc::transform<4, 4>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 4, 4, j, ket_range, bra_eq_ket);
        }
    }
}

} // ovlrec namespace

#endif /* OverlapGeom100RecGG_hpp */
