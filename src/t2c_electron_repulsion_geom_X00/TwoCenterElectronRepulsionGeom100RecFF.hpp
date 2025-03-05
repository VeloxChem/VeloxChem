#ifndef TwoCenterElectronRepulsionGeom100RecFF_hpp
#define TwoCenterElectronRepulsionGeom100RecFF_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "TwoCenterElectronRepulsionPrimRecSS.hpp"
#include "TwoCenterElectronRepulsionPrimRecSP.hpp"
#include "TwoCenterElectronRepulsionPrimRecSD.hpp"
#include "TwoCenterElectronRepulsionPrimRecSF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPS.hpp"
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecFD.hpp"
#include "TwoCenterElectronRepulsionPrimRecFF.hpp"
#include "TwoCenterElectronRepulsionPrimRecGF.hpp"
#include "GeometricalDerivatives1X0ForFY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)F|1/|r-r'||F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_ff(T& distributor,
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

    CSimdArray<double> pbuffer(1108, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(300, 1);

    CSimdArray<double> sbuffer(147, 1);

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

                t2cfunc::comp_distances_pa(factors, 8, 5, a_exp);

                t2cfunc::comp_distances_pb(factors, 11, 5, a_exp);

                t2cfunc::comp_boys_args_with_rho(bf_data, 7, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 7);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 7, 0, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 10, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 13, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 16, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 19, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 22, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 25, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 28, 0, 1, 13, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 34, 1, 2, 16, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 40, 2, 3, 19, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 46, 3, 4, 22, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 52, 4, 5, 25, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 58, 7, 10, 28, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 68, 10, 13, 34, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 78, 13, 16, 40, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 88, 16, 19, 46, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 98, 19, 22, 52, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 108, 3, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 111, 3, 19, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 120, 13, 34, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 138, 16, 40, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 156, 19, 46, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 174, 34, 78, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 204, 40, 88, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 234, 46, 98, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 264, 13, 16, 108, 111, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 282, 34, 40, 111, 156, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 318, 58, 68, 120, 174, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 378, 68, 78, 138, 204, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 438, 78, 88, 156, 234, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 498, 120, 138, 264, 282, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 558, 174, 204, 282, 438, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 658, 318, 378, 498, 558, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_fx(pbuffer, 808, 318, 658, 1, 10, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 808, ket_width, ket_npgtos);
            }

            t2cfunc::transform<3, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 3, 3, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecFF_hpp */
