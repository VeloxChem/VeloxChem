#ifndef TwoCenterElectronRepulsionGeom100RecHF_hpp
#define TwoCenterElectronRepulsionGeom100RecHF_hpp

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
#include "TwoCenterElectronRepulsionPrimRecDS.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecFS.hpp"
#include "TwoCenterElectronRepulsionPrimRecFP.hpp"
#include "TwoCenterElectronRepulsionPrimRecFD.hpp"
#include "TwoCenterElectronRepulsionPrimRecFF.hpp"
#include "TwoCenterElectronRepulsionPrimRecGP.hpp"
#include "TwoCenterElectronRepulsionPrimRecGD.hpp"
#include "TwoCenterElectronRepulsionPrimRecGF.hpp"
#include "TwoCenterElectronRepulsionPrimRecHD.hpp"
#include "TwoCenterElectronRepulsionPrimRecHF.hpp"
#include "TwoCenterElectronRepulsionPrimRecIF.hpp"
#include "GeometricalDerivatives1X0ForHY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)H|1/|r-r'||F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_hf(T& distributor,
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

    CSimdArray<double> pbuffer(3243, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(630, 1);

    CSimdArray<double> sbuffer(231, 1);

    // setup Boys function data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_data(11, ket_npgtos);

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

                t2cfunc::comp_boys_args_with_rho(bf_data, 9, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 9);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 7, bf_data, 8, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 8, bf_data, 9, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 9, 0, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 12, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 15, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 18, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 21, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 24, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 27, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 30, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 33, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 36, 0, 1, 15, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 42, 1, 2, 18, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 48, 2, 3, 21, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 54, 3, 4, 24, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 60, 4, 5, 27, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 66, 5, 6, 30, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 72, 6, 7, 33, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 78, 9, 12, 36, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 88, 12, 15, 42, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 98, 15, 18, 48, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 108, 18, 21, 54, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 118, 21, 24, 60, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 128, 24, 27, 66, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 138, 27, 30, 72, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 148, 3, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 151, 4, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 154, 5, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 157, 3, 21, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 166, 4, 24, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 175, 5, 27, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 184, 15, 42, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 202, 18, 48, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 220, 21, 54, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 238, 24, 60, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 256, 27, 66, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 274, 42, 98, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 304, 48, 108, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 334, 54, 118, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 364, 60, 128, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 394, 66, 138, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 424, 3, 4, 154, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 430, 15, 18, 148, 157, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 448, 18, 21, 151, 166, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 466, 21, 24, 154, 175, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 484, 42, 48, 157, 220, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 520, 48, 54, 166, 238, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 556, 54, 60, 175, 256, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 592, 78, 88, 184, 274, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 652, 88, 98, 202, 304, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 712, 98, 108, 220, 334, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 772, 108, 118, 238, 364, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 832, 118, 128, 256, 394, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 892, 148, 151, 424, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 902, 157, 166, 424, 466, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 932, 184, 202, 430, 484, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 992, 202, 220, 448, 520, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 1052, 220, 238, 466, 556, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1112, 274, 304, 484, 712, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1212, 304, 334, 520, 772, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1312, 334, 364, 556, 832, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 1412, 430, 448, 892, 902, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 1457, 484, 520, 902, 1052, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 1547, 592, 652, 932, 1112, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 1697, 652, 712, 992, 1212, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 1847, 712, 772, 1052, 1312, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hd(pbuffer, 1997, 932, 992, 1412, 1457, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hf(pbuffer, 2123, 1112, 1212, 1457, 1847, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_if(pbuffer, 2333, 1547, 1697, 1997, 2123, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_hx(pbuffer, 2613, 1547, 2333, 1, 10, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 2613, ket_width, ket_npgtos);
            }

            t2cfunc::transform<5, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 5, 3, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecHF_hpp */
