#ifndef TwoCenterElectronRepulsionGeom100RecIF_hpp
#define TwoCenterElectronRepulsionGeom100RecIF_hpp

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
#include "TwoCenterElectronRepulsionPrimRecGS.hpp"
#include "TwoCenterElectronRepulsionPrimRecGP.hpp"
#include "TwoCenterElectronRepulsionPrimRecGD.hpp"
#include "TwoCenterElectronRepulsionPrimRecGF.hpp"
#include "TwoCenterElectronRepulsionPrimRecHP.hpp"
#include "TwoCenterElectronRepulsionPrimRecHD.hpp"
#include "TwoCenterElectronRepulsionPrimRecHF.hpp"
#include "TwoCenterElectronRepulsionPrimRecID.hpp"
#include "TwoCenterElectronRepulsionPrimRecIF.hpp"
#include "TwoCenterElectronRepulsionPrimRecKF.hpp"
#include "GeometricalDerivatives1X0ForIY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)I|1/|r-r'||F)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_if(T& distributor,
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

    CSimdArray<double> pbuffer(5088, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(840, 1);

    CSimdArray<double> sbuffer(273, 1);

    // setup Boys function data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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

                t2cfunc::comp_boys_args_with_rho(bf_data, 10, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 10);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 0, bf_data, 1, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 1, bf_data, 2, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 2, bf_data, 3, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 3, bf_data, 4, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 4, bf_data, 5, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 5, bf_data, 6, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 6, bf_data, 7, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 7, bf_data, 8, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 8, bf_data, 9, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 9, bf_data, 10, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 10, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 13, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 16, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 19, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 22, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 25, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 28, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 31, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 34, 9, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 37, 0, 1, 13, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 43, 1, 2, 16, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 49, 2, 3, 19, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 55, 3, 4, 22, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 61, 4, 5, 25, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 67, 5, 6, 28, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 73, 6, 7, 31, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 79, 7, 8, 34, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 85, 10, 13, 43, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 95, 13, 16, 49, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 105, 16, 19, 55, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 115, 19, 22, 61, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 125, 22, 25, 67, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 135, 25, 28, 73, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 145, 28, 31, 79, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 155, 4, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 158, 5, factors, 8);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 161, 6, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 164, 2, 16, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 173, 3, 19, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 182, 4, 22, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 191, 5, 25, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 200, 6, 28, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 209, 16, 49, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 227, 19, 55, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 245, 22, 61, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 263, 25, 67, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 281, 28, 73, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 299, 37, 85, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 329, 43, 95, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 359, 49, 105, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 389, 55, 115, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 419, 61, 125, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 449, 67, 135, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 479, 73, 145, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 509, 2, 3, 155, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 515, 3, 4, 158, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ds(pbuffer, 521, 4, 5, 161, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 527, 16, 19, 155, 182, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 545, 19, 22, 158, 191, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 563, 22, 25, 161, 200, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 581, 37, 43, 164, 209, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 617, 43, 49, 173, 227, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 653, 49, 55, 182, 245, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 689, 55, 61, 191, 263, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 725, 61, 67, 200, 281, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 761, 85, 95, 209, 359, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 821, 95, 105, 227, 389, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 881, 105, 115, 245, 419, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 941, 115, 125, 263, 449, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 1001, 125, 135, 281, 479, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fs(pbuffer, 1061, 155, 158, 521, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 1071, 164, 173, 509, 527, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 1101, 173, 182, 515, 545, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fp(pbuffer, 1131, 182, 191, 521, 563, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 1161, 209, 227, 527, 653, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 1221, 227, 245, 545, 689, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 1281, 245, 263, 563, 725, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1341, 299, 329, 581, 761, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1441, 329, 359, 617, 821, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1541, 359, 389, 653, 881, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1641, 389, 419, 689, 941, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 1741, 419, 449, 725, 1001, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gs(pbuffer, 1841, 509, 515, 1061, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gp(pbuffer, 1856, 527, 545, 1061, 1131, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 1901, 581, 617, 1071, 1161, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 1991, 617, 653, 1101, 1221, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gd(pbuffer, 2081, 653, 689, 1131, 1281, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 2171, 761, 821, 1161, 1541, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 2321, 821, 881, 1221, 1641, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 2471, 881, 941, 1281, 1741, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hp(pbuffer, 2621, 1071, 1101, 1841, 1856, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hd(pbuffer, 2684, 1161, 1221, 1856, 2081, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hf(pbuffer, 2810, 1341, 1441, 1901, 2171, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hf(pbuffer, 3020, 1441, 1541, 1991, 2321, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hf(pbuffer, 3230, 1541, 1641, 2081, 2471, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_id(pbuffer, 3440, 1901, 1991, 2621, 2684, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_if(pbuffer, 3608, 2171, 2321, 2684, 3230, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_kf(pbuffer, 3888, 2810, 3020, 3440, 3608, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_ix(pbuffer, 4248, 2810, 3888, 1, 10, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 4248, ket_width, ket_npgtos);
            }

            t2cfunc::transform<6, 3>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 6, 3, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecIF_hpp */
