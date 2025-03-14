#ifndef TwoCenterElectronRepulsionGeom100RecII_hpp
#define TwoCenterElectronRepulsionGeom100RecII_hpp

#include <cstddef>
#include <array>
#include <utility>

#include "GtoBlock.hpp"
#include "SimdArray.hpp"
#include "TwoCenterElectronRepulsionPrimRecSS.hpp"
#include "TwoCenterElectronRepulsionPrimRecSP.hpp"
#include "TwoCenterElectronRepulsionPrimRecSD.hpp"
#include "TwoCenterElectronRepulsionPrimRecSF.hpp"
#include "TwoCenterElectronRepulsionPrimRecSG.hpp"
#include "TwoCenterElectronRepulsionPrimRecSH.hpp"
#include "TwoCenterElectronRepulsionPrimRecSI.hpp"
#include "TwoCenterElectronRepulsionPrimRecPS.hpp"
#include "TwoCenterElectronRepulsionPrimRecPP.hpp"
#include "TwoCenterElectronRepulsionPrimRecPD.hpp"
#include "TwoCenterElectronRepulsionPrimRecPF.hpp"
#include "TwoCenterElectronRepulsionPrimRecPG.hpp"
#include "TwoCenterElectronRepulsionPrimRecPH.hpp"
#include "TwoCenterElectronRepulsionPrimRecPI.hpp"
#include "TwoCenterElectronRepulsionPrimRecDP.hpp"
#include "TwoCenterElectronRepulsionPrimRecDD.hpp"
#include "TwoCenterElectronRepulsionPrimRecDF.hpp"
#include "TwoCenterElectronRepulsionPrimRecDG.hpp"
#include "TwoCenterElectronRepulsionPrimRecDH.hpp"
#include "TwoCenterElectronRepulsionPrimRecDI.hpp"
#include "TwoCenterElectronRepulsionPrimRecFD.hpp"
#include "TwoCenterElectronRepulsionPrimRecFF.hpp"
#include "TwoCenterElectronRepulsionPrimRecFG.hpp"
#include "TwoCenterElectronRepulsionPrimRecFH.hpp"
#include "TwoCenterElectronRepulsionPrimRecFI.hpp"
#include "TwoCenterElectronRepulsionPrimRecGF.hpp"
#include "TwoCenterElectronRepulsionPrimRecGG.hpp"
#include "TwoCenterElectronRepulsionPrimRecGH.hpp"
#include "TwoCenterElectronRepulsionPrimRecGI.hpp"
#include "TwoCenterElectronRepulsionPrimRecHG.hpp"
#include "TwoCenterElectronRepulsionPrimRecHH.hpp"
#include "TwoCenterElectronRepulsionPrimRecHI.hpp"
#include "TwoCenterElectronRepulsionPrimRecIH.hpp"
#include "TwoCenterElectronRepulsionPrimRecII.hpp"
#include "TwoCenterElectronRepulsionPrimRecKI.hpp"
#include "GeometricalDerivatives1X0ForIY.hpp"

#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BatchFunc.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes (d^(1)/dA^(1)I|1/|r-r'||I)  integrals for pair of basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
comp_electron_repulsion_geom_10_ii(T& distributor,
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

    CSimdArray<double> pbuffer(16444, ket_npgtos);

    // allocate aligned contracted integrals

    CSimdArray<double> cbuffer(2352, 1);

    CSimdArray<double> sbuffer(507, 1);

    // setup Boys function data

    const CBoysFunc<13> bf_table;

    CSimdArray<double> bf_data(15, ket_npgtos);

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

                t2cfunc::comp_boys_args_with_rho(bf_data, 13, factors, 5, a_exp);

                bf_table.compute(bf_data, 0, 13);

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

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 10, bf_data, 11, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 11, bf_data, 12, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_ss(pbuffer, 12, bf_data, 13, factors, a_exp, a_norm);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 13, 1, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 16, 2, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 19, 3, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 22, 4, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 25, 5, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 28, 6, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 31, 7, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 34, 8, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 37, 9, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 40, 10, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 43, 11, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sp(pbuffer, 46, 12, factors, 11);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 49, 0, 1, 16, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 55, 1, 2, 19, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 61, 2, 3, 22, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 67, 3, 4, 25, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 73, 4, 5, 28, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 79, 5, 6, 31, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 85, 6, 7, 34, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 91, 7, 8, 37, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 97, 8, 9, 40, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 103, 9, 10, 43, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sd(pbuffer, 109, 10, 11, 46, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 115, 13, 16, 55, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 125, 16, 19, 61, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 135, 19, 22, 67, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 145, 22, 25, 73, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 155, 25, 28, 79, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 165, 28, 31, 85, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 175, 31, 34, 91, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 185, 34, 37, 97, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 195, 37, 40, 103, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sf(pbuffer, 205, 40, 43, 109, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 215, 49, 55, 125, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 230, 55, 61, 135, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 245, 61, 67, 145, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 260, 67, 73, 155, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 275, 73, 79, 165, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 290, 79, 85, 175, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 305, 85, 91, 185, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 320, 91, 97, 195, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sg(pbuffer, 335, 97, 103, 205, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 350, 115, 125, 230, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 371, 125, 135, 245, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 392, 135, 145, 260, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 413, 145, 155, 275, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 434, 155, 165, 290, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 455, 165, 175, 305, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 476, 175, 185, 320, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_sh(pbuffer, 497, 185, 195, 335, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 518, 215, 230, 371, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 546, 230, 245, 392, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 574, 245, 260, 413, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 602, 260, 275, 434, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 630, 275, 290, 455, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 658, 290, 305, 476, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_si(pbuffer, 686, 305, 320, 497, factors, 11, a_exp);

                t2ceri::comp_prim_electron_repulsion_ps(pbuffer, 714, 6, factors, 8);

                t2ceri::comp_prim_electron_repulsion_pp(pbuffer, 717, 6, 31, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 726, 25, 73, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 744, 28, 79, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pd(pbuffer, 762, 31, 85, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 780, 73, 155, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 810, 79, 165, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pf(pbuffer, 840, 85, 175, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 870, 135, 245, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 915, 145, 260, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 960, 155, 275, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 1005, 165, 290, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pg(pbuffer, 1050, 175, 305, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1095, 245, 392, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1158, 260, 413, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1221, 275, 434, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1284, 290, 455, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ph(pbuffer, 1347, 305, 476, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1410, 350, 518, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1494, 371, 546, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1578, 392, 574, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1662, 413, 602, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1746, 434, 630, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1830, 455, 658, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_pi(pbuffer, 1914, 476, 686, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dp(pbuffer, 1998, 25, 28, 714, 717, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dd(pbuffer, 2016, 73, 79, 717, 762, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 2052, 135, 145, 726, 780, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 2112, 145, 155, 744, 810, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_df(pbuffer, 2172, 155, 165, 762, 840, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 2232, 245, 260, 780, 960, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 2322, 260, 275, 810, 1005, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dg(pbuffer, 2412, 275, 290, 840, 1050, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2502, 350, 371, 870, 1095, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2628, 371, 392, 915, 1158, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2754, 392, 413, 960, 1221, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 2880, 413, 434, 1005, 1284, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_dh(pbuffer, 3006, 434, 455, 1050, 1347, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 3132, 518, 546, 1095, 1578, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 3300, 546, 574, 1158, 1662, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 3468, 574, 602, 1221, 1746, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 3636, 602, 630, 1284, 1830, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_di(pbuffer, 3804, 630, 658, 1347, 1914, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fd(pbuffer, 3972, 726, 744, 1998, 2016, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ff(pbuffer, 4032, 780, 810, 2016, 2172, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 4132, 870, 915, 2052, 2232, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 4282, 915, 960, 2112, 2322, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fg(pbuffer, 4432, 960, 1005, 2172, 2412, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 4582, 1095, 1158, 2232, 2754, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 4792, 1158, 1221, 2322, 2880, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fh(pbuffer, 5002, 1221, 1284, 2412, 3006, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 5212, 1410, 1494, 2502, 3132, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 5492, 1494, 1578, 2628, 3300, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 5772, 1578, 1662, 2754, 3468, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 6052, 1662, 1746, 2880, 3636, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_fi(pbuffer, 6332, 1746, 1830, 3006, 3804, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gf(pbuffer, 6612, 2052, 2112, 3972, 4032, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gg(pbuffer, 6762, 2232, 2322, 4032, 4432, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 6987, 2502, 2628, 4132, 4582, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 7302, 2628, 2754, 4282, 4792, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gh(pbuffer, 7617, 2754, 2880, 4432, 5002, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gi(pbuffer, 7932, 3132, 3300, 4582, 5772, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gi(pbuffer, 8352, 3300, 3468, 4792, 6052, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_gi(pbuffer, 8772, 3468, 3636, 5002, 6332, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hg(pbuffer, 9192, 4132, 4282, 6612, 6762, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hh(pbuffer, 9507, 4582, 4792, 6762, 7617, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hi(pbuffer, 9948, 5212, 5492, 6987, 7932, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hi(pbuffer, 10536, 5492, 5772, 7302, 8352, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_hi(pbuffer, 11124, 5772, 6052, 7617, 8772, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ih(pbuffer, 11712, 6987, 7302, 9192, 9507, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ii(pbuffer, 12300, 7932, 8352, 9507, 11124, factors, 8, a_exp);

                t2ceri::comp_prim_electron_repulsion_ki(pbuffer, 13084, 9948, 10536, 11712, 12300, factors, 8, a_exp);

                t2cgeom::comp_prim_op_geom_10_ix(pbuffer, 14092, 9948, 13084, 1, 28, a_exp);

                t2cfunc::reduce(cbuffer, pbuffer, 14092, ket_width, ket_npgtos);
            }

            t2cfunc::transform<6, 6>(sbuffer, cbuffer);

            distributor.distribute(sbuffer, bra_gto_indices, ket_gto_indices, 6, 6, j, ket_range, bra_eq_ket);
        }
    }
}

} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionGeom100RecII_hpp */
