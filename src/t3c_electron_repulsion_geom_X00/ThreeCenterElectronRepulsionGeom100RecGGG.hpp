#ifndef ThreeCenterElectronRepulsionGeom100RecGGG_hpp
#define ThreeCenterElectronRepulsionGeom100RecGGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionGeom100ContrRecGXX.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T3CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"
#include "GtoBlock.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes d^(1)/dA^(1)(G|1/|r-r'||GG)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_ggg(T& distributor,
                                    const CGtoBlock& bra_gto_block,
                                    const CGtoPairBlock& ket_gto_pair_block,
                                    const std::pair<size_t, size_t>& bra_range) -> void
{
    // intialize GTOs data on bra side

    const auto bra_gto_coords = bra_gto_block.coordinates();

    const auto bra_gto_exps = bra_gto_block.exponents();

    const auto bra_gto_norms = bra_gto_block.normalization_factors();

    const auto bra_gto_indices = bra_gto_block.orbital_indices();

    const auto bra_ncgtos = bra_gto_block.number_of_basis_functions();

    const auto bra_npgtos = bra_gto_block.number_of_primitives();

    // intialize GTOs data on ket side

    const auto c_coords = ket_gto_pair_block.bra_coordinates();

    const auto d_coords = ket_gto_pair_block.ket_coordinates();

    const auto c_vec_exps = ket_gto_pair_block.bra_exponents();

    const auto d_vec_exps = ket_gto_pair_block.ket_exponents();

    const auto cd_vec_norms = ket_gto_pair_block.normalization_factors();

    const auto cd_vec_ovls = ket_gto_pair_block.overlap_factors();

    const auto c_indices = ket_gto_pair_block.bra_orbital_indices();

    const auto d_indices = ket_gto_pair_block.ket_orbital_indices();

    const auto ket_npgtos = ket_gto_pair_block.number_of_primitive_pairs();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> pfactors(29, ket_npgtos);

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(16210, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(11020, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(38178, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2187, 1);

    // setup Boys fuction data

    const CBoysFunc<13> bf_table;

    CSimdArray<double> bf_data(15, ket_npgtos);

    // set up ket partitioning

    const auto ket_dim = ket_gto_pair_block.number_of_contracted_pairs();

    const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());

    for (size_t i = 0; i < ket_blocks; i++)
    {
        auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), size_t{0});

        pfactors.load(c_vec_exps, ket_range, 0, ket_npgtos);

        pfactors.load(d_vec_exps, ket_range, 1, ket_npgtos);

        pfactors.load(cd_vec_ovls, ket_range, 2, ket_npgtos);

        pfactors.load(cd_vec_norms, ket_range, 3, ket_npgtos);

        pfactors.replicate_points(c_coords, ket_range, 4, ket_npgtos);

        pfactors.replicate_points(d_coords, ket_range, 7, ket_npgtos);

        cfactors.replicate_points(c_coords, ket_range, 0, 1);

        cfactors.replicate_points(d_coords, ket_range, 3, 1);

        t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // loop over basis function pairs on bra side

        for (auto j = bra_range.first; j < bra_range.second; j++)
        {
            // zero integral buffers

            cbuffer.zero();

            skbuffer.zero();

            sbuffer.zero();

            // set up coordinates on bra side

            const auto r_a = bra_gto_coords[j];

            for (int k = 0; k < bra_npgtos; k++)
            {
                const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];

                const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];

                t4cfunc::comp_coordinates_q(pfactors, 10, 4, 7);

                t3cfunc::comp_distances_aq(pfactors, 13, 10, r_a);

                t3cfunc::comp_coordinates_w(pfactors, 17, 10, r_a, a_exp);

                t4cfunc::comp_distances_qd(pfactors, 20, 10, 7);

                t4cfunc::comp_distances_wq(pfactors, 23, 17, 10);

                t4cfunc::comp_distances_wp(pfactors, 26, 17, r_a);

                t3cfunc::comp_boys_args(bf_data, 14, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 14);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 6);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 6, pfactors, 16, bf_data, 7);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 7, pfactors, 16, bf_data, 8);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 8, pfactors, 16, bf_data, 9);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 9, pfactors, 16, bf_data, 10);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 10, pfactors, 16, bf_data, 11);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 11, pfactors, 16, bf_data, 12);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 12, pfactors, 16, bf_data, 13);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 13, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 16, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 19, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 22, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 25, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 28, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 31, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 34, 7, 8, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 37, 8, 9, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 40, 9, 10, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 43, 10, 11, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 46, 11, 12, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 49, 0, 1, 13, 16, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 55, 1, 2, 16, 19, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 61, 2, 3, 19, 22, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 67, 3, 4, 22, 25, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 73, 4, 5, 25, 28, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 79, 5, 6, 28, 31, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 85, 6, 7, 31, 34, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 91, 7, 8, 34, 37, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 97, 8, 9, 37, 40, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 103, 9, 10, 40, 43, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 109, 10, 11, 43, 46, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 115, 13, 16, 49, 55, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 125, 16, 19, 55, 61, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 135, 19, 22, 61, 67, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 145, 22, 25, 67, 73, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 155, 25, 28, 73, 79, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 165, 28, 31, 79, 85, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 175, 31, 34, 85, 91, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 185, 34, 37, 91, 97, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 195, 37, 40, 97, 103, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 205, 40, 43, 103, 109, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 215, 49, 55, 115, 125, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 230, 55, 61, 125, 135, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 245, 61, 67, 135, 145, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 260, 67, 73, 145, 155, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 275, 73, 79, 155, 165, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 290, 79, 85, 165, 175, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 305, 85, 91, 175, 185, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 320, 91, 97, 185, 195, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 335, 97, 103, 195, 205, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 350, 115, 125, 215, 230, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 371, 125, 135, 230, 245, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 392, 135, 145, 245, 260, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 413, 145, 155, 260, 275, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 434, 155, 165, 275, 290, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 455, 165, 175, 290, 305, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 476, 175, 185, 305, 320, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 497, 185, 195, 320, 335, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 518, 215, 230, 350, 371, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 546, 230, 245, 371, 392, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 574, 245, 260, 392, 413, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 602, 260, 275, 413, 434, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 630, 275, 290, 434, 455, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 658, 290, 305, 455, 476, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 686, 305, 320, 476, 497, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 714, 350, 371, 518, 546, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 750, 371, 392, 546, 574, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 786, 392, 413, 574, 602, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 822, 413, 434, 602, 630, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 858, 434, 455, 630, 658, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 894, 455, 476, 658, 686, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 930, 518, 546, 714, 750, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 975, 546, 574, 750, 786, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1020, 574, 602, 786, 822, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1065, 602, 630, 822, 858, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1110, 630, 658, 858, 894, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 1155, 4, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 1158, 4, 25, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1167, 19, 61, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1185, 22, 67, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1203, 25, 73, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1221, 61, 135, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1251, 67, 145, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1281, 73, 155, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1311, 115, 215, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1356, 125, 230, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1401, 135, 245, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1446, 145, 260, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1491, 155, 275, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1536, 215, 350, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1599, 230, 371, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1662, 245, 392, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1725, 260, 413, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1788, 275, 434, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1851, 350, 518, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1935, 371, 546, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2019, 392, 574, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2103, 413, 602, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2187, 434, 630, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 2271, 518, 714, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 2379, 546, 750, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 2487, 574, 786, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 2595, 602, 822, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 2703, 630, 858, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 2811, 714, 930, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 2946, 750, 975, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 3081, 786, 1020, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 3216, 822, 1065, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 3351, 858, 1110, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 3486, 19, 22, 1155, 1158, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 3504, 61, 67, 1158, 1203, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 3540, 115, 125, 1167, 1221, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 3600, 125, 135, 1185, 1251, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 3660, 135, 145, 1203, 1281, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 3720, 215, 230, 1221, 1401, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 3810, 230, 245, 1251, 1446, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 3900, 245, 260, 1281, 1491, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 3990, 350, 371, 1401, 1662, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 4116, 371, 392, 1446, 1725, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 4242, 392, 413, 1491, 1788, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 4368, 518, 546, 1662, 2019, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 4536, 546, 574, 1725, 2103, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 4704, 574, 602, 1788, 2187, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 4872, 714, 750, 2019, 2487, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 5088, 750, 786, 2103, 2595, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 5304, 786, 822, 2187, 2703, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 5520, 930, 975, 2487, 3081, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 5790, 975, 1020, 2595, 3216, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 6060, 1020, 1065, 2703, 3351, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 6330, 1167, 1185, 3486, 3504, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 6390, 1221, 1251, 3504, 3660, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 6490, 1311, 1356, 3540, 3720, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 6640, 1356, 1401, 3600, 3810, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 6790, 1401, 1446, 3660, 3900, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 6940, 1536, 1599, 3720, 3990, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 7150, 1599, 1662, 3810, 4116, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 7360, 1662, 1725, 3900, 4242, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 7570, 1851, 1935, 3990, 4368, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 7850, 1935, 2019, 4116, 4536, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 8130, 2019, 2103, 4242, 4704, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 8410, 2271, 2379, 4368, 4872, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 8770, 2379, 2487, 4536, 5088, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 9130, 2487, 2595, 4704, 5304, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 9490, 2811, 2946, 4872, 5520, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 9940, 2946, 3081, 5088, 5790, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 10390, 3081, 3216, 5304, 6060, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 10840, 3540, 3600, 6330, 6390, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 10990, 3720, 3810, 6390, 6790, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 11215, 3990, 4116, 6790, 7360, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 11530, 4368, 4536, 7360, 8130, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsk(pbuffer, 11950, 4872, 5088, 8130, 9130, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsl(pbuffer, 12490, 5520, 5790, 9130, 10390, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsg(pbuffer, 13165, 6490, 6640, 10840, 10990, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsh(pbuffer, 13480, 6940, 7150, 10990, 11215, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsi(pbuffer, 13921, 7570, 7850, 11215, 11530, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsk(pbuffer, 14509, 8410, 8770, 11530, 11950, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsl(pbuffer, 15265, 9490, 9940, 11950, 12490, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 6490, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 150, pbuffer, 6940, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 360, pbuffer, 7570, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 640, pbuffer, 8410, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1000, pbuffer, 9490, 450, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {13165, 13480});

                pbuffer.scale(2.0 * a_exp, {13480, 13921});

                pbuffer.scale(2.0 * a_exp, {13921, 14509});

                pbuffer.scale(2.0 * a_exp, {14509, 15265});

                pbuffer.scale(2.0 * a_exp, {15265, 16210});

                t2cfunc::reduce(cbuffer, 7975, pbuffer, 13165, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8290, pbuffer, 13480, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8731, pbuffer, 13921, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9319, pbuffer, 14509, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10075, pbuffer, 15265, 945, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_gxx(cbuffer, 1450, 0, 7975, 0, 4);

            t3ceri::comp_bra_geom1_electron_repulsion_gxx(cbuffer, 2125, 150, 8290, 0, 5);

            t3ceri::comp_bra_geom1_electron_repulsion_gxx(cbuffer, 3070, 360, 8731, 0, 6);

            t3ceri::comp_bra_geom1_electron_repulsion_gxx(cbuffer, 4330, 640, 9319, 0, 7);

            t3ceri::comp_bra_geom1_electron_repulsion_gxx(cbuffer, 5950, 1000, 10075, 0, 8);

            t3cfunc::bra_transform<4>(skbuffer, 0, cbuffer, 1450, 0, 4);

            t3cfunc::bra_transform<4>(skbuffer, 135, cbuffer, 1675, 0, 4);

            t3cfunc::bra_transform<4>(skbuffer, 270, cbuffer, 1900, 0, 4);

            t3cfunc::bra_transform<4>(skbuffer, 405, cbuffer, 2125, 0, 5);

            t3cfunc::bra_transform<4>(skbuffer, 594, cbuffer, 2440, 0, 5);

            t3cfunc::bra_transform<4>(skbuffer, 783, cbuffer, 2755, 0, 5);

            t3cfunc::bra_transform<4>(skbuffer, 972, cbuffer, 3070, 0, 6);

            t3cfunc::bra_transform<4>(skbuffer, 1224, cbuffer, 3490, 0, 6);

            t3cfunc::bra_transform<4>(skbuffer, 1476, cbuffer, 3910, 0, 6);

            t3cfunc::bra_transform<4>(skbuffer, 1728, cbuffer, 4330, 0, 7);

            t3cfunc::bra_transform<4>(skbuffer, 2052, cbuffer, 4870, 0, 7);

            t3cfunc::bra_transform<4>(skbuffer, 2376, cbuffer, 5410, 0, 7);

            t3cfunc::bra_transform<4>(skbuffer, 2700, cbuffer, 5950, 0, 8);

            t3cfunc::bra_transform<4>(skbuffer, 3105, cbuffer, 6625, 0, 8);

            t3cfunc::bra_transform<4>(skbuffer, 3510, cbuffer, 7300, 0, 8);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 3915, 0, 405, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 4320, 135, 594, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 4725, 270, 783, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 5130, 405, 972, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 5697, 594, 1224, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 6264, 783, 1476, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 6831, 972, 1728, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 7587, 1224, 2052, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 8343, 1476, 2376, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 9099, 1728, 2700, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 10071, 2052, 3105, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 11043, 2376, 3510, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 12015, 3915, 5130, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 12825, 4320, 5697, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 13635, 4725, 6264, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 14445, 5130, 6831, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 15579, 5697, 7587, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 16713, 6264, 8343, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 17847, 6831, 9099, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 19359, 7587, 10071, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 20871, 8343, 11043, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 22383, 12015, 14445, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 23733, 12825, 15579, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 25083, 13635, 16713, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 26433, 14445, 17847, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 28323, 15579, 19359, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 30213, 16713, 20871, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 32103, 22383, 26433, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 34128, 23733, 28323, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 36153, 25083, 30213, cfactors, 6, 4);

            t3cfunc::ket_transform<4, 4>(sbuffer, 0, skbuffer, 32103, 4);

            t3cfunc::ket_transform<4, 4>(sbuffer, 729, skbuffer, 34128, 4);

            t3cfunc::ket_transform<4, 4>(sbuffer, 1458, skbuffer, 36153, 4);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 4, 4, 4, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecGGG_hpp */
