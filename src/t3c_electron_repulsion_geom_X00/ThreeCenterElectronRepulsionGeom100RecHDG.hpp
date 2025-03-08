#ifndef ThreeCenterElectronRepulsionGeom100RecHDG_hpp
#define ThreeCenterElectronRepulsionGeom100RecHDG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionGeom100ContrRecHXX.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSI.hpp"
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

/// @brief Computes d^(1)/dA^(1)(H|1/|r-r'||DG)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_hdg(T& distributor,
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

    CSimdArray<double> pbuffer(12916, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(8576, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(8646, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1485, 1);

    // setup Boys fuction data

    const CBoysFunc<12> bf_table;

    CSimdArray<double> bf_data(14, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 13, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 13);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 0);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 6, pfactors, 16, bf_data, 6);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 7, pfactors, 16, bf_data, 7);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 8, pfactors, 16, bf_data, 8);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 9, pfactors, 16, bf_data, 9);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 10, pfactors, 16, bf_data, 10);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 11, pfactors, 16, bf_data, 11);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 12, pfactors, 16, bf_data, 12);

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

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 714, 6, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 717, 4, 25, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 726, 5, 28, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 735, 6, 31, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 744, 25, 73, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 762, 28, 79, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 780, 31, 85, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 798, 61, 135, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 828, 67, 145, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 858, 73, 155, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 888, 79, 165, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 918, 85, 175, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 948, 135, 245, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 993, 145, 260, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1038, 155, 275, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1083, 165, 290, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1128, 175, 305, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1173, 245, 392, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1236, 260, 413, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1299, 275, 434, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1362, 290, 455, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1425, 305, 476, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1488, 392, 574, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1572, 413, 602, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1656, 434, 630, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1740, 455, 658, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1824, 476, 686, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 1908, 4, 5, 714, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 1914, 25, 28, 714, 735, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 1932, 61, 67, 717, 744, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 1968, 67, 73, 726, 762, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 2004, 73, 79, 735, 780, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 2040, 135, 145, 744, 858, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 2100, 145, 155, 762, 888, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 2160, 155, 165, 780, 918, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 2220, 215, 230, 798, 948, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 2310, 230, 245, 828, 993, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 2400, 245, 260, 858, 1038, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 2490, 260, 275, 888, 1083, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 2580, 275, 290, 918, 1128, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 2670, 350, 371, 948, 1173, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 2796, 371, 392, 993, 1236, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 2922, 392, 413, 1038, 1299, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 3048, 413, 434, 1083, 1362, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 3174, 434, 455, 1128, 1425, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 3300, 518, 546, 1173, 1488, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 3468, 546, 574, 1236, 1572, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 3636, 574, 602, 1299, 1656, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 3804, 602, 630, 1362, 1740, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 3972, 630, 658, 1425, 1824, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 4140, 717, 726, 1908, 1914, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 4170, 744, 762, 1914, 2004, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 4230, 798, 828, 1932, 2040, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 4330, 828, 858, 1968, 2100, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 4430, 858, 888, 2004, 2160, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 4530, 948, 993, 2040, 2400, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 4680, 993, 1038, 2100, 2490, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 4830, 1038, 1083, 2160, 2580, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 4980, 1173, 1236, 2400, 2922, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 5190, 1236, 1299, 2490, 3048, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 5400, 1299, 1362, 2580, 3174, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 5610, 1488, 1572, 2922, 3636, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 5890, 1572, 1656, 3048, 3804, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 6170, 1656, 1740, 3174, 3972, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 6450, 1932, 1968, 4140, 4170, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 6540, 2040, 2100, 4170, 4430, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 6690, 2220, 2310, 4230, 4530, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 6915, 2310, 2400, 4330, 4680, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 7140, 2400, 2490, 4430, 4830, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 7365, 2670, 2796, 4530, 4980, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 7680, 2796, 2922, 4680, 5190, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 7995, 2922, 3048, 4830, 5400, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 8310, 3300, 3468, 4980, 5610, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 8730, 3468, 3636, 5190, 5890, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 9150, 3636, 3804, 5400, 6170, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsf(pbuffer, 9570, 4230, 4330, 6450, 6540, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsg(pbuffer, 9780, 4530, 4680, 6540, 7140, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsh(pbuffer, 10095, 4980, 5190, 7140, 7995, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsi(pbuffer, 10536, 5610, 5890, 7995, 9150, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isg(pbuffer, 11124, 6690, 6915, 9570, 9780, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ish(pbuffer, 11544, 7365, 7680, 9780, 10095, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isi(pbuffer, 12132, 8310, 8730, 10095, 10536, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 6690, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 225, pbuffer, 7365, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 540, pbuffer, 8310, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4992, pbuffer, 11124, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5412, pbuffer, 11544, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6000, pbuffer, 12132, 784, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {11124, 11544});

                pbuffer.scale(2.0 * a_exp, {11544, 12132});

                pbuffer.scale(2.0 * a_exp, {12132, 12916});

                t2cfunc::reduce(cbuffer, 6784, pbuffer, 11124, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7204, pbuffer, 11544, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7792, pbuffer, 12132, 784, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_hxx(cbuffer, 960, 0, 6784, 0, 4);

            t3ceri::comp_bra_geom1_electron_repulsion_hxx(cbuffer, 1905, 225, 7204, 0, 5);

            t3ceri::comp_bra_geom1_electron_repulsion_hxx(cbuffer, 3228, 540, 7792, 0, 6);

            t3cfunc::bra_transform<5>(skbuffer, 0, cbuffer, 960, 0, 4);

            t3cfunc::bra_transform<5>(skbuffer, 165, cbuffer, 1275, 0, 4);

            t3cfunc::bra_transform<5>(skbuffer, 330, cbuffer, 1590, 0, 4);

            t3cfunc::bra_transform<5>(skbuffer, 495, cbuffer, 1905, 0, 5);

            t3cfunc::bra_transform<5>(skbuffer, 726, cbuffer, 2346, 0, 5);

            t3cfunc::bra_transform<5>(skbuffer, 957, cbuffer, 2787, 0, 5);

            t3cfunc::bra_transform<5>(skbuffer, 1188, cbuffer, 3228, 0, 6);

            t3cfunc::bra_transform<5>(skbuffer, 1496, cbuffer, 3816, 0, 6);

            t3cfunc::bra_transform<5>(skbuffer, 1804, cbuffer, 4404, 0, 6);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 2112, 0, 495, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 2607, 165, 726, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 3102, 330, 957, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 3597, 495, 1188, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 4290, 726, 1496, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 4983, 957, 1804, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 5676, 2112, 3597, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 6666, 2607, 4290, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 7656, 3102, 4983, cfactors, 6, 5);

            t3cfunc::ket_transform<2, 4>(sbuffer, 0, skbuffer, 5676, 5);

            t3cfunc::ket_transform<2, 4>(sbuffer, 495, skbuffer, 6666, 5);

            t3cfunc::ket_transform<2, 4>(sbuffer, 990, skbuffer, 7656, 5);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 5, 2, 4, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecHDG_hpp */
