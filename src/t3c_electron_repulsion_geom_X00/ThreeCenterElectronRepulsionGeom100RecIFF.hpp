#ifndef ThreeCenterElectronRepulsionGeom100RecIFF_hpp
#define ThreeCenterElectronRepulsionGeom100RecIFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionGeom100ContrRecIXX.hpp"
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
#include "ThreeCenterElectronRepulsionPrimRecFSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecKSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecKSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecKSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecKSI.hpp"
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

/// @brief Computes d^(1)/dA^(1)(I|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_iff(T& distributor,
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

    CSimdArray<double> pbuffer(22279, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(13098, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(18018, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1911, 1);

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

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 714, 4, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 717, 5, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 720, 6, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 723, 2, 19, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 732, 3, 22, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 741, 4, 25, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 750, 5, 28, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 759, 6, 31, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 768, 19, 61, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 786, 22, 67, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 804, 25, 73, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 822, 28, 79, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 840, 31, 85, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 858, 49, 115, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 888, 55, 125, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 918, 61, 135, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 948, 67, 145, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 978, 73, 155, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1008, 79, 165, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1038, 85, 175, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1068, 115, 215, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1113, 125, 230, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1158, 135, 245, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1203, 145, 260, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1248, 155, 275, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1293, 165, 290, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1338, 175, 305, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1383, 215, 350, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1446, 230, 371, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1509, 245, 392, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1572, 260, 413, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1635, 275, 434, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1698, 290, 455, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1761, 305, 476, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1824, 350, 518, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1908, 371, 546, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1992, 392, 574, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2076, 413, 602, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2160, 434, 630, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2244, 455, 658, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2328, 476, 686, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 2412, 2, 3, 714, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 2418, 3, 4, 717, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 2424, 4, 5, 720, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 2430, 19, 22, 714, 741, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 2448, 22, 25, 717, 750, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 2466, 25, 28, 720, 759, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 2484, 49, 55, 723, 768, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 2520, 55, 61, 732, 786, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 2556, 61, 67, 741, 804, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 2592, 67, 73, 750, 822, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 2628, 73, 79, 759, 840, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 2664, 115, 125, 768, 918, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 2724, 125, 135, 786, 948, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 2784, 135, 145, 804, 978, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 2844, 145, 155, 822, 1008, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 2904, 155, 165, 840, 1038, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 2964, 215, 230, 918, 1158, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 3054, 230, 245, 948, 1203, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 3144, 245, 260, 978, 1248, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 3234, 260, 275, 1008, 1293, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 3324, 275, 290, 1038, 1338, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 3414, 350, 371, 1158, 1509, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 3540, 371, 392, 1203, 1572, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 3666, 392, 413, 1248, 1635, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 3792, 413, 434, 1293, 1698, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 3918, 434, 455, 1338, 1761, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 4044, 518, 546, 1509, 1992, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 4212, 546, 574, 1572, 2076, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 4380, 574, 602, 1635, 2160, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 4548, 602, 630, 1698, 2244, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 4716, 630, 658, 1761, 2328, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 4884, 714, 717, 2424, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 4894, 723, 732, 2412, 2430, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 4924, 732, 741, 2418, 2448, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 4954, 741, 750, 2424, 2466, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 4984, 768, 786, 2430, 2556, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 5044, 786, 804, 2448, 2592, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 5104, 804, 822, 2466, 2628, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 5164, 858, 888, 2484, 2664, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 5264, 888, 918, 2520, 2724, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 5364, 918, 948, 2556, 2784, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 5464, 948, 978, 2592, 2844, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 5564, 978, 1008, 2628, 2904, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 5664, 1068, 1113, 2664, 2964, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 5814, 1113, 1158, 2724, 3054, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 5964, 1158, 1203, 2784, 3144, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 6114, 1203, 1248, 2844, 3234, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 6264, 1248, 1293, 2904, 3324, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 6414, 1383, 1446, 2964, 3414, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 6624, 1446, 1509, 3054, 3540, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 6834, 1509, 1572, 3144, 3666, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 7044, 1572, 1635, 3234, 3792, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 7254, 1635, 1698, 3324, 3918, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 7464, 1824, 1908, 3414, 4044, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 7744, 1908, 1992, 3540, 4212, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 8024, 1992, 2076, 3666, 4380, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 8304, 2076, 2160, 3792, 4548, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 8584, 2160, 2244, 3918, 4716, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gss(pbuffer, 8864, 2412, 2418, 4884, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsp(pbuffer, 8879, 2430, 2448, 4884, 4954, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 8924, 2484, 2520, 4894, 4984, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 9014, 2520, 2556, 4924, 5044, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 9104, 2556, 2592, 4954, 5104, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 9194, 2664, 2724, 4984, 5364, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 9344, 2724, 2784, 5044, 5464, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 9494, 2784, 2844, 5104, 5564, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 9644, 2964, 3054, 5364, 5964, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 9869, 3054, 3144, 5464, 6114, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 10094, 3144, 3234, 5564, 6264, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 10319, 3414, 3540, 5964, 6834, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 10634, 3540, 3666, 6114, 7044, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 10949, 3666, 3792, 6264, 7254, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 11264, 4044, 4212, 6834, 8024, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 11684, 4212, 4380, 7044, 8304, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 12104, 4380, 4548, 7254, 8584, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsp(pbuffer, 12524, 4894, 4924, 8864, 8879, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsd(pbuffer, 12587, 4984, 5044, 8879, 9104, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsf(pbuffer, 12713, 5164, 5264, 8924, 9194, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsf(pbuffer, 12923, 5264, 5364, 9014, 9344, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsf(pbuffer, 13133, 5364, 5464, 9104, 9494, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsg(pbuffer, 13343, 5664, 5814, 9194, 9644, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsg(pbuffer, 13658, 5814, 5964, 9344, 9869, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsg(pbuffer, 13973, 5964, 6114, 9494, 10094, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsh(pbuffer, 14288, 6414, 6624, 9644, 10319, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsh(pbuffer, 14729, 6624, 6834, 9869, 10634, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsh(pbuffer, 15170, 6834, 7044, 10094, 10949, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsi(pbuffer, 15611, 7464, 7744, 10319, 11264, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsi(pbuffer, 16199, 7744, 8024, 10634, 11684, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsi(pbuffer, 16787, 8024, 8304, 10949, 12104, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isd(pbuffer, 17375, 8924, 9014, 12524, 12587, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isf(pbuffer, 17543, 9194, 9344, 12587, 13133, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isg(pbuffer, 17823, 9644, 9869, 13133, 13973, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ish(pbuffer, 18243, 10319, 10634, 13973, 15170, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isi(pbuffer, 18831, 11264, 11684, 15170, 16787, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ksf(pbuffer, 19615, 12713, 12923, 17375, 17543, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ksg(pbuffer, 19975, 13343, 13658, 17543, 17823, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ksh(pbuffer, 20515, 14288, 14729, 17823, 18243, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ksi(pbuffer, 21271, 15611, 16199, 18243, 18831, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 12713, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 210, pbuffer, 13343, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 525, pbuffer, 14288, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 966, pbuffer, 15611, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7770, pbuffer, 19615, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8130, pbuffer, 19975, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8670, pbuffer, 20515, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9426, pbuffer, 21271, 1008, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {19615, 19975});

                pbuffer.scale(2.0 * a_exp, {19975, 20515});

                pbuffer.scale(2.0 * a_exp, {20515, 21271});

                pbuffer.scale(2.0 * a_exp, {21271, 22279});

                t2cfunc::reduce(cbuffer, 10434, pbuffer, 19615, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10794, pbuffer, 19975, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11334, pbuffer, 20515, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12090, pbuffer, 21271, 1008, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_ixx(cbuffer, 1554, 0, 10434, 0, 3);

            t3ceri::comp_bra_geom1_electron_repulsion_ixx(cbuffer, 2394, 210, 10794, 0, 4);

            t3ceri::comp_bra_geom1_electron_repulsion_ixx(cbuffer, 3654, 525, 11334, 0, 5);

            t3ceri::comp_bra_geom1_electron_repulsion_ixx(cbuffer, 5418, 966, 12090, 0, 6);

            t3cfunc::bra_transform<6>(skbuffer, 0, cbuffer, 1554, 0, 3);

            t3cfunc::bra_transform<6>(skbuffer, 130, cbuffer, 1834, 0, 3);

            t3cfunc::bra_transform<6>(skbuffer, 260, cbuffer, 2114, 0, 3);

            t3cfunc::bra_transform<6>(skbuffer, 390, cbuffer, 2394, 0, 4);

            t3cfunc::bra_transform<6>(skbuffer, 585, cbuffer, 2814, 0, 4);

            t3cfunc::bra_transform<6>(skbuffer, 780, cbuffer, 3234, 0, 4);

            t3cfunc::bra_transform<6>(skbuffer, 975, cbuffer, 3654, 0, 5);

            t3cfunc::bra_transform<6>(skbuffer, 1248, cbuffer, 4242, 0, 5);

            t3cfunc::bra_transform<6>(skbuffer, 1521, cbuffer, 4830, 0, 5);

            t3cfunc::bra_transform<6>(skbuffer, 1794, cbuffer, 5418, 0, 6);

            t3cfunc::bra_transform<6>(skbuffer, 2158, cbuffer, 6202, 0, 6);

            t3cfunc::bra_transform<6>(skbuffer, 2522, cbuffer, 6986, 0, 6);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 2886, 0, 390, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 3276, 130, 585, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 3666, 260, 780, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 4056, 390, 975, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 4641, 585, 1248, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 5226, 780, 1521, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 5811, 975, 1794, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 6630, 1248, 2158, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 7449, 1521, 2522, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdf(skbuffer, 8268, 2886, 4056, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdf(skbuffer, 9048, 3276, 4641, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdf(skbuffer, 9828, 3666, 5226, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 10608, 4056, 5811, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 11778, 4641, 6630, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 12948, 5226, 7449, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xff(skbuffer, 14118, 8268, 10608, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xff(skbuffer, 15418, 9048, 11778, cfactors, 6, 6);

            t3ceri::comp_hrr_electron_repulsion_xff(skbuffer, 16718, 9828, 12948, cfactors, 6, 6);

            t3cfunc::ket_transform<3, 3>(sbuffer, 0, skbuffer, 14118, 6);

            t3cfunc::ket_transform<3, 3>(sbuffer, 637, skbuffer, 15418, 6);

            t3cfunc::ket_transform<3, 3>(sbuffer, 1274, skbuffer, 16718, 6);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 6, 3, 3, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecIFF_hpp */
