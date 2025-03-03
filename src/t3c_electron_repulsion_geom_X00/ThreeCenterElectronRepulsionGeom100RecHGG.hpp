#ifndef ThreeCenterElectronRepulsionGeom100RecHGG_hpp
#define ThreeCenterElectronRepulsionGeom100RecHGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDI.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXGG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPI.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPK.hpp"
#include "ThreeCenterElectronRepulsionGeom100ContrRecHXX.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISL.hpp"
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

/// @brief Computes d^(1)/dA^(1)(H|1/|r-r'||GG)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_hgg(T& distributor,
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

    CSimdArray<double> pbuffer(27376, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(19430, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(46662, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2673, 1);

    // setup Boys fuction data

    const CBoysFunc<14> bf_table;

    CSimdArray<double> bf_data(16, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 15, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 15);

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

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 13, pfactors, 16, bf_data, 13);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 14, pfactors, 16, bf_data, 14);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 15, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 18, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 21, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 24, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 27, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 30, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 33, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 36, 7, 8, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 39, 8, 9, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 42, 9, 10, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 45, 10, 11, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 48, 11, 12, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 51, 12, 13, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 54, 13, 14, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 57, 0, 1, 15, 18, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 63, 1, 2, 18, 21, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 69, 2, 3, 21, 24, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 75, 3, 4, 24, 27, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 81, 4, 5, 27, 30, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 87, 5, 6, 30, 33, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 93, 6, 7, 33, 36, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 99, 7, 8, 36, 39, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 105, 8, 9, 39, 42, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 111, 9, 10, 42, 45, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 117, 10, 11, 45, 48, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 123, 11, 12, 48, 51, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 129, 12, 13, 51, 54, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 135, 15, 18, 57, 63, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 145, 18, 21, 63, 69, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 155, 21, 24, 69, 75, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 165, 24, 27, 75, 81, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 175, 27, 30, 81, 87, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 185, 30, 33, 87, 93, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 195, 33, 36, 93, 99, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 205, 36, 39, 99, 105, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 215, 39, 42, 105, 111, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 225, 42, 45, 111, 117, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 235, 45, 48, 117, 123, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 245, 48, 51, 123, 129, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 255, 57, 63, 135, 145, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 270, 63, 69, 145, 155, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 285, 69, 75, 155, 165, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 300, 75, 81, 165, 175, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 315, 81, 87, 175, 185, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 330, 87, 93, 185, 195, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 345, 93, 99, 195, 205, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 360, 99, 105, 205, 215, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 375, 105, 111, 215, 225, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 390, 111, 117, 225, 235, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 405, 117, 123, 235, 245, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 420, 135, 145, 255, 270, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 441, 145, 155, 270, 285, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 462, 155, 165, 285, 300, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 483, 165, 175, 300, 315, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 504, 175, 185, 315, 330, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 525, 185, 195, 330, 345, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 546, 195, 205, 345, 360, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 567, 205, 215, 360, 375, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 588, 215, 225, 375, 390, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 609, 225, 235, 390, 405, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 630, 255, 270, 420, 441, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 658, 270, 285, 441, 462, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 686, 285, 300, 462, 483, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 714, 300, 315, 483, 504, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 742, 315, 330, 504, 525, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 770, 330, 345, 525, 546, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 798, 345, 360, 546, 567, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 826, 360, 375, 567, 588, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 854, 375, 390, 588, 609, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 882, 420, 441, 630, 658, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 918, 441, 462, 658, 686, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 954, 462, 483, 686, 714, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 990, 483, 504, 714, 742, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 1026, 504, 525, 742, 770, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 1062, 525, 546, 770, 798, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 1098, 546, 567, 798, 826, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 1134, 567, 588, 826, 854, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1170, 630, 658, 882, 918, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1215, 658, 686, 918, 954, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1260, 686, 714, 954, 990, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1305, 714, 742, 990, 1026, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1350, 742, 770, 1026, 1062, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1395, 770, 798, 1062, 1098, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 1440, 798, 826, 1098, 1134, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 1485, 6, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 1488, 4, 27, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 1497, 5, 30, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 1506, 6, 33, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1515, 27, 81, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1533, 30, 87, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1551, 33, 93, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1569, 69, 155, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1599, 75, 165, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1629, 81, 175, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1659, 87, 185, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1689, 93, 195, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1719, 155, 285, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1764, 165, 300, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1809, 175, 315, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1854, 185, 330, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1899, 195, 345, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1944, 285, 462, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 2007, 300, 483, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 2070, 315, 504, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 2133, 330, 525, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 2196, 345, 546, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2259, 462, 686, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2343, 483, 714, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2427, 504, 742, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2511, 525, 770, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 2595, 546, 798, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 2679, 686, 954, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 2787, 714, 990, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 2895, 742, 1026, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 3003, 770, 1062, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 3111, 798, 1098, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 3219, 954, 1260, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 3354, 990, 1305, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 3489, 1026, 1350, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 3624, 1062, 1395, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 3759, 1098, 1440, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 3894, 4, 5, 1485, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 3900, 27, 30, 1485, 1506, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 3918, 69, 75, 1488, 1515, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 3954, 75, 81, 1497, 1533, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 3990, 81, 87, 1506, 1551, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 4026, 155, 165, 1515, 1629, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 4086, 165, 175, 1533, 1659, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 4146, 175, 185, 1551, 1689, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 4206, 255, 270, 1569, 1719, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 4296, 270, 285, 1599, 1764, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 4386, 285, 300, 1629, 1809, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 4476, 300, 315, 1659, 1854, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 4566, 315, 330, 1689, 1899, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 4656, 420, 441, 1719, 1944, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 4782, 441, 462, 1764, 2007, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 4908, 462, 483, 1809, 2070, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 5034, 483, 504, 1854, 2133, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 5160, 504, 525, 1899, 2196, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 5286, 630, 658, 1944, 2259, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 5454, 658, 686, 2007, 2343, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 5622, 686, 714, 2070, 2427, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 5790, 714, 742, 2133, 2511, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 5958, 742, 770, 2196, 2595, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 6126, 882, 918, 2259, 2679, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 6342, 918, 954, 2343, 2787, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 6558, 954, 990, 2427, 2895, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 6774, 990, 1026, 2511, 3003, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 6990, 1026, 1062, 2595, 3111, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 7206, 1170, 1215, 2679, 3219, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 7476, 1215, 1260, 2787, 3354, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 7746, 1260, 1305, 2895, 3489, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 8016, 1305, 1350, 3003, 3624, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 8286, 1350, 1395, 3111, 3759, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 8556, 1488, 1497, 3894, 3900, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 8586, 1515, 1533, 3900, 3990, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 8646, 1569, 1599, 3918, 4026, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 8746, 1599, 1629, 3954, 4086, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 8846, 1629, 1659, 3990, 4146, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 8946, 1719, 1764, 4026, 4386, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 9096, 1764, 1809, 4086, 4476, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 9246, 1809, 1854, 4146, 4566, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 9396, 1944, 2007, 4386, 4908, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 9606, 2007, 2070, 4476, 5034, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 9816, 2070, 2133, 4566, 5160, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 10026, 2259, 2343, 4908, 5622, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 10306, 2343, 2427, 5034, 5790, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 10586, 2427, 2511, 5160, 5958, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 10866, 2679, 2787, 5622, 6558, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 11226, 2787, 2895, 5790, 6774, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 11586, 2895, 3003, 5958, 6990, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 11946, 3219, 3354, 6558, 7746, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 12396, 3354, 3489, 6774, 8016, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 12846, 3489, 3624, 6990, 8286, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 13296, 3918, 3954, 8556, 8586, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 13386, 4026, 4086, 8586, 8846, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 13536, 4206, 4296, 8646, 8946, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 13761, 4296, 4386, 8746, 9096, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 13986, 4386, 4476, 8846, 9246, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 14211, 4656, 4782, 8946, 9396, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 14526, 4782, 4908, 9096, 9606, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 14841, 4908, 5034, 9246, 9816, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 15156, 5286, 5454, 9396, 10026, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 15576, 5454, 5622, 9606, 10306, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 15996, 5622, 5790, 9816, 10586, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsk(pbuffer, 16416, 6126, 6342, 10026, 10866, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsk(pbuffer, 16956, 6342, 6558, 10306, 11226, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsk(pbuffer, 17496, 6558, 6774, 10586, 11586, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsl(pbuffer, 18036, 7206, 7476, 10866, 11946, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsl(pbuffer, 18711, 7476, 7746, 11226, 12396, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsl(pbuffer, 19386, 7746, 8016, 11586, 12846, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsf(pbuffer, 20061, 8646, 8746, 13296, 13386, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsg(pbuffer, 20271, 8946, 9096, 13386, 13986, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsh(pbuffer, 20586, 9396, 9606, 13986, 14841, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsi(pbuffer, 21027, 10026, 10306, 14841, 15996, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsk(pbuffer, 21615, 10866, 11226, 15996, 17496, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsl(pbuffer, 22371, 11946, 12396, 17496, 19386, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isg(pbuffer, 23316, 13536, 13761, 20061, 20271, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_ish(pbuffer, 23736, 14211, 14526, 20271, 20586, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isi(pbuffer, 24324, 15156, 15576, 20586, 21027, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isk(pbuffer, 25108, 16416, 16956, 21027, 21615, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isl(pbuffer, 26116, 18036, 18711, 21615, 22371, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 13536, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 225, pbuffer, 14211, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 540, pbuffer, 15156, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 960, pbuffer, 16416, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1500, pbuffer, 18036, 675, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11310, pbuffer, 23316, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11730, pbuffer, 23736, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12318, pbuffer, 24324, 784, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13102, pbuffer, 25108, 1008, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14110, pbuffer, 26116, 1260, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {23316, 23736});

                pbuffer.scale(2.0 * a_exp, {23736, 24324});

                pbuffer.scale(2.0 * a_exp, {24324, 25108});

                pbuffer.scale(2.0 * a_exp, {25108, 26116});

                pbuffer.scale(2.0 * a_exp, {26116, 27376});

                t2cfunc::reduce(cbuffer, 15370, pbuffer, 23316, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15790, pbuffer, 23736, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16378, pbuffer, 24324, 784, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 17162, pbuffer, 25108, 1008, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18170, pbuffer, 26116, 1260, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_hxx(cbuffer, 2175, 0, 15370, 0, 4);

            t3ceri::comp_bra_geom1_electron_repulsion_hxx(cbuffer, 3120, 225, 15790, 0, 5);

            t3ceri::comp_bra_geom1_electron_repulsion_hxx(cbuffer, 4443, 540, 16378, 0, 6);

            t3ceri::comp_bra_geom1_electron_repulsion_hxx(cbuffer, 6207, 960, 17162, 0, 7);

            t3ceri::comp_bra_geom1_electron_repulsion_hxx(cbuffer, 8475, 1500, 18170, 0, 8);

            t3cfunc::bra_transform<5>(skbuffer, 0, cbuffer, 2175, 0, 4);

            t3cfunc::bra_transform<5>(skbuffer, 165, cbuffer, 2490, 0, 4);

            t3cfunc::bra_transform<5>(skbuffer, 330, cbuffer, 2805, 0, 4);

            t3cfunc::bra_transform<5>(skbuffer, 495, cbuffer, 3120, 0, 5);

            t3cfunc::bra_transform<5>(skbuffer, 726, cbuffer, 3561, 0, 5);

            t3cfunc::bra_transform<5>(skbuffer, 957, cbuffer, 4002, 0, 5);

            t3cfunc::bra_transform<5>(skbuffer, 1188, cbuffer, 4443, 0, 6);

            t3cfunc::bra_transform<5>(skbuffer, 1496, cbuffer, 5031, 0, 6);

            t3cfunc::bra_transform<5>(skbuffer, 1804, cbuffer, 5619, 0, 6);

            t3cfunc::bra_transform<5>(skbuffer, 2112, cbuffer, 6207, 0, 7);

            t3cfunc::bra_transform<5>(skbuffer, 2508, cbuffer, 6963, 0, 7);

            t3cfunc::bra_transform<5>(skbuffer, 2904, cbuffer, 7719, 0, 7);

            t3cfunc::bra_transform<5>(skbuffer, 3300, cbuffer, 8475, 0, 8);

            t3cfunc::bra_transform<5>(skbuffer, 3795, cbuffer, 9420, 0, 8);

            t3cfunc::bra_transform<5>(skbuffer, 4290, cbuffer, 10365, 0, 8);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 4785, 0, 495, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 5280, 165, 726, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 5775, 330, 957, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 6270, 495, 1188, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 6963, 726, 1496, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 7656, 957, 1804, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 8349, 1188, 2112, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 9273, 1496, 2508, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 10197, 1804, 2904, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 11121, 2112, 3300, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 12309, 2508, 3795, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 13497, 2904, 4290, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 14685, 4785, 6270, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 15675, 5280, 6963, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 16665, 5775, 7656, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 17655, 6270, 8349, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 19041, 6963, 9273, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 20427, 7656, 10197, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 21813, 8349, 11121, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 23661, 9273, 12309, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 25509, 10197, 13497, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 27357, 14685, 17655, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 29007, 15675, 19041, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 30657, 16665, 20427, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 32307, 17655, 21813, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 34617, 19041, 23661, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 36927, 20427, 25509, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 39237, 27357, 32307, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 41712, 29007, 34617, cfactors, 6, 5);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 44187, 30657, 36927, cfactors, 6, 5);

            t3cfunc::ket_transform<4, 4>(sbuffer, 0, skbuffer, 39237, 5);

            t3cfunc::ket_transform<4, 4>(sbuffer, 891, skbuffer, 41712, 5);

            t3cfunc::ket_transform<4, 4>(sbuffer, 1782, skbuffer, 44187, 5);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 5, 4, 4, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecHGG_hpp */
