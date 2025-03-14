#ifndef ThreeCenterElectronRepulsionGeom100RecDGG_hpp
#define ThreeCenterElectronRepulsionGeom100RecDGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionGeom100ContrRecDXX.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSL.hpp"
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

/// @brief Computes d^(1)/dA^(1)(D|1/|r-r'||GG)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_dgg(T& distributor,
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

    CSimdArray<double> pbuffer(4558, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(4495, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(21210, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1215, 1);

    // setup Boys fuction data

    const CBoysFunc<11> bf_table;

    CSimdArray<double> bf_data(13, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 12, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 12);

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

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 11, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 14, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 17, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 20, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 23, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 26, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 29, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 32, 7, 8, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 35, 8, 9, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 38, 9, 10, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 41, 0, 1, 11, 14, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 47, 1, 2, 14, 17, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 53, 2, 3, 17, 20, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 59, 3, 4, 20, 23, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 65, 4, 5, 23, 26, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 71, 5, 6, 26, 29, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 77, 6, 7, 29, 32, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 83, 7, 8, 32, 35, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 89, 8, 9, 35, 38, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 95, 11, 14, 41, 47, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 105, 14, 17, 47, 53, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 115, 17, 20, 53, 59, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 125, 20, 23, 59, 65, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 135, 23, 26, 65, 71, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 145, 26, 29, 71, 77, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 155, 29, 32, 77, 83, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 165, 32, 35, 83, 89, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 175, 41, 47, 95, 105, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 190, 47, 53, 105, 115, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 205, 53, 59, 115, 125, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 220, 59, 65, 125, 135, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 235, 65, 71, 135, 145, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 250, 71, 77, 145, 155, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 265, 77, 83, 155, 165, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 280, 95, 105, 175, 190, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 301, 105, 115, 190, 205, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 322, 115, 125, 205, 220, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 343, 125, 135, 220, 235, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 364, 135, 145, 235, 250, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 385, 145, 155, 250, 265, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 406, 175, 190, 280, 301, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 434, 190, 205, 301, 322, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 462, 205, 220, 322, 343, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 490, 220, 235, 343, 364, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 518, 235, 250, 364, 385, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 546, 280, 301, 406, 434, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 582, 301, 322, 434, 462, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 618, 322, 343, 462, 490, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 654, 343, 364, 490, 518, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 690, 406, 434, 546, 582, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 735, 434, 462, 582, 618, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 780, 462, 490, 618, 654, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 825, 17, 53, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 843, 53, 115, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 873, 95, 175, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 918, 105, 190, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 963, 115, 205, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1008, 175, 280, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1071, 190, 301, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1134, 205, 322, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1197, 280, 406, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1281, 301, 434, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1365, 322, 462, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 1449, 406, 546, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 1557, 434, 582, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 1665, 462, 618, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 1773, 546, 690, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 1908, 582, 735, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 2043, 618, 780, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 2178, 95, 105, 825, 843, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 2238, 175, 190, 843, 963, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 2328, 280, 301, 963, 1134, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 2454, 406, 434, 1134, 1365, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 2622, 546, 582, 1365, 1665, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 2838, 690, 735, 1665, 2043, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 3108, 873, 918, 2178, 2238, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 3258, 1008, 1071, 2238, 2328, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 3468, 1197, 1281, 2328, 2454, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 3748, 1449, 1557, 2454, 2622, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 4108, 1773, 1908, 2622, 2838, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 873, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 45, pbuffer, 1008, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 108, pbuffer, 1197, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 192, pbuffer, 1449, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 300, pbuffer, 1773, 135, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {3108, 3258});

                pbuffer.scale(2.0 * a_exp, {3258, 3468});

                pbuffer.scale(2.0 * a_exp, {3468, 3748});

                pbuffer.scale(2.0 * a_exp, {3748, 4108});

                pbuffer.scale(2.0 * a_exp, {4108, 4558});

                t2cfunc::reduce(cbuffer, 3045, pbuffer, 3108, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3195, pbuffer, 3258, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3405, pbuffer, 3468, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3685, pbuffer, 3748, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4045, pbuffer, 4108, 450, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_dxx(cbuffer, 435, 0, 3045, 0, 4);

            t3ceri::comp_bra_geom1_electron_repulsion_dxx(cbuffer, 705, 45, 3195, 0, 5);

            t3ceri::comp_bra_geom1_electron_repulsion_dxx(cbuffer, 1083, 108, 3405, 0, 6);

            t3ceri::comp_bra_geom1_electron_repulsion_dxx(cbuffer, 1587, 192, 3685, 0, 7);

            t3ceri::comp_bra_geom1_electron_repulsion_dxx(cbuffer, 2235, 300, 4045, 0, 8);

            t3cfunc::bra_transform<2>(skbuffer, 0, cbuffer, 435, 0, 4);

            t3cfunc::bra_transform<2>(skbuffer, 75, cbuffer, 525, 0, 4);

            t3cfunc::bra_transform<2>(skbuffer, 150, cbuffer, 615, 0, 4);

            t3cfunc::bra_transform<2>(skbuffer, 225, cbuffer, 705, 0, 5);

            t3cfunc::bra_transform<2>(skbuffer, 330, cbuffer, 831, 0, 5);

            t3cfunc::bra_transform<2>(skbuffer, 435, cbuffer, 957, 0, 5);

            t3cfunc::bra_transform<2>(skbuffer, 540, cbuffer, 1083, 0, 6);

            t3cfunc::bra_transform<2>(skbuffer, 680, cbuffer, 1251, 0, 6);

            t3cfunc::bra_transform<2>(skbuffer, 820, cbuffer, 1419, 0, 6);

            t3cfunc::bra_transform<2>(skbuffer, 960, cbuffer, 1587, 0, 7);

            t3cfunc::bra_transform<2>(skbuffer, 1140, cbuffer, 1803, 0, 7);

            t3cfunc::bra_transform<2>(skbuffer, 1320, cbuffer, 2019, 0, 7);

            t3cfunc::bra_transform<2>(skbuffer, 1500, cbuffer, 2235, 0, 8);

            t3cfunc::bra_transform<2>(skbuffer, 1725, cbuffer, 2505, 0, 8);

            t3cfunc::bra_transform<2>(skbuffer, 1950, cbuffer, 2775, 0, 8);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 2175, 0, 225, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 2400, 75, 330, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 2625, 150, 435, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 2850, 225, 540, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 3165, 330, 680, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 3480, 435, 820, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 3795, 540, 960, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 4215, 680, 1140, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 4635, 820, 1320, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 5055, 960, 1500, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 5595, 1140, 1725, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 6135, 1320, 1950, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 6675, 2175, 2850, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 7125, 2400, 3165, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 7575, 2625, 3480, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 8025, 2850, 3795, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 8655, 3165, 4215, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 9285, 3480, 4635, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 9915, 3795, 5055, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 10755, 4215, 5595, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdi(skbuffer, 11595, 4635, 6135, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 12435, 6675, 8025, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 13185, 7125, 8655, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 13935, 7575, 9285, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 14685, 8025, 9915, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 15735, 8655, 10755, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xfh(skbuffer, 16785, 9285, 11595, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 17835, 12435, 14685, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 18960, 13185, 15735, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xgg(skbuffer, 20085, 13935, 16785, cfactors, 6, 2);

            t3cfunc::ket_transform<4, 4>(sbuffer, 0, skbuffer, 17835, 2);

            t3cfunc::ket_transform<4, 4>(sbuffer, 405, skbuffer, 18960, 2);

            t3cfunc::ket_transform<4, 4>(sbuffer, 810, skbuffer, 20085, 2);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 2, 4, 4, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecDGG_hpp */
