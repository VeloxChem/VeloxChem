#ifndef ThreeCenterElectronRepulsionGeom100RecDFF_hpp
#define ThreeCenterElectronRepulsionGeom100RecDFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionGeom100ContrRecDXX.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
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

/// @brief Computes d^(1)/dA^(1)(D|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_dff(T& distributor,
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

    CSimdArray<double> pbuffer(2291, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2294, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(6930, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(735, 1);

    // setup Boys fuction data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_data(11, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 10);

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

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 9, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 12, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 15, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 18, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 21, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 24, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 27, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 30, 7, 8, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 33, 0, 1, 9, 12, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 39, 1, 2, 12, 15, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 45, 2, 3, 15, 18, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 51, 3, 4, 18, 21, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 57, 4, 5, 21, 24, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 63, 5, 6, 24, 27, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 69, 6, 7, 27, 30, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 75, 9, 12, 33, 39, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 85, 12, 15, 39, 45, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 95, 15, 18, 45, 51, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 105, 18, 21, 51, 57, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 115, 21, 24, 57, 63, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 125, 24, 27, 63, 69, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 135, 33, 39, 75, 85, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 150, 39, 45, 85, 95, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 165, 45, 51, 95, 105, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 180, 51, 57, 105, 115, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 195, 57, 63, 115, 125, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 210, 75, 85, 135, 150, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 231, 85, 95, 150, 165, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 252, 95, 105, 165, 180, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 273, 105, 115, 180, 195, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 294, 135, 150, 210, 231, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 322, 150, 165, 231, 252, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 350, 165, 180, 252, 273, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 378, 2, 15, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 387, 15, 45, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 405, 33, 75, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 435, 39, 85, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 465, 45, 95, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 495, 75, 135, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 540, 85, 150, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 585, 95, 165, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 630, 135, 210, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 693, 150, 231, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 756, 165, 252, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 819, 210, 294, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 903, 231, 322, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 987, 252, 350, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 1071, 33, 39, 378, 387, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 1107, 75, 85, 387, 465, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1167, 135, 150, 465, 585, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 1257, 210, 231, 585, 756, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 1383, 294, 322, 756, 987, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 1551, 405, 435, 1071, 1107, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 1651, 495, 540, 1107, 1167, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 1801, 630, 693, 1167, 1257, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 2011, 819, 903, 1257, 1383, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 405, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 495, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 75, pbuffer, 630, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 138, pbuffer, 819, 84, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1551, 1651});

                pbuffer.scale(2.0 * a_exp, {1651, 1801});

                pbuffer.scale(2.0 * a_exp, {1801, 2011});

                pbuffer.scale(2.0 * a_exp, {2011, 2291});

                t2cfunc::reduce(cbuffer, 1554, pbuffer, 1551, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1654, pbuffer, 1651, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1804, pbuffer, 1801, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2014, pbuffer, 2011, 280, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_dxx(cbuffer, 222, 0, 1554, 0, 3);

            t3ceri::comp_bra_geom1_electron_repulsion_dxx(cbuffer, 402, 30, 1654, 0, 4);

            t3ceri::comp_bra_geom1_electron_repulsion_dxx(cbuffer, 672, 75, 1804, 0, 5);

            t3ceri::comp_bra_geom1_electron_repulsion_dxx(cbuffer, 1050, 138, 2014, 0, 6);

            t3cfunc::bra_transform<2>(skbuffer, 0, cbuffer, 222, 0, 3);

            t3cfunc::bra_transform<2>(skbuffer, 50, cbuffer, 282, 0, 3);

            t3cfunc::bra_transform<2>(skbuffer, 100, cbuffer, 342, 0, 3);

            t3cfunc::bra_transform<2>(skbuffer, 150, cbuffer, 402, 0, 4);

            t3cfunc::bra_transform<2>(skbuffer, 225, cbuffer, 492, 0, 4);

            t3cfunc::bra_transform<2>(skbuffer, 300, cbuffer, 582, 0, 4);

            t3cfunc::bra_transform<2>(skbuffer, 375, cbuffer, 672, 0, 5);

            t3cfunc::bra_transform<2>(skbuffer, 480, cbuffer, 798, 0, 5);

            t3cfunc::bra_transform<2>(skbuffer, 585, cbuffer, 924, 0, 5);

            t3cfunc::bra_transform<2>(skbuffer, 690, cbuffer, 1050, 0, 6);

            t3cfunc::bra_transform<2>(skbuffer, 830, cbuffer, 1218, 0, 6);

            t3cfunc::bra_transform<2>(skbuffer, 970, cbuffer, 1386, 0, 6);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 1110, 0, 150, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 1260, 50, 225, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 1410, 100, 300, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 1560, 150, 375, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 1785, 225, 480, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 2010, 300, 585, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 2235, 375, 690, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 2550, 480, 830, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 2865, 585, 970, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdf(skbuffer, 3180, 1110, 1560, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdf(skbuffer, 3480, 1260, 1785, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdf(skbuffer, 3780, 1410, 2010, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 4080, 1560, 2235, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 4530, 1785, 2550, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 4980, 2010, 2865, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xff(skbuffer, 5430, 3180, 4080, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xff(skbuffer, 5930, 3480, 4530, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xff(skbuffer, 6430, 3780, 4980, cfactors, 6, 2);

            t3cfunc::ket_transform<3, 3>(sbuffer, 0, skbuffer, 5430, 2);

            t3cfunc::ket_transform<3, 3>(sbuffer, 245, skbuffer, 5930, 2);

            t3cfunc::ket_transform<3, 3>(sbuffer, 490, skbuffer, 6430, 2);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 2, 3, 3, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecDFF_hpp */
