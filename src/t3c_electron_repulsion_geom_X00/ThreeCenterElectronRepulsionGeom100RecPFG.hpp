#ifndef ThreeCenterElectronRepulsionGeom100RecPFG_hpp
#define ThreeCenterElectronRepulsionGeom100RecPFG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionGeom100ContrRecPXX.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSK.hpp"
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

/// @brief Computes d^(1)/dA^(1)(P|1/|r-r'||FG)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_pfg(T& distributor,
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

    CSimdArray<double> pbuffer(1500, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1600, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(5922, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(567, 1);

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

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 10, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 13, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 16, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 19, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 22, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 25, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 28, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 31, 7, 8, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 34, 8, 9, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 37, 0, 1, 10, 13, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 43, 1, 2, 13, 16, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 49, 2, 3, 16, 19, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 55, 3, 4, 19, 22, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 61, 4, 5, 22, 25, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 67, 5, 6, 25, 28, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 73, 6, 7, 28, 31, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 79, 7, 8, 31, 34, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 85, 10, 13, 37, 43, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 95, 13, 16, 43, 49, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 105, 16, 19, 49, 55, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 115, 19, 22, 55, 61, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 125, 22, 25, 61, 67, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 135, 25, 28, 67, 73, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 145, 28, 31, 73, 79, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 155, 37, 43, 85, 95, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 170, 43, 49, 95, 105, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 185, 49, 55, 105, 115, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 200, 55, 61, 115, 125, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 215, 61, 67, 125, 135, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 230, 67, 73, 135, 145, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 245, 85, 95, 155, 170, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 266, 95, 105, 170, 185, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 287, 105, 115, 185, 200, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 308, 115, 125, 200, 215, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 329, 125, 135, 215, 230, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 350, 155, 170, 245, 266, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 378, 170, 185, 266, 287, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 406, 185, 200, 287, 308, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 434, 200, 215, 308, 329, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 462, 245, 266, 350, 378, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 498, 266, 287, 378, 406, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 534, 287, 308, 406, 434, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 570, 49, 105, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 600, 105, 185, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 645, 185, 287, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 708, 287, 406, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 792, 406, 534, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 900, 155, 170, 570, 600, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 990, 245, 266, 600, 645, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 1116, 350, 378, 645, 708, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 1284, 462, 498, 708, 792, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15, pbuffer, 245, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 350, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 64, pbuffer, 462, 36, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {900, 990});

                pbuffer.scale(2.0 * a_exp, {990, 1116});

                pbuffer.scale(2.0 * a_exp, {1116, 1284});

                pbuffer.scale(2.0 * a_exp, {1284, 1500});

                t2cfunc::reduce(cbuffer, 1000, pbuffer, 900, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1090, pbuffer, 990, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1216, pbuffer, 1116, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1384, pbuffer, 1284, 216, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_pxx(cbuffer, 100, 0, 1000, 0, 4);

            t3ceri::comp_bra_geom1_electron_repulsion_pxx(cbuffer, 235, 15, 1090, 0, 5);

            t3ceri::comp_bra_geom1_electron_repulsion_pxx(cbuffer, 424, 36, 1216, 0, 6);

            t3ceri::comp_bra_geom1_electron_repulsion_pxx(cbuffer, 676, 64, 1384, 0, 7);

            t3cfunc::bra_transform<1>(skbuffer, 0, cbuffer, 100, 0, 4);

            t3cfunc::bra_transform<1>(skbuffer, 45, cbuffer, 145, 0, 4);

            t3cfunc::bra_transform<1>(skbuffer, 90, cbuffer, 190, 0, 4);

            t3cfunc::bra_transform<1>(skbuffer, 135, cbuffer, 235, 0, 5);

            t3cfunc::bra_transform<1>(skbuffer, 198, cbuffer, 298, 0, 5);

            t3cfunc::bra_transform<1>(skbuffer, 261, cbuffer, 361, 0, 5);

            t3cfunc::bra_transform<1>(skbuffer, 324, cbuffer, 424, 0, 6);

            t3cfunc::bra_transform<1>(skbuffer, 408, cbuffer, 508, 0, 6);

            t3cfunc::bra_transform<1>(skbuffer, 492, cbuffer, 592, 0, 6);

            t3cfunc::bra_transform<1>(skbuffer, 576, cbuffer, 676, 0, 7);

            t3cfunc::bra_transform<1>(skbuffer, 684, cbuffer, 784, 0, 7);

            t3cfunc::bra_transform<1>(skbuffer, 792, cbuffer, 892, 0, 7);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 900, 0, 135, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 1035, 45, 198, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 1170, 90, 261, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 1305, 135, 324, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 1494, 198, 408, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 1683, 261, 492, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 1872, 324, 576, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 2124, 408, 684, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 2376, 492, 792, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 2628, 900, 1305, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 2898, 1035, 1494, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 3168, 1170, 1683, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 3438, 1305, 1872, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 3816, 1494, 2124, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 4194, 1683, 2376, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 4572, 2628, 3438, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 5022, 2898, 3816, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 5472, 3168, 4194, cfactors, 6, 1);

            t3cfunc::ket_transform<3, 4>(sbuffer, 0, skbuffer, 4572, 1);

            t3cfunc::ket_transform<3, 4>(sbuffer, 189, skbuffer, 5022, 1);

            t3cfunc::ket_transform<3, 4>(sbuffer, 378, skbuffer, 5472, 1);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 1, 3, 4, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecPFG_hpp */
