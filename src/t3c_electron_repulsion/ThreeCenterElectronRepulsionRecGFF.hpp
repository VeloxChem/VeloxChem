#ifndef ThreeCenterElectronRepulsionRecGFF_hpp
#define ThreeCenterElectronRepulsionRecGFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSI.hpp"
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
#include "GtoPairBlock.hpp"
#include "GtoBlock.hpp"
#include "BatchFunc.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes (G|1/|r-r'||FF)  integrals for basis functions block and basis function pairs block.
/// @param distributor The pointer to integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_gff(T& distributor,
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

    CSimdArray<double> pbuffer(4574, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1110, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(4158, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(441, 1);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 11);

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

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 546, 4, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 549, 4, 23, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 558, 17, 53, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 576, 20, 59, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 594, 23, 65, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 612, 53, 115, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 642, 59, 125, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 672, 65, 135, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 702, 115, 205, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 747, 125, 220, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 792, 135, 235, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 837, 205, 322, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 900, 220, 343, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 963, 235, 364, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1026, 322, 462, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1110, 343, 490, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1194, 364, 518, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 1278, 17, 20, 546, 549, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 1296, 53, 59, 549, 594, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 1332, 95, 105, 558, 612, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 1392, 105, 115, 576, 642, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 1452, 115, 125, 594, 672, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1512, 175, 190, 612, 702, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1602, 190, 205, 642, 747, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1692, 205, 220, 672, 792, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 1782, 280, 301, 702, 837, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 1908, 301, 322, 747, 900, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 2034, 322, 343, 792, 963, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 2160, 406, 434, 837, 1026, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 2328, 434, 462, 900, 1110, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 2496, 462, 490, 963, 1194, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 2664, 558, 576, 1278, 1296, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 2724, 612, 642, 1296, 1452, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 2824, 702, 747, 1452, 1692, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 2974, 837, 900, 1692, 2034, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 3184, 1026, 1110, 2034, 2496, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 3464, 1332, 1392, 2664, 2724, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 3614, 1512, 1602, 2724, 2824, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsh(pbuffer, 3839, 1782, 1908, 2824, 2974, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsi(pbuffer, 4154, 2160, 2328, 2974, 3184, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 3464, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 150, pbuffer, 3614, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 375, pbuffer, 3839, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 690, pbuffer, 4154, 420, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<4>(skbuffer, 0, cbuffer, 0, 0, 3);

            t3cfunc::bra_transform<4>(skbuffer, 90, cbuffer, 150, 0, 4);

            t3cfunc::bra_transform<4>(skbuffer, 225, cbuffer, 375, 0, 5);

            t3cfunc::bra_transform<4>(skbuffer, 414, cbuffer, 690, 0, 6);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 666, 0, 90, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 936, 90, 225, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 1341, 225, 414, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xdf(skbuffer, 1908, 666, 936, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 2448, 936, 1341, cfactors, 6, 4);

            t3ceri::comp_hrr_electron_repulsion_xff(skbuffer, 3258, 1908, 2448, cfactors, 6, 4);

            t3cfunc::ket_transform<3, 3>(sbuffer, 0, skbuffer, 3258, 4);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 4, 3, 3, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionRecGFF_hpp */
