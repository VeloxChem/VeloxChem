#ifndef ThreeCenterElectronRepulsionGeom010RecDGD_hpp
#define ThreeCenterElectronRepulsionGeom010RecDGD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXGD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSH.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
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

/// @brief Computes d^(1)/dC^(1)(D|1/|r-r'||GD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_dgd(T& distributor,
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

    CSimdArray<double> pbuffer(1623, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1008, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(12165, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(675, 1);

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

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 570, 2, 16, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 579, 16, 49, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 597, 49, 105, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 627, 105, 185, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 672, 185, 287, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 735, 287, 406, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 819, 406, 534, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 927, 37, 43, 570, 579, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 963, 85, 95, 579, 597, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1023, 155, 170, 597, 627, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 1113, 245, 266, 627, 672, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 1239, 350, 378, 672, 735, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 1407, 462, 498, 735, 819, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 927, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 963, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 96, pbuffer, 1023, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 186, pbuffer, 1113, 126, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {927, 963});

                pbuffer.scale(pfactors, 0, 2.0, {963, 1023});

                pbuffer.scale(pfactors, 0, 2.0, {1023, 1113});

                pbuffer.scale(pfactors, 0, 2.0, {1113, 1239});

                pbuffer.scale(pfactors, 0, 2.0, {1239, 1407});

                pbuffer.scale(pfactors, 0, 2.0, {1407, 1623});

                t2cfunc::reduce(cbuffer, 312, pbuffer, 927, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 348, pbuffer, 963, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 408, pbuffer, 1023, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 498, pbuffer, 1113, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 624, pbuffer, 1239, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 792, pbuffer, 1407, 216, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<2>(skbuffer, 0, cbuffer, 0, 0, 2);

            t3cfunc::bra_transform<2>(skbuffer, 30, cbuffer, 36, 0, 3);

            t3cfunc::bra_transform<2>(skbuffer, 80, cbuffer, 96, 0, 4);

            t3cfunc::bra_transform<2>(skbuffer, 155, cbuffer, 186, 0, 5);

            t3cfunc::bra_transform<2>(skbuffer, 10385, cbuffer, 312, 0, 2);

            t3cfunc::bra_transform<2>(skbuffer, 10415, cbuffer, 348, 0, 3);

            t3cfunc::bra_transform<2>(skbuffer, 10465, cbuffer, 408, 0, 4);

            t3cfunc::bra_transform<2>(skbuffer, 10540, cbuffer, 498, 0, 5);

            t3cfunc::bra_transform<2>(skbuffer, 10645, cbuffer, 624, 0, 6);

            t3cfunc::bra_transform<2>(skbuffer, 10785, cbuffer, 792, 0, 7);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 260, 0, 30, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 620, 30, 80, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 1220, 80, 155, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdd(skbuffer, 3065, 260, 620, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdf(skbuffer, 3785, 620, 1220, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xfd(skbuffer, 6335, 3065, 3785, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xsd(skbuffer, 10965, 10385, 10415, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xsf(skbuffer, 11055, 10415, 10465, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xsg(skbuffer, 11205, 10465, 10540, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xsh(skbuffer, 11430, 10540, 10645, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xsi(skbuffer, 11745, 10645, 10785, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xpd(skbuffer, 350, 0, 10965, 11055, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xpf(skbuffer, 770, 30, 11055, 11205, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xpg(skbuffer, 1445, 80, 11205, 11430, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xph(skbuffer, 2120, 155, 11430, 11745, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xdd(skbuffer, 3245, 260, 350, 770, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xdf(skbuffer, 4085, 620, 770, 1445, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xdg(skbuffer, 4985, 1220, 1445, 2120, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xfd(skbuffer, 6635, 3065, 3245, 4085, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xff(skbuffer, 7535, 3785, 4085, 4985, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xgd(skbuffer, 9035, 6335, 6635, 7535, cfactors, 6, 2);

            t3cfunc::ket_transform<4, 2>(sbuffer, 0, skbuffer, 9035, 2);

            t3cfunc::ket_transform<4, 2>(sbuffer, 225, skbuffer, 9485, 2);

            t3cfunc::ket_transform<4, 2>(sbuffer, 450, skbuffer, 9935, 2);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 2, 4, 2, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecDGD_hpp */
