#ifndef ThreeCenterElectronRepulsionGeom010RecDGG_hpp
#define ThreeCenterElectronRepulsionGeom010RecDGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPI.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDH.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDI.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFH.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXGG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPI.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPK.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSH.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSI.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSK.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSM.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSM.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSM.hpp"
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

/// @brief Computes d^(1)/dC^(1)(D|1/|r-r'||GG)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_dgg(T& distributor,
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

    CSimdArray<double> pbuffer(2985, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1800, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(25500, 1);

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

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 12, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 15, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 18, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 21, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 24, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 27, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 30, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 33, 7, 8, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 36, 8, 9, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 39, 9, 10, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 42, 10, 11, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 45, 0, 1, 12, 15, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 51, 1, 2, 15, 18, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 57, 2, 3, 18, 21, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 63, 3, 4, 21, 24, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 69, 4, 5, 24, 27, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 75, 5, 6, 27, 30, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 81, 6, 7, 30, 33, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 87, 7, 8, 33, 36, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 93, 8, 9, 36, 39, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 99, 9, 10, 39, 42, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 105, 12, 15, 45, 51, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 115, 15, 18, 51, 57, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 125, 18, 21, 57, 63, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 135, 21, 24, 63, 69, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 145, 24, 27, 69, 75, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 155, 27, 30, 75, 81, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 165, 30, 33, 81, 87, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 175, 33, 36, 87, 93, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 185, 36, 39, 93, 99, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 195, 45, 51, 105, 115, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 210, 51, 57, 115, 125, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 225, 57, 63, 125, 135, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 240, 63, 69, 135, 145, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 255, 69, 75, 145, 155, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 270, 75, 81, 155, 165, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 285, 81, 87, 165, 175, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 300, 87, 93, 175, 185, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 315, 105, 115, 195, 210, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 336, 115, 125, 210, 225, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 357, 125, 135, 225, 240, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 378, 135, 145, 240, 255, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 399, 145, 155, 255, 270, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 420, 155, 165, 270, 285, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 441, 165, 175, 285, 300, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 462, 195, 210, 315, 336, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 490, 210, 225, 336, 357, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 518, 225, 240, 357, 378, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 546, 240, 255, 378, 399, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 574, 255, 270, 399, 420, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 602, 270, 285, 420, 441, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 630, 315, 336, 462, 490, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 666, 336, 357, 490, 518, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 702, 357, 378, 518, 546, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 738, 378, 399, 546, 574, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 774, 399, 420, 574, 602, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 810, 462, 490, 630, 666, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 855, 490, 518, 666, 702, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 900, 518, 546, 702, 738, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 945, 546, 574, 738, 774, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssm(pbuffer, 990, 630, 666, 810, 855, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssm(pbuffer, 1045, 666, 702, 855, 900, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssm(pbuffer, 1100, 702, 738, 900, 945, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1155, 57, 125, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1185, 125, 225, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1230, 225, 357, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1293, 357, 518, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 1377, 518, 702, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 1485, 702, 900, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psm(pbuffer, 1620, 900, 1100, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1785, 195, 210, 1155, 1185, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 1875, 315, 336, 1185, 1230, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 2001, 462, 490, 1230, 1293, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 2169, 630, 666, 1293, 1377, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 2385, 810, 855, 1377, 1485, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsm(pbuffer, 2655, 990, 1045, 1485, 1620, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1785, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 1875, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 216, pbuffer, 2001, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 384, pbuffer, 2169, 216, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {1785, 1875});

                pbuffer.scale(pfactors, 0, 2.0, {1875, 2001});

                pbuffer.scale(pfactors, 0, 2.0, {2001, 2169});

                pbuffer.scale(pfactors, 0, 2.0, {2169, 2385});

                pbuffer.scale(pfactors, 0, 2.0, {2385, 2655});

                pbuffer.scale(pfactors, 0, 2.0, {2655, 2985});

                t2cfunc::reduce(cbuffer, 600, pbuffer, 1785, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 690, pbuffer, 1875, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 816, pbuffer, 2001, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 984, pbuffer, 2169, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1200, pbuffer, 2385, 270, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1470, pbuffer, 2655, 330, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<2>(skbuffer, 0, cbuffer, 0, 0, 4);

            t3cfunc::bra_transform<2>(skbuffer, 75, cbuffer, 90, 0, 5);

            t3cfunc::bra_transform<2>(skbuffer, 180, cbuffer, 216, 0, 6);

            t3cfunc::bra_transform<2>(skbuffer, 320, cbuffer, 384, 0, 7);

            t3cfunc::bra_transform<2>(skbuffer, 22325, cbuffer, 600, 0, 4);

            t3cfunc::bra_transform<2>(skbuffer, 22400, cbuffer, 690, 0, 5);

            t3cfunc::bra_transform<2>(skbuffer, 22505, cbuffer, 816, 0, 6);

            t3cfunc::bra_transform<2>(skbuffer, 22645, cbuffer, 984, 0, 7);

            t3cfunc::bra_transform<2>(skbuffer, 22825, cbuffer, 1200, 0, 8);

            t3cfunc::bra_transform<2>(skbuffer, 23050, cbuffer, 1470, 0, 9);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 500, 0, 75, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 1400, 75, 180, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 2660, 180, 320, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 5960, 500, 1400, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 7760, 1400, 2660, cfactors, 6, 2);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 12800, 5960, 7760, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xsg(skbuffer, 23325, 22325, 22400, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xsh(skbuffer, 23550, 22400, 22505, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xsi(skbuffer, 23865, 22505, 22645, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xsk(skbuffer, 24285, 22645, 22825, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xsl(skbuffer, 24825, 22825, 23050, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xpg(skbuffer, 725, 0, 23325, 23550, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xph(skbuffer, 1715, 75, 23550, 23865, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xpi(skbuffer, 3080, 180, 23865, 24285, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xpk(skbuffer, 4340, 320, 24285, 24825, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xdg(skbuffer, 6410, 500, 725, 1715, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xdh(skbuffer, 8390, 1400, 1715, 3080, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xdi(skbuffer, 10280, 2660, 3080, 4340, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xfg(skbuffer, 13550, 5960, 6410, 8390, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xfh(skbuffer, 15800, 7760, 8390, 10280, cfactors, 6, 2);

            t3ceri::comp_ket_geom010_electron_repulsion_xgg(skbuffer, 18950, 12800, 13550, 15800, cfactors, 6, 2);

            t3cfunc::ket_transform<4, 4>(sbuffer, 0, skbuffer, 18950, 2);

            t3cfunc::ket_transform<4, 4>(sbuffer, 405, skbuffer, 20075, 2);

            t3cfunc::ket_transform<4, 4>(sbuffer, 810, skbuffer, 21200, 2);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 2, 4, 4, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecDGG_hpp */
