#ifndef ThreeCenterElectronRepulsionGeom010RecFGG_hpp
#define ThreeCenterElectronRepulsionGeom010RecFGG_hpp

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
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSM.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSI.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSK.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSL.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSM.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
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

/// @brief Computes d^(1)/dC^(1)(F|1/|r-r'||GG)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_fgg(T& distributor,
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

    CSimdArray<double> pbuffer(6263, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(3000, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(35700, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1701, 1);

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

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 1155, 18, 57, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 1173, 57, 125, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1203, 105, 195, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1248, 115, 210, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 1293, 125, 225, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1338, 195, 315, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1401, 210, 336, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 1464, 225, 357, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1527, 315, 462, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1611, 336, 490, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psi(pbuffer, 1695, 357, 518, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 1779, 462, 630, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 1887, 490, 666, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psk(pbuffer, 1995, 518, 702, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 2103, 630, 810, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 2238, 666, 855, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psl(pbuffer, 2373, 702, 900, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psm(pbuffer, 2508, 810, 990, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psm(pbuffer, 2673, 855, 1045, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psm(pbuffer, 2838, 900, 1100, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 3003, 105, 115, 1155, 1173, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 3063, 195, 210, 1173, 1293, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 3153, 315, 336, 1293, 1464, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsi(pbuffer, 3279, 462, 490, 1464, 1695, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsk(pbuffer, 3447, 630, 666, 1695, 1995, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsl(pbuffer, 3663, 810, 855, 1995, 2373, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsm(pbuffer, 3933, 990, 1045, 2373, 2838, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 4263, 1203, 1248, 3003, 3063, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 4413, 1338, 1401, 3063, 3153, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsi(pbuffer, 4623, 1527, 1611, 3153, 3279, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsk(pbuffer, 4903, 1779, 1887, 3279, 3447, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsl(pbuffer, 5263, 2103, 2238, 3447, 3663, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsm(pbuffer, 5713, 2508, 2673, 3663, 3933, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 4263, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 150, pbuffer, 4413, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 360, pbuffer, 4623, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 640, pbuffer, 4903, 360, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {4263, 4413});

                pbuffer.scale(pfactors, 0, 2.0, {4413, 4623});

                pbuffer.scale(pfactors, 0, 2.0, {4623, 4903});

                pbuffer.scale(pfactors, 0, 2.0, {4903, 5263});

                pbuffer.scale(pfactors, 0, 2.0, {5263, 5713});

                pbuffer.scale(pfactors, 0, 2.0, {5713, 6263});

                t2cfunc::reduce(cbuffer, 1000, pbuffer, 4263, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1150, pbuffer, 4413, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1360, pbuffer, 4623, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1640, pbuffer, 4903, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2000, pbuffer, 5263, 450, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2450, pbuffer, 5713, 550, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<3>(skbuffer, 0, cbuffer, 0, 0, 4);

            t3cfunc::bra_transform<3>(skbuffer, 105, cbuffer, 150, 0, 5);

            t3cfunc::bra_transform<3>(skbuffer, 252, cbuffer, 360, 0, 6);

            t3cfunc::bra_transform<3>(skbuffer, 448, cbuffer, 640, 0, 7);

            t3cfunc::bra_transform<3>(skbuffer, 31255, cbuffer, 1000, 0, 4);

            t3cfunc::bra_transform<3>(skbuffer, 31360, cbuffer, 1150, 0, 5);

            t3cfunc::bra_transform<3>(skbuffer, 31507, cbuffer, 1360, 0, 6);

            t3cfunc::bra_transform<3>(skbuffer, 31703, cbuffer, 1640, 0, 7);

            t3cfunc::bra_transform<3>(skbuffer, 31955, cbuffer, 2000, 0, 8);

            t3cfunc::bra_transform<3>(skbuffer, 32270, cbuffer, 2450, 0, 9);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 700, 0, 105, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 1960, 105, 252, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 3724, 252, 448, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 8344, 700, 1960, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 10864, 1960, 3724, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 17920, 8344, 10864, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsg(skbuffer, 32655, 31255, 31360, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsh(skbuffer, 32970, 31360, 31507, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsi(skbuffer, 33411, 31507, 31703, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsk(skbuffer, 33999, 31703, 31955, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsl(skbuffer, 34755, 31955, 32270, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpg(skbuffer, 1015, 0, 32655, 32970, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xph(skbuffer, 2401, 105, 32970, 33411, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpi(skbuffer, 4312, 252, 33411, 33999, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpk(skbuffer, 6076, 448, 33999, 34755, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xdg(skbuffer, 8974, 700, 1015, 2401, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xdh(skbuffer, 11746, 1960, 2401, 4312, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xdi(skbuffer, 14392, 3724, 4312, 6076, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xfg(skbuffer, 18970, 8344, 8974, 11746, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xfh(skbuffer, 22120, 10864, 11746, 14392, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xgg(skbuffer, 26530, 17920, 18970, 22120, cfactors, 6, 3);

            t3cfunc::ket_transform<4, 4>(sbuffer, 0, skbuffer, 26530, 3);

            t3cfunc::ket_transform<4, 4>(sbuffer, 567, skbuffer, 28105, 3);

            t3cfunc::ket_transform<4, 4>(sbuffer, 1134, skbuffer, 29680, 3);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 3, 4, 4, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecFGG_hpp */
