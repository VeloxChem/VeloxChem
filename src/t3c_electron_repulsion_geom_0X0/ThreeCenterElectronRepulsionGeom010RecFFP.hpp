#ifndef ThreeCenterElectronRepulsionGeom010RecFFP_hpp
#define ThreeCenterElectronRepulsionGeom010RecFFP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSH.hpp"
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

/// @brief Computes d^(1)/dC^(1)(F|1/|r-r'||FP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_ffp(T& distributor,
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

    CSimdArray<double> pbuffer(1622, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(740, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(4508, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(441, 1);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 9);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 6);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 6, pfactors, 16, bf_data, 7);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 7, pfactors, 16, bf_data, 8);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 8, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 11, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 14, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 17, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 20, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 23, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 26, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 29, 0, 1, 8, 11, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 35, 1, 2, 11, 14, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 41, 2, 3, 14, 17, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 47, 3, 4, 17, 20, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 53, 4, 5, 20, 23, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 59, 5, 6, 23, 26, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 65, 8, 11, 29, 35, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 75, 11, 14, 35, 41, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 85, 14, 17, 41, 47, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 95, 17, 20, 47, 53, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 105, 20, 23, 53, 59, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 115, 29, 35, 65, 75, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 130, 35, 41, 75, 85, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 145, 41, 47, 85, 95, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 160, 47, 53, 95, 105, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 175, 65, 75, 115, 130, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 196, 75, 85, 130, 145, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 217, 85, 95, 145, 160, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 238, 2, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 241, 0, 8, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 250, 1, 11, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 259, 2, 14, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 268, 8, 29, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 286, 11, 35, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 304, 14, 41, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 322, 29, 65, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 352, 35, 75, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 382, 41, 85, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 412, 65, 115, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 457, 75, 130, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 502, 85, 145, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 547, 115, 175, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 610, 130, 196, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psh(pbuffer, 673, 145, 217, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 736, 0, 1, 238, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 742, 8, 11, 238, 259, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 760, 29, 35, 259, 304, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 796, 65, 75, 304, 382, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 856, 115, 130, 382, 502, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsh(pbuffer, 946, 175, 196, 502, 673, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 1072, 241, 250, 736, 742, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 1102, 268, 286, 742, 760, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 1162, 322, 352, 760, 796, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 1262, 412, 457, 796, 856, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsh(pbuffer, 1412, 547, 610, 856, 946, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1072, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 1102, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 1162, 100, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {1072, 1102});

                pbuffer.scale(pfactors, 0, 2.0, {1102, 1162});

                pbuffer.scale(pfactors, 0, 2.0, {1162, 1262});

                pbuffer.scale(pfactors, 0, 2.0, {1262, 1412});

                pbuffer.scale(pfactors, 0, 2.0, {1412, 1622});

                t2cfunc::reduce(cbuffer, 190, pbuffer, 1072, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 220, pbuffer, 1102, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 280, pbuffer, 1162, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 380, pbuffer, 1262, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 530, pbuffer, 1412, 210, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<3>(skbuffer, 0, cbuffer, 0, 0, 1);

            t3cfunc::bra_transform<3>(skbuffer, 21, cbuffer, 30, 0, 2);

            t3cfunc::bra_transform<3>(skbuffer, 63, cbuffer, 90, 0, 3);

            t3cfunc::bra_transform<3>(skbuffer, 3409, cbuffer, 190, 0, 1);

            t3cfunc::bra_transform<3>(skbuffer, 3430, cbuffer, 220, 0, 2);

            t3cfunc::bra_transform<3>(skbuffer, 3472, cbuffer, 280, 0, 3);

            t3cfunc::bra_transform<3>(skbuffer, 3542, cbuffer, 380, 0, 4);

            t3cfunc::bra_transform<3>(skbuffer, 3647, cbuffer, 530, 0, 5);

            t3ceri::comp_hrr_electron_repulsion_xpp(skbuffer, 133, 0, 21, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 385, 21, 63, cfactors, 6, 3);

            t3ceri::comp_hrr_electron_repulsion_xdp(skbuffer, 1519, 133, 385, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsp(skbuffer, 3794, 3409, 3430, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsd(skbuffer, 3857, 3430, 3472, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsf(skbuffer, 3983, 3472, 3542, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xsg(skbuffer, 4193, 3542, 3647, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpp(skbuffer, 196, 0, 3794, 3857, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpd(skbuffer, 511, 21, 3857, 3983, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xpf(skbuffer, 889, 63, 3983, 4193, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xdp(skbuffer, 1645, 133, 196, 511, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xdd(skbuffer, 2023, 385, 511, 889, cfactors, 6, 3);

            t3ceri::comp_ket_geom010_electron_repulsion_xfp(skbuffer, 2779, 1519, 1645, 2023, cfactors, 6, 3);

            t3cfunc::ket_transform<3, 1>(sbuffer, 0, skbuffer, 2779, 3);

            t3cfunc::ket_transform<3, 1>(sbuffer, 147, skbuffer, 2989, 3);

            t3cfunc::ket_transform<3, 1>(sbuffer, 294, skbuffer, 3199, 3);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 3, 3, 1, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecFFP_hpp */
