#ifndef ThreeCenterElectronRepulsionGeom010RecSGG_hpp
#define ThreeCenterElectronRepulsionGeom010RecSGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPI.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPK.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPL.hpp"
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

/// @brief Computes d^(1)/dC^(1)(S|1/|r-r'||GG)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_sgg(T& distributor,
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

    CSimdArray<double> pfactors(26, ket_npgtos);

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(715, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(300, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(5100, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(243, 1);

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

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 570, 350, 378, 462, 498, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssl(pbuffer, 615, 378, 406, 498, 534, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssm(pbuffer, 660, 462, 498, 570, 615, pfactors, 20, 23, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15, pbuffer, 245, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 350, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 64, pbuffer, 462, 36, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {155, 170});

                pbuffer.scale(pfactors, 0, 2.0, {245, 266});

                pbuffer.scale(pfactors, 0, 2.0, {350, 378});

                pbuffer.scale(pfactors, 0, 2.0, {462, 498});

                pbuffer.scale(pfactors, 0, 2.0, {570, 615});

                pbuffer.scale(pfactors, 0, 2.0, {660, 715});

                t2cfunc::reduce(cbuffer, 100, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 115, pbuffer, 245, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 136, pbuffer, 350, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 164, pbuffer, 462, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 200, pbuffer, 570, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 245, pbuffer, 660, 55, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<0>(skbuffer, 0, cbuffer, 0, 0, 4);

            t3cfunc::bra_transform<0>(skbuffer, 15, cbuffer, 15, 0, 5);

            t3cfunc::bra_transform<0>(skbuffer, 36, cbuffer, 36, 0, 6);

            t3cfunc::bra_transform<0>(skbuffer, 64, cbuffer, 64, 0, 7);

            t3cfunc::bra_transform<0>(skbuffer, 4465, cbuffer, 100, 0, 4);

            t3cfunc::bra_transform<0>(skbuffer, 4480, cbuffer, 115, 0, 5);

            t3cfunc::bra_transform<0>(skbuffer, 4501, cbuffer, 136, 0, 6);

            t3cfunc::bra_transform<0>(skbuffer, 4529, cbuffer, 164, 0, 7);

            t3cfunc::bra_transform<0>(skbuffer, 4565, cbuffer, 200, 0, 8);

            t3cfunc::bra_transform<0>(skbuffer, 4610, cbuffer, 245, 0, 9);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 100, 0, 15, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 280, 15, 36, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 532, 36, 64, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xdg(skbuffer, 1192, 100, 280, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xdh(skbuffer, 1552, 280, 532, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xfg(skbuffer, 2560, 1192, 1552, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 4665, 4465, 4480, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 4710, 4480, 4501, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 4773, 4501, 4529, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpk(skbuffer, 4857, 4529, 4565, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpl(skbuffer, 4965, 4565, 4610, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpg(skbuffer, 145, 0, 4665, 4710, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xph(skbuffer, 343, 15, 4710, 4773, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpi(skbuffer, 616, 36, 4773, 4857, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpk(skbuffer, 868, 64, 4857, 4965, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xdg(skbuffer, 1282, 100, 145, 343, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xdh(skbuffer, 1678, 280, 343, 616, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xdi(skbuffer, 2056, 532, 616, 868, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xfg(skbuffer, 2710, 1192, 1282, 1678, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xfh(skbuffer, 3160, 1552, 1678, 2056, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xgg(skbuffer, 3790, 2560, 2710, 3160, cfactors, 6, 0);

            t3cfunc::ket_transform<4, 4>(sbuffer, 0, skbuffer, 3790, 0);

            t3cfunc::ket_transform<4, 4>(sbuffer, 81, skbuffer, 4015, 0);

            t3cfunc::ket_transform<4, 4>(sbuffer, 162, skbuffer, 4240, 0);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 0, 4, 4, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecSGG_hpp */
