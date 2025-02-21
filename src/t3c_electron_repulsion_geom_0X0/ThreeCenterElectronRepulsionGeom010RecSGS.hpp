#ifndef ThreeCenterElectronRepulsionGeom010RecSGS_hpp
#define ThreeCenterElectronRepulsionGeom010RecSGS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDS.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFS.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXGS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPS.hpp"
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

/// @brief Computes d^(1)/dC^(1)(S|1/|r-r'||GS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_sgs(T& distributor,
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

    CSimdArray<double> pbuffer(126, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(76, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(770, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(27, 1);

    // setup Boys fuction data

    const CBoysFunc<5> bf_table;

    CSimdArray<double> bf_data(7, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 6, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 6);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 0);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 6, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 9, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 12, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 15, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 18, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 21, 0, 1, 6, 9, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 27, 1, 2, 9, 12, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 33, 2, 3, 12, 15, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 39, 3, 4, 15, 18, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 45, 6, 9, 21, 27, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 55, 9, 12, 27, 33, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 65, 12, 15, 33, 39, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 75, 21, 27, 45, 55, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 90, 27, 33, 55, 65, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 105, 45, 55, 75, 90, pfactors, 20, 23, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 45, 10, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {0, 1});

                pbuffer.scale(pfactors, 0, 2.0, {6, 9});

                pbuffer.scale(pfactors, 0, 2.0, {21, 27});

                pbuffer.scale(pfactors, 0, 2.0, {45, 55});

                pbuffer.scale(pfactors, 0, 2.0, {75, 90});

                pbuffer.scale(pfactors, 0, 2.0, {105, 126});

                t2cfunc::reduce(cbuffer, 20, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 21, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 24, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 45, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 40, pbuffer, 75, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 55, pbuffer, 105, 21, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<0>(skbuffer, 0, cbuffer, 0, 0, 0);

            t3cfunc::bra_transform<0>(skbuffer, 1, cbuffer, 1, 0, 1);

            t3cfunc::bra_transform<0>(skbuffer, 4, cbuffer, 4, 0, 2);

            t3cfunc::bra_transform<0>(skbuffer, 10, cbuffer, 10, 0, 3);

            t3cfunc::bra_transform<0>(skbuffer, 609, cbuffer, 20, 0, 0);

            t3cfunc::bra_transform<0>(skbuffer, 610, cbuffer, 21, 0, 1);

            t3cfunc::bra_transform<0>(skbuffer, 613, cbuffer, 24, 0, 2);

            t3cfunc::bra_transform<0>(skbuffer, 619, cbuffer, 30, 0, 3);

            t3cfunc::bra_transform<0>(skbuffer, 629, cbuffer, 40, 0, 4);

            t3cfunc::bra_transform<0>(skbuffer, 644, cbuffer, 55, 0, 5);

            t3ceri::comp_hrr_electron_repulsion_xps(skbuffer, 20, 0, 1, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpp(skbuffer, 32, 1, 4, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 68, 4, 10, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xds(skbuffer, 230, 20, 32, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xdp(skbuffer, 254, 32, 68, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xfs(skbuffer, 434, 230, 254, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xps(skbuffer, 665, 609, 610, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpp(skbuffer, 668, 610, 613, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 677, 613, 619, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 695, 619, 629, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 725, 629, 644, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xps(skbuffer, 23, 0, 665, 668, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpp(skbuffer, 41, 1, 668, 677, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpd(skbuffer, 86, 4, 677, 695, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpf(skbuffer, 140, 10, 695, 725, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xds(skbuffer, 236, 20, 23, 41, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xdp(skbuffer, 272, 32, 41, 86, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xdd(skbuffer, 326, 68, 86, 140, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xfs(skbuffer, 444, 230, 236, 272, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xfp(skbuffer, 474, 254, 272, 326, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xgs(skbuffer, 564, 434, 444, 474, cfactors, 6, 0);

            t3cfunc::ket_transform<4, 0>(sbuffer, 0, skbuffer, 564, 0);

            t3cfunc::ket_transform<4, 0>(sbuffer, 9, skbuffer, 579, 0);

            t3cfunc::ket_transform<4, 0>(sbuffer, 18, skbuffer, 594, 0);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 0, 4, 0, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecSGS_hpp */
