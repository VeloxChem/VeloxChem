#ifndef ThreeCenterElectronRepulsionGeom010RecSGP_hpp
#define ThreeCenterElectronRepulsionGeom010RecSGP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFP.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXGP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPP.hpp"
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

/// @brief Computes d^(1)/dC^(1)(S|1/|r-r'||GP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_sgp(T& distributor,
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

    CSimdArray<double> pbuffer(210, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(117, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(1476, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(81, 1);

    // setup Boys fuction data

    const CBoysFunc<6> bf_table;

    CSimdArray<double> bf_data(8, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 7, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 7);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 0);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 6, pfactors, 16, bf_data, 6);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 7, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 10, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 13, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 16, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 19, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 22, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 25, 0, 1, 7, 10, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 31, 1, 2, 10, 13, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 37, 2, 3, 13, 16, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 43, 3, 4, 16, 19, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 49, 4, 5, 19, 22, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 55, 7, 10, 25, 31, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 65, 10, 13, 31, 37, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 75, 13, 16, 37, 43, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 85, 16, 19, 43, 49, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 95, 25, 31, 55, 65, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 110, 31, 37, 65, 75, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 125, 37, 43, 75, 85, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 140, 55, 65, 95, 110, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssh(pbuffer, 161, 65, 75, 110, 125, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 182, 95, 110, 140, 161, pfactors, 20, 23, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 19, pbuffer, 95, 15, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {7, 10});

                pbuffer.scale(pfactors, 0, 2.0, {25, 31});

                pbuffer.scale(pfactors, 0, 2.0, {55, 65});

                pbuffer.scale(pfactors, 0, 2.0, {95, 110});

                pbuffer.scale(pfactors, 0, 2.0, {140, 161});

                pbuffer.scale(pfactors, 0, 2.0, {182, 210});

                t2cfunc::reduce(cbuffer, 34, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 37, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 43, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 53, pbuffer, 95, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 68, pbuffer, 140, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 89, pbuffer, 182, 28, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<0>(skbuffer, 0, cbuffer, 0, 0, 1);

            t3cfunc::bra_transform<0>(skbuffer, 3, cbuffer, 3, 0, 2);

            t3cfunc::bra_transform<0>(skbuffer, 9, cbuffer, 9, 0, 3);

            t3cfunc::bra_transform<0>(skbuffer, 19, cbuffer, 19, 0, 4);

            t3cfunc::bra_transform<0>(skbuffer, 1228, cbuffer, 34, 0, 1);

            t3cfunc::bra_transform<0>(skbuffer, 1231, cbuffer, 37, 0, 2);

            t3cfunc::bra_transform<0>(skbuffer, 1237, cbuffer, 43, 0, 3);

            t3cfunc::bra_transform<0>(skbuffer, 1247, cbuffer, 53, 0, 4);

            t3cfunc::bra_transform<0>(skbuffer, 1262, cbuffer, 68, 0, 5);

            t3cfunc::bra_transform<0>(skbuffer, 1283, cbuffer, 89, 0, 6);

            t3ceri::comp_hrr_electron_repulsion_xpp(skbuffer, 34, 0, 3, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 70, 3, 9, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 142, 9, 19, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xdp(skbuffer, 397, 34, 70, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xdd(skbuffer, 469, 70, 142, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xfp(skbuffer, 793, 397, 469, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpp(skbuffer, 1311, 1228, 1231, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 1320, 1231, 1237, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 1338, 1237, 1247, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 1368, 1247, 1262, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 1413, 1262, 1283, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpp(skbuffer, 43, 0, 1311, 1320, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpd(skbuffer, 88, 3, 1320, 1338, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpf(skbuffer, 172, 9, 1338, 1368, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpg(skbuffer, 262, 19, 1368, 1413, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xdp(skbuffer, 415, 34, 43, 88, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xdd(skbuffer, 505, 70, 88, 172, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xdf(skbuffer, 613, 142, 172, 262, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xfp(skbuffer, 823, 397, 415, 505, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xfd(skbuffer, 913, 469, 505, 613, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xgp(skbuffer, 1093, 793, 823, 913, cfactors, 6, 0);

            t3cfunc::ket_transform<4, 1>(sbuffer, 0, skbuffer, 1093, 0);

            t3cfunc::ket_transform<4, 1>(sbuffer, 27, skbuffer, 1138, 0);

            t3cfunc::ket_transform<4, 1>(sbuffer, 54, skbuffer, 1183, 0);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 0, 4, 1, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecSGP_hpp */
