#ifndef ThreeCenterElectronRepulsionGeom010RecSGD_hpp
#define ThreeCenterElectronRepulsionGeom010RecSGD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXDF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXFD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPG.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPH.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPI.hpp"
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

/// @brief Computes d^(1)/dC^(1)(S|1/|r-r'||GD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_sgd(T& distributor,
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

    CSimdArray<double> pbuffer(330, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(168, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(2433, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(135, 1);

    // setup Boys fuction data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 8);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 0);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 6, pfactors, 16, bf_data, 6);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 7, pfactors, 16, bf_data, 7);

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

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 238, 115, 130, 175, 196, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssi(pbuffer, 266, 130, 145, 196, 217, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssk(pbuffer, 294, 175, 196, 238, 266, pfactors, 20, 23, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 65, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 115, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 31, pbuffer, 175, 21, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {29, 35});

                pbuffer.scale(pfactors, 0, 2.0, {65, 75});

                pbuffer.scale(pfactors, 0, 2.0, {115, 130});

                pbuffer.scale(pfactors, 0, 2.0, {175, 196});

                pbuffer.scale(pfactors, 0, 2.0, {238, 266});

                pbuffer.scale(pfactors, 0, 2.0, {294, 330});

                t2cfunc::reduce(cbuffer, 52, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 58, pbuffer, 65, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 68, pbuffer, 115, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 83, pbuffer, 175, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 104, pbuffer, 238, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 132, pbuffer, 294, 36, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<0>(skbuffer, 0, cbuffer, 0, 0, 2);

            t3cfunc::bra_transform<0>(skbuffer, 6, cbuffer, 6, 0, 3);

            t3cfunc::bra_transform<0>(skbuffer, 16, cbuffer, 16, 0, 4);

            t3cfunc::bra_transform<0>(skbuffer, 31, cbuffer, 31, 0, 5);

            t3cfunc::bra_transform<0>(skbuffer, 2077, cbuffer, 52, 0, 2);

            t3cfunc::bra_transform<0>(skbuffer, 2083, cbuffer, 58, 0, 3);

            t3cfunc::bra_transform<0>(skbuffer, 2093, cbuffer, 68, 0, 4);

            t3cfunc::bra_transform<0>(skbuffer, 2108, cbuffer, 83, 0, 5);

            t3cfunc::bra_transform<0>(skbuffer, 2129, cbuffer, 104, 0, 6);

            t3cfunc::bra_transform<0>(skbuffer, 2157, cbuffer, 132, 0, 7);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 52, 0, 6, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 124, 6, 16, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 244, 16, 31, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xdd(skbuffer, 613, 52, 124, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xdf(skbuffer, 757, 124, 244, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xfd(skbuffer, 1267, 613, 757, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpd(skbuffer, 2193, 2077, 2083, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 2211, 2083, 2093, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpg(skbuffer, 2241, 2093, 2108, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xph(skbuffer, 2286, 2108, 2129, cfactors, 6, 0);

            t3ceri::comp_hrr_electron_repulsion_xpi(skbuffer, 2349, 2129, 2157, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpd(skbuffer, 70, 0, 2193, 2211, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpf(skbuffer, 154, 6, 2211, 2241, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xpg(skbuffer, 289, 16, 2241, 2286, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xph(skbuffer, 424, 31, 2286, 2349, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xdd(skbuffer, 649, 52, 70, 154, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xdf(skbuffer, 817, 124, 154, 289, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xdg(skbuffer, 997, 244, 289, 424, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xfd(skbuffer, 1327, 613, 649, 817, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xff(skbuffer, 1507, 757, 817, 997, cfactors, 6, 0);

            t3ceri::comp_ket_geom010_electron_repulsion_xgd(skbuffer, 1807, 1267, 1327, 1507, cfactors, 6, 0);

            t3cfunc::ket_transform<4, 2>(sbuffer, 0, skbuffer, 1807, 0);

            t3cfunc::ket_transform<4, 2>(sbuffer, 45, skbuffer, 1897, 0);

            t3cfunc::ket_transform<4, 2>(sbuffer, 90, skbuffer, 1987, 0);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 0, 4, 2, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecSGD_hpp */
