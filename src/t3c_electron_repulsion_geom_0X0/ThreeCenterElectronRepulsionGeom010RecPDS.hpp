#ifndef ThreeCenterElectronRepulsionGeom010RecPDS_hpp
#define ThreeCenterElectronRepulsionGeom010RecPDS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionGeom010ContrRecXDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPS.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
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

/// @brief Computes d^(1)/dC^(1)(P|1/|r-r'||DS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_pds(T& distributor,
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

    CSimdArray<double> pbuffer(95, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(72, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(333, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(45, 1);

    // setup Boys fuction data

    const CBoysFunc<4> bf_table;

    CSimdArray<double> bf_data(6, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 5, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 5);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 4, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 7, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 10, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 13, 0, 1, 4, 7, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 19, 1, 2, 7, 10, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 25, 4, 7, 13, 19, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 35, 0, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 38, 0, 4, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 47, 4, 13, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 65, 13, 25, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 35, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 38, 9, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {35, 38});

                pbuffer.scale(pfactors, 0, 2.0, {38, 47});

                pbuffer.scale(pfactors, 0, 2.0, {47, 65});

                pbuffer.scale(pfactors, 0, 2.0, {65, 95});

                t2cfunc::reduce(cbuffer, 12, pbuffer, 35, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15, pbuffer, 38, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 24, pbuffer, 47, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 42, pbuffer, 65, 30, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<1>(skbuffer, 0, cbuffer, 0, 0, 0);

            t3cfunc::bra_transform<1>(skbuffer, 3, cbuffer, 3, 0, 1);

            t3cfunc::bra_transform<1>(skbuffer, 183, cbuffer, 12, 0, 0);

            t3cfunc::bra_transform<1>(skbuffer, 186, cbuffer, 15, 0, 1);

            t3cfunc::bra_transform<1>(skbuffer, 195, cbuffer, 24, 0, 2);

            t3cfunc::bra_transform<1>(skbuffer, 213, cbuffer, 42, 0, 3);

            t3ceri::comp_hrr_electron_repulsion_xps(skbuffer, 12, 0, 3, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xss(skbuffer, 243, 183, 186, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xsp(skbuffer, 252, 186, 195, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xsd(skbuffer, 279, 195, 213, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xps(skbuffer, 21, 0, 243, 252, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xpp(skbuffer, 48, 3, 252, 279, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xds(skbuffer, 129, 12, 21, 48, cfactors, 6, 1);

            t3cfunc::ket_transform<2, 0>(sbuffer, 0, skbuffer, 129, 1);

            t3cfunc::ket_transform<2, 0>(sbuffer, 15, skbuffer, 147, 1);

            t3cfunc::ket_transform<2, 0>(sbuffer, 30, skbuffer, 165, 1);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 1, 2, 0, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecPDS_hpp */
