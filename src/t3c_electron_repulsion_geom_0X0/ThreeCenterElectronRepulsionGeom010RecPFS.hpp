#ifndef ThreeCenterElectronRepulsionGeom010RecPFS_hpp
#define ThreeCenterElectronRepulsionGeom010RecPFS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXDS.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionContrRecXPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXFS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010ContrRecXSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSG.hpp"
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

/// @brief Computes d^(1)/dC^(1)(P|1/|r-r'||FS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom010_pfs(T& distributor,
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

    CSimdArray<double> pbuffer(175, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(135, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(945, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(63, 1);

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

                t4cfunc::comp_distances_wp(pfactors, 26, 17, r_a);

                t3cfunc::comp_boys_args(bf_data, 6, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 6);

                t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 5, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 8, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 11, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 14, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 17, 0, 1, 5, 8, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 23, 1, 2, 8, 11, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 29, 2, 3, 11, 14, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 35, 5, 8, 17, 23, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 45, 8, 11, 23, 29, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 55, 17, 23, 35, 45, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 70, 0, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 73, 0, 5, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 82, 5, 17, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 100, 17, 35, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 130, 35, 55, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 70, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 73, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12, pbuffer, 82, 18, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {70, 73});

                pbuffer.scale(pfactors, 0, 2.0, {73, 82});

                pbuffer.scale(pfactors, 0, 2.0, {82, 100});

                pbuffer.scale(pfactors, 0, 2.0, {100, 130});

                pbuffer.scale(pfactors, 0, 2.0, {130, 175});

                t2cfunc::reduce(cbuffer, 30, pbuffer, 70, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 33, pbuffer, 73, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 42, pbuffer, 82, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 100, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 130, 45, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<1>(skbuffer, 0, cbuffer, 0, 0, 0);

            t3cfunc::bra_transform<1>(skbuffer, 3, cbuffer, 3, 0, 1);

            t3cfunc::bra_transform<1>(skbuffer, 12, cbuffer, 12, 0, 2);

            t3cfunc::bra_transform<1>(skbuffer, 660, cbuffer, 30, 0, 0);

            t3cfunc::bra_transform<1>(skbuffer, 663, cbuffer, 33, 0, 1);

            t3cfunc::bra_transform<1>(skbuffer, 672, cbuffer, 42, 0, 2);

            t3cfunc::bra_transform<1>(skbuffer, 690, cbuffer, 60, 0, 3);

            t3cfunc::bra_transform<1>(skbuffer, 720, cbuffer, 90, 0, 4);

            t3ceri::comp_hrr_electron_repulsion_xps(skbuffer, 30, 0, 3, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xpp(skbuffer, 66, 3, 12, cfactors, 6, 1);

            t3ceri::comp_hrr_electron_repulsion_xds(skbuffer, 336, 30, 66, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xss(skbuffer, 765, 660, 663, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xsp(skbuffer, 774, 663, 672, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xsd(skbuffer, 801, 672, 690, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xsf(skbuffer, 855, 690, 720, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xps(skbuffer, 39, 0, 765, 774, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xpp(skbuffer, 93, 3, 774, 801, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xpd(skbuffer, 174, 12, 801, 855, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xds(skbuffer, 354, 30, 39, 93, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xdp(skbuffer, 408, 66, 93, 174, cfactors, 6, 1);

            t3ceri::comp_ket_geom010_electron_repulsion_xfs(skbuffer, 570, 336, 354, 408, cfactors, 6, 1);

            t3cfunc::ket_transform<3, 0>(sbuffer, 0, skbuffer, 570, 1);

            t3cfunc::ket_transform<3, 0>(sbuffer, 21, skbuffer, 600, 1);

            t3cfunc::ket_transform<3, 0>(sbuffer, 42, skbuffer, 630, 1);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 1, 3, 0, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom010RecPFS_hpp */
