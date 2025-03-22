#ifndef ThreeCenterElectronRepulsionGeom100RecHSS_hpp
#define ThreeCenterElectronRepulsionGeom100RecHSS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionGeom100ContrRecHXX.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)(H|1/|r-r'||SS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_geom100_hss(T& distributor,
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

    CSimdArray<double> pfactors(23, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(176, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(134, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(33, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(33, 1);

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

                t4cfunc::comp_distances_wp(pfactors, 20, 17, r_a);

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

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 7, 2, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 10, 3, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 13, 4, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 16, 5, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 19, 6, pfactors, 20);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 22, 0, 1, 7, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 28, 1, 2, 10, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 34, 2, 3, 13, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 40, 3, 4, 16, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 46, 4, 5, 19, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 52, 7, 10, 34, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 62, 10, 13, 40, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 72, 13, 16, 46, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_gss(pbuffer, 82, 22, 28, 52, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_gss(pbuffer, 97, 28, 34, 62, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_gss(pbuffer, 112, 34, 40, 72, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_hss(pbuffer, 127, 52, 62, 112, pfactors, 20, a_exp);

                t3ceri::comp_prim_electron_repulsion_iss(pbuffer, 148, 82, 97, 127, pfactors, 20, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 82, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 78, pbuffer, 148, 28, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {148, 176});

                t2cfunc::reduce(cbuffer, 106, pbuffer, 148, 28, ket_width, ket_npgtos);

            }

            t3ceri::comp_bra_geom1_electron_repulsion_hxx(cbuffer, 15, 0, 106, 0, 0);

            t3cfunc::bra_transform<5>(skbuffer, 0, cbuffer, 15, 0, 0);

            t3cfunc::bra_transform<5>(skbuffer, 11, cbuffer, 36, 0, 0);

            t3cfunc::bra_transform<5>(skbuffer, 22, cbuffer, 57, 0, 0);

            t3cfunc::ket_transform<0, 0>(sbuffer, 0, skbuffer, 0, 5);

            t3cfunc::ket_transform<0, 0>(sbuffer, 11, skbuffer, 11, 5);

            t3cfunc::ket_transform<0, 0>(sbuffer, 22, skbuffer, 22, 5);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 5, 0, 0, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionGeom100RecHSS_hpp */
