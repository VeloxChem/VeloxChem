#ifndef ThreeCenterElectronRepulsionRecISP_hpp
#define ThreeCenterElectronRepulsionRecISP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T3CUtils.hpp"
#include "T2CUtils.hpp"
#include "GtoPairBlock.hpp"
#include "GtoBlock.hpp"
#include "BatchFunc.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes (I|1/|r-r'||SP)  integrals for basis functions block and basis function pairs block.
/// @param distributor The pointer to integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_isp(T& distributor,
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

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(635, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(84, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(39, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(39, 1);

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

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 29, 2, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 32, 3, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 35, 4, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 38, 5, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 41, 6, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 44, 2, 14, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 53, 3, 17, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 62, 4, 20, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 71, 5, 23, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 80, 6, 26, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 89, 2, 3, 35, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 95, 3, 4, 38, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 101, 4, 5, 41, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 107, 8, 11, 29, 44, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 125, 11, 14, 32, 53, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 143, 14, 17, 35, 62, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 161, 17, 20, 38, 71, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 179, 20, 23, 41, 80, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 197, 29, 32, 89, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 207, 32, 35, 95, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 217, 35, 38, 101, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 227, 44, 53, 89, 143, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 257, 53, 62, 95, 161, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 287, 62, 71, 101, 179, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gss(pbuffer, 317, 89, 95, 217, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsp(pbuffer, 332, 107, 125, 197, 227, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsp(pbuffer, 377, 125, 143, 207, 257, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsp(pbuffer, 422, 143, 161, 217, 287, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hss(pbuffer, 467, 197, 207, 317, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsp(pbuffer, 488, 227, 257, 317, 422, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isp(pbuffer, 551, 332, 377, 467, 488, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 551, 84, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<6>(skbuffer, 0, cbuffer, 0, 0, 1);

            t3cfunc::ket_transform<0, 1>(sbuffer, 0, skbuffer, 0, 6);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 6, 0, 1, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionRecISP_hpp */
