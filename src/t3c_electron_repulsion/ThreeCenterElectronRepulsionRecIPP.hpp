#ifndef ThreeCenterElectronRepulsionRecIPP_hpp
#define ThreeCenterElectronRepulsionRecIPP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXPP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecSSD.hpp"
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

/// @brief Computes (I|1/|r-r'||PP)  integrals for basis functions block and basis function pairs block.
/// @param distributor The pointer to integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_ipp(T& distributor,
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

    CSimdArray<double> pbuffer(1695, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(252, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(234, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(117, 1);

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

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 0, pfactors, 16, bf_data, 0);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 1, pfactors, 16, bf_data, 1);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 2, pfactors, 16, bf_data, 2);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 3, pfactors, 16, bf_data, 3);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 4, pfactors, 16, bf_data, 4);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 5, pfactors, 16, bf_data, 5);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 6, pfactors, 16, bf_data, 6);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 7, pfactors, 16, bf_data, 7);

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 8, pfactors, 16, bf_data, 8);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 9, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 12, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 15, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 18, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 21, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 24, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 27, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 30, 7, 8, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 33, 0, 1, 9, 12, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 39, 1, 2, 12, 15, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 45, 2, 3, 15, 18, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 51, 3, 4, 18, 21, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 57, 4, 5, 21, 24, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 63, 5, 6, 24, 27, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 69, 6, 7, 27, 30, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 75, 2, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 78, 3, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 81, 4, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 84, 5, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 87, 6, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 90, 2, 15, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 99, 3, 18, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 108, 4, 21, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 117, 5, 24, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 126, 6, 27, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 135, 15, 45, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 153, 18, 51, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 171, 21, 57, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 189, 24, 63, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 207, 27, 69, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 225, 2, 3, 81, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 231, 3, 4, 84, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 237, 4, 5, 87, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 243, 9, 12, 75, 90, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 261, 12, 15, 78, 99, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 279, 15, 18, 81, 108, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 297, 18, 21, 84, 117, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 315, 21, 24, 87, 126, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 333, 33, 39, 90, 135, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 369, 39, 45, 99, 153, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 405, 45, 51, 108, 171, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 441, 51, 57, 117, 189, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 477, 57, 63, 126, 207, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 513, 75, 78, 225, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 523, 78, 81, 231, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 533, 81, 84, 237, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 543, 90, 99, 225, 279, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 573, 99, 108, 231, 297, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 603, 108, 117, 237, 315, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 633, 135, 153, 279, 405, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 693, 153, 171, 297, 441, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 753, 171, 189, 315, 477, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gss(pbuffer, 813, 225, 231, 533, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsp(pbuffer, 828, 243, 261, 513, 543, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsp(pbuffer, 873, 261, 279, 523, 573, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsp(pbuffer, 918, 279, 297, 533, 603, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 963, 333, 369, 543, 633, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 1053, 369, 405, 573, 693, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 1143, 405, 441, 603, 753, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hss(pbuffer, 1233, 513, 523, 813, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsp(pbuffer, 1254, 543, 573, 813, 918, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsd(pbuffer, 1317, 633, 693, 918, 1143, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isp(pbuffer, 1443, 828, 873, 1233, 1254, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isd(pbuffer, 1527, 963, 1053, 1254, 1317, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1443, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 84, pbuffer, 1527, 168, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<6>(skbuffer, 0, cbuffer, 0, 0, 1);

            t3cfunc::bra_transform<6>(skbuffer, 39, cbuffer, 84, 0, 2);

            t3ceri::comp_hrr_electron_repulsion_xpp(skbuffer, 117, 0, 39, cfactors, 6, 6);

            t3cfunc::ket_transform<1, 1>(sbuffer, 0, skbuffer, 117, 6);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 6, 1, 1, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionRecIPP_hpp */
