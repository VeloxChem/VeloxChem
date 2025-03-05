#ifndef ThreeCenterElectronRepulsionRecIPF_hpp
#define ThreeCenterElectronRepulsionRecIPF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ThreeCenterElectronRepulsionContrRecXPF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecDSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecFSS.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecGSP.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSD.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecHSG.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISF.hpp"
#include "ThreeCenterElectronRepulsionPrimRecISG.hpp"
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
#include "GtoPairBlock.hpp"
#include "GtoBlock.hpp"
#include "BatchFunc.hpp"

namespace t3ceri { // t3ceri namespace

/// @brief Computes (I|1/|r-r'||PF)  integrals for basis functions block and basis function pairs block.
/// @param distributor The pointer to integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.
template <class T>
inline auto
comp_electron_repulsion_ipf(T& distributor,
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

    CSimdArray<double> pbuffer(5280, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(700, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(715, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(273, 1);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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

                t3cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp);

                bf_table.compute(bf_data, 0, 11);

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

                t3ceri::comp_prim_electron_repulsion_sss(pbuffer, 10, pfactors, 16, bf_data, 10);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 11, 0, 1, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 14, 1, 2, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 17, 2, 3, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 20, 3, 4, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 23, 4, 5, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 26, 5, 6, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 29, 6, 7, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 32, 7, 8, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 35, 8, 9, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssp(pbuffer, 38, 9, 10, pfactors, 20, 23);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 41, 0, 1, 11, 14, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 47, 1, 2, 14, 17, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 53, 2, 3, 17, 20, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 59, 3, 4, 20, 23, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 65, 4, 5, 23, 26, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 71, 5, 6, 26, 29, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 77, 6, 7, 29, 32, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 83, 7, 8, 32, 35, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssd(pbuffer, 89, 8, 9, 35, 38, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 95, 11, 14, 41, 47, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 105, 14, 17, 47, 53, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 115, 17, 20, 53, 59, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 125, 20, 23, 59, 65, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 135, 23, 26, 65, 71, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 145, 26, 29, 71, 77, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 155, 29, 32, 77, 83, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssf(pbuffer, 165, 32, 35, 83, 89, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 175, 41, 47, 95, 105, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 190, 47, 53, 105, 115, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 205, 53, 59, 115, 125, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 220, 59, 65, 125, 135, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 235, 65, 71, 135, 145, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 250, 71, 77, 145, 155, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_ssg(pbuffer, 265, 77, 83, 155, 165, pfactors, 20, 23, a_exp);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 280, 4, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 283, 5, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_pss(pbuffer, 286, 6, pfactors, 26);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 289, 4, 23, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 298, 5, 26, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psp(pbuffer, 307, 6, 29, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 316, 17, 53, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 334, 20, 59, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 352, 23, 65, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 370, 26, 71, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psd(pbuffer, 388, 29, 77, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 406, 53, 115, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 436, 59, 125, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 466, 65, 135, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 496, 71, 145, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psf(pbuffer, 526, 77, 155, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 556, 115, 205, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 601, 125, 220, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 646, 135, 235, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 691, 145, 250, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_psg(pbuffer, 736, 155, 265, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dss(pbuffer, 781, 4, 5, 286, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 787, 17, 20, 280, 289, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 805, 20, 23, 283, 298, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsp(pbuffer, 823, 23, 26, 286, 307, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 841, 53, 59, 289, 352, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 877, 59, 65, 298, 370, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsd(pbuffer, 913, 65, 71, 307, 388, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 949, 95, 105, 316, 406, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 1009, 105, 115, 334, 436, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 1069, 115, 125, 352, 466, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 1129, 125, 135, 370, 496, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsf(pbuffer, 1189, 135, 145, 388, 526, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1249, 175, 190, 406, 556, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1339, 190, 205, 436, 601, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1429, 205, 220, 466, 646, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1519, 220, 235, 496, 691, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_dsg(pbuffer, 1609, 235, 250, 526, 736, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fss(pbuffer, 1699, 280, 283, 781, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsp(pbuffer, 1709, 289, 298, 781, 823, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 1739, 316, 334, 787, 841, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 1799, 334, 352, 805, 877, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsd(pbuffer, 1859, 352, 370, 823, 913, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 1919, 406, 436, 841, 1069, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 2019, 436, 466, 877, 1129, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsf(pbuffer, 2119, 466, 496, 913, 1189, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 2219, 556, 601, 1069, 1429, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 2369, 601, 646, 1129, 1519, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_fsg(pbuffer, 2519, 646, 691, 1189, 1609, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsp(pbuffer, 2669, 787, 805, 1699, 1709, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsd(pbuffer, 2714, 841, 877, 1709, 1859, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 2804, 949, 1009, 1739, 1919, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 2954, 1009, 1069, 1799, 2019, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsf(pbuffer, 3104, 1069, 1129, 1859, 2119, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 3254, 1249, 1339, 1919, 2219, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 3479, 1339, 1429, 2019, 2369, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_gsg(pbuffer, 3704, 1429, 1519, 2119, 2519, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsd(pbuffer, 3929, 1739, 1799, 2669, 2714, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsf(pbuffer, 4055, 1919, 2019, 2714, 3104, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_hsg(pbuffer, 4265, 2219, 2369, 3104, 3704, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isf(pbuffer, 4580, 2804, 2954, 3929, 4055, pfactors, 26, a_exp);

                t3ceri::comp_prim_electron_repulsion_isg(pbuffer, 4860, 3254, 3479, 4055, 4265, pfactors, 26, a_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 4580, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 280, pbuffer, 4860, 420, ket_width, ket_npgtos);

            }

            t3cfunc::bra_transform<6>(skbuffer, 0, cbuffer, 0, 0, 3);

            t3cfunc::bra_transform<6>(skbuffer, 130, cbuffer, 280, 0, 4);

            t3ceri::comp_hrr_electron_repulsion_xpf(skbuffer, 325, 0, 130, cfactors, 6, 6);

            t3cfunc::ket_transform<1, 3>(sbuffer, 0, skbuffer, 325, 6);

            distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, 6, 1, 3, j, ket_range);
        }
    }

}

} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionRecIPF_hpp */
