#ifndef ElectronRepulsionGeom1100RecFFSF_hpp
#define ElectronRepulsionGeom1100RecFFSF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSIXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSKXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecFFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPHXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSIXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSF.hpp"
#include "ElectronRepulsionPrimRecSLSF.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FF|1/|r-r'||SF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_ffsf(T& distributor,
                                      const CGtoPairBlock& bra_gto_pair_block,
                                      const CGtoPairBlock& ket_gto_pair_block,
                                      const std::pair<size_t, size_t>& bra_indices,
                                      const std::pair<size_t, size_t>& ket_indices) -> void
{
    // intialize GTOs pair data on bra side

    const auto a_coords = bra_gto_pair_block.bra_coordinates();

    const auto b_coords = bra_gto_pair_block.ket_coordinates();

    const auto a_vec_exps = bra_gto_pair_block.bra_exponents();

    const auto b_vec_exps = bra_gto_pair_block.ket_exponents();

    const auto ab_vec_norms = bra_gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = bra_gto_pair_block.overlap_factors();

    const auto a_indices = bra_gto_pair_block.bra_orbital_indices();

    const auto b_indices = bra_gto_pair_block.ket_orbital_indices();

    const auto bra_ncgtos = bra_gto_pair_block.number_of_contracted_pairs();

    const auto bra_npgtos = bra_gto_pair_block.number_of_primitive_pairs();

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

    CSimdArray<double> pbuffer(7716, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(3200, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(41468, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3087, 1);

    // setup Boys fuction data

    const CBoysFunc<11> bf_table;

    CSimdArray<double> bf_data(13, ket_npgtos);

    // set up range seperation factor

    const auto use_rs = distributor.need_omega();

    const auto omega = distributor.get_omega();

    // set up ket partitioning

    const auto ket_dim = ket_indices.second - ket_indices.first;

    const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());

    for (size_t i = 0; i < ket_blocks; i++)
    {
        auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);

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

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            // zero integral buffers

            cbuffer.zero();

            skbuffer.zero();

            sbuffer.zero();

            // set up coordinates on bra side

            const auto r_a = a_coords[j];

            const auto r_b = b_coords[j];

            const auto a_xyz = r_a.coordinates();

            const auto b_xyz = r_b.coordinates();

            const auto r_ab = TPoint<double>({a_xyz[0] - b_xyz[0], a_xyz[1] - b_xyz[1], a_xyz[2] - b_xyz[2]});

            for (int k = 0; k < bra_npgtos; k++)
            {
                const auto a_exp = a_vec_exps[k * bra_ncgtos + j];

                const auto b_exp = b_vec_exps[k * bra_ncgtos + j];

                const auto ab_norm = ab_vec_norms[k * bra_ncgtos + j];

                const auto ab_ovl = ab_vec_ovls[k * bra_ncgtos + j];

                const auto p_x = (a_xyz[0] * a_exp + b_xyz[0] * b_exp) / (a_exp + b_exp);

                const auto p_y = (a_xyz[1] * a_exp + b_xyz[1] * b_exp) / (a_exp + b_exp);

                const auto p_z = (a_xyz[2] * a_exp + b_xyz[2] * b_exp) / (a_exp + b_exp);

                const auto r_p = TPoint<double>({p_x, p_y, p_z});

                const auto pb_x = p_x - b_xyz[0];

                const auto pb_y = p_y - b_xyz[1];

                const auto pb_z = p_z - b_xyz[2];

                const auto r_pb = TPoint<double>({pb_x, pb_y, pb_z});

                t4cfunc::comp_coordinates_q(pfactors, 10, 4, 7);

                t4cfunc::comp_distances_pq(pfactors, 13, 10, r_p);

                t4cfunc::comp_coordinates_w(pfactors, 17, 10, r_p, a_exp, b_exp);

                t4cfunc::comp_distances_qd(pfactors, 20, 10, 7);

                t4cfunc::comp_distances_wq(pfactors, 23, 17, 10);

                t4cfunc::comp_distances_wp(pfactors, 26, 17, r_p);

                if (use_rs)
                {
                    t4cfunc::comp_boys_args(bf_data, 12, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 12, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 12, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 12);
                }

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 6, pfactors, 16, bf_data, 6);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 7, pfactors, 16, bf_data, 7);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 8, pfactors, 16, bf_data, 8);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 9, pfactors, 16, bf_data, 9);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 10, pfactors, 16, bf_data, 10);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 11, pfactors, 16, bf_data, 11);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 12, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 21, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 24, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 27, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 30, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 33, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 36, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 39, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 42, 10, 11, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 45, 0, 1, 12, 15, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 51, 1, 2, 15, 18, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 57, 2, 3, 18, 21, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 63, 3, 4, 21, 24, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 69, 4, 5, 24, 27, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 75, 5, 6, 27, 30, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 81, 6, 7, 30, 33, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 87, 7, 8, 33, 36, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 93, 8, 9, 36, 39, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 99, 9, 10, 39, 42, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 12, 15, 45, 51, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 15, 18, 51, 57, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 18, 21, 57, 63, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 21, 24, 63, 69, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 24, 27, 69, 75, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 27, 30, 75, 81, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 30, 33, 81, 87, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 175, 33, 36, 87, 93, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 185, 36, 39, 93, 99, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 195, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 198, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 201, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 204, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 207, 7, 8, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 210, 3, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 219, 4, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 228, 5, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 237, 6, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 246, 7, 30, 33, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 255, 8, 33, 36, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 264, 18, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 282, 21, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 300, 24, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 318, 27, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 336, 30, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 354, 33, 81, 87, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 372, 36, 87, 93, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 390, 51, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 420, 57, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 450, 63, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 480, 69, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 510, 75, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 540, 81, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 570, 87, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 600, 93, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 630, 3, 4, 195, 198, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 636, 4, 5, 198, 201, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 642, 5, 6, 201, 204, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 648, 6, 7, 204, 207, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 654, 18, 21, 195, 210, 219, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 672, 21, 24, 198, 219, 228, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 690, 24, 27, 201, 228, 237, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 708, 27, 30, 204, 237, 246, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 726, 30, 33, 207, 246, 255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 744, 51, 57, 210, 264, 282, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 780, 57, 63, 219, 282, 300, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 816, 63, 69, 228, 300, 318, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 852, 69, 75, 237, 318, 336, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 888, 75, 81, 246, 336, 354, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 924, 81, 87, 255, 354, 372, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 960, 105, 115, 264, 390, 420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1020, 115, 125, 282, 420, 450, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1080, 125, 135, 300, 450, 480, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1140, 135, 145, 318, 480, 510, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1200, 145, 155, 336, 510, 540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1260, 155, 165, 354, 540, 570, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1320, 165, 175, 372, 570, 600, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1380, 195, 198, 630, 636, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1390, 198, 201, 636, 642, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1400, 201, 204, 642, 648, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1410, 210, 219, 630, 654, 672, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1440, 219, 228, 636, 672, 690, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1470, 228, 237, 642, 690, 708, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1500, 237, 246, 648, 708, 726, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1530, 264, 282, 654, 744, 780, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1590, 282, 300, 672, 780, 816, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1650, 300, 318, 690, 816, 852, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1710, 318, 336, 708, 852, 888, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1770, 336, 354, 726, 888, 924, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1830, 390, 420, 744, 960, 1020, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1930, 420, 450, 780, 1020, 1080, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2030, 450, 480, 816, 1080, 1140, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2130, 480, 510, 852, 1140, 1200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2230, 510, 540, 888, 1200, 1260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2330, 540, 570, 924, 1260, 1320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2430, 630, 636, 1380, 1390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2445, 636, 642, 1390, 1400, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2460, 654, 672, 1380, 1410, 1440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2505, 672, 690, 1390, 1440, 1470, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2550, 690, 708, 1400, 1470, 1500, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2595, 744, 780, 1410, 1530, 1590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2685, 780, 816, 1440, 1590, 1650, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2775, 816, 852, 1470, 1650, 1710, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2865, 852, 888, 1500, 1710, 1770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2955, 960, 1020, 1530, 1830, 1930, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3105, 1020, 1080, 1590, 1930, 2030, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3255, 1080, 1140, 1650, 2030, 2130, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3405, 1140, 1200, 1710, 2130, 2230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3555, 1200, 1260, 1770, 2230, 2330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 3705, 1380, 1390, 2430, 2445, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3726, 1410, 1440, 2430, 2460, 2505, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3789, 1440, 1470, 2445, 2505, 2550, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3852, 1530, 1590, 2460, 2595, 2685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3978, 1590, 1650, 2505, 2685, 2775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 4104, 1650, 1710, 2550, 2775, 2865, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4230, 1830, 1930, 2595, 2955, 3105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4440, 1930, 2030, 2685, 3105, 3255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4650, 2030, 2130, 2775, 3255, 3405, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4860, 2130, 2230, 2865, 3405, 3555, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 5070, 2460, 2505, 3705, 3726, 3789, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 5154, 2595, 2685, 3726, 3852, 3978, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 5322, 2685, 2775, 3789, 3978, 4104, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 5490, 2955, 3105, 3852, 4230, 4440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 5770, 3105, 3255, 3978, 4440, 4650, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 6050, 3255, 3405, 4104, 4650, 4860, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 6330, 3852, 3978, 5070, 5154, 5322, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 6546, 4230, 4440, 5154, 5490, 5770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 6906, 4440, 4650, 5322, 5770, 6050, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsf(pbuffer, 7266, 5490, 5770, 6330, 6546, 6906, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 960, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 1830, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 2955, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {2955, 3105});

                pbuffer.scale(2.0 * b_exp, {4230, 4440});

                pbuffer.scale(2.0 * b_exp, {5490, 5770});

                t2cfunc::reduce(cbuffer, 310, pbuffer, 2955, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 460, pbuffer, 4230, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 670, pbuffer, 5490, 280, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {960, 1020});

                pbuffer.scale(2.0 * a_exp, {1830, 1930});

                pbuffer.scale(a_exp / b_exp, {2955, 3105});

                pbuffer.scale(a_exp / b_exp, {4230, 4440});

                pbuffer.scale(a_exp / b_exp, {5490, 5770});

                t2cfunc::reduce(cbuffer, 950, pbuffer, 960, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1010, pbuffer, 1830, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1110, pbuffer, 2955, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1260, pbuffer, 4230, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1470, pbuffer, 5490, 280, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {2955, 3105});

                pbuffer.scale(2.0 * b_exp, {4230, 4440});

                pbuffer.scale(2.0 * b_exp, {5490, 5770});

                pbuffer.scale(4.0 * a_exp * b_exp, {6546, 6906});

                pbuffer.scale(4.0 * a_exp * b_exp, {7266, 7716});

                t2cfunc::reduce(cbuffer, 1750, pbuffer, 2955, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1900, pbuffer, 4230, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2110, pbuffer, 5490, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2390, pbuffer, 6546, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2750, pbuffer, 7266, 450, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 3>(skbuffer, 0, cbuffer, 0, 0, 2);

            t4cfunc::ket_transform<0, 3>(skbuffer, 42, cbuffer, 60, 0, 3);

            t4cfunc::ket_transform<0, 3>(skbuffer, 952, cbuffer, 160, 0, 4);

            t4cfunc::ket_transform<0, 3>(skbuffer, 36169, cbuffer, 310, 0, 4);

            t4cfunc::ket_transform<0, 3>(skbuffer, 36274, cbuffer, 460, 0, 5);

            t4cfunc::ket_transform<0, 3>(skbuffer, 36421, cbuffer, 670, 0, 6);

            t4cfunc::ket_transform<0, 3>(skbuffer, 36617, cbuffer, 950, 0, 2);

            t4cfunc::ket_transform<0, 3>(skbuffer, 36659, cbuffer, 1010, 0, 3);

            t4cfunc::ket_transform<0, 3>(skbuffer, 36939, cbuffer, 1110, 0, 4);

            t4cfunc::ket_transform<0, 3>(skbuffer, 37359, cbuffer, 1260, 0, 5);

            t4cfunc::ket_transform<0, 3>(skbuffer, 37947, cbuffer, 1470, 0, 6);

            t4cfunc::ket_transform<0, 3>(skbuffer, 40453, cbuffer, 1750, 0, 4);

            t4cfunc::ket_transform<0, 3>(skbuffer, 40558, cbuffer, 1900, 0, 5);

            t4cfunc::ket_transform<0, 3>(skbuffer, 40705, cbuffer, 2110, 0, 6);

            t4cfunc::ket_transform<0, 3>(skbuffer, 40901, cbuffer, 2390, 0, 7);

            t4cfunc::ket_transform<0, 3>(skbuffer, 41153, cbuffer, 2750, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 5845, 42, 952, r_ab, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 39487, 36659, 36939, r_ab, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 39697, 36939, 37359, r_ab, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 40012, 37359, 37947, r_ab, 0, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 6685, 42, 39487, 39697, r_ab, 0, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 10150, 952, 39697, 40012, r_ab, 0, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 19159, 5845, 6685, 10150, r_ab, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 112, 0, 36169, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 1057, 42, 36274, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 2317, 952, 36421, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 6055, 42, 112, 1057, r_ab, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pgxx(skbuffer, 9205, 952, 1057, 2317, r_ab, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dfxx(skbuffer, 17899, 5845, 6055, 9205, r_ab, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 36729, 36617, 40453, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 37044, 36659, 40558, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 37506, 36939, 40705, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 38143, 37359, 40901, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_skxx(skbuffer, 38731, 37947, 41153, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 322, 36659, 36729, 37044, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 1372, 36939, 37044, 37506, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 2758, 37359, 37506, 38143, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sixx(skbuffer, 4081, 37947, 38143, 38731, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 7315, 112, 39487, 322, 1372, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 11095, 1057, 39697, 1372, 2758, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_phxx(skbuffer, 13930, 2317, 40012, 2758, 4081, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 20419, 6055, 6685, 7315, 11095, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dgxx(skbuffer, 24199, 9205, 10150, 11095, 13930, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ffxx(skbuffer, 29869, 17899, 19159, 20419, 24199, r_ab, 0, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 29869, 0, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 343, skbuffer, 30569, 0, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 686, skbuffer, 31269, 0, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1029, skbuffer, 31969, 0, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1372, skbuffer, 32669, 0, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1715, skbuffer, 33369, 0, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2058, skbuffer, 34069, 0, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2401, skbuffer, 34769, 0, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2744, skbuffer, 35469, 0, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 0, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFFSF_hpp */
