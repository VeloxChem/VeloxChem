#ifndef ElectronRepulsionGeom1100RecFFPD_hpp
#define ElectronRepulsionGeom1100RecFFPD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
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
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSF.hpp"
#include "ElectronRepulsionPrimRecSKSP.hpp"
#include "ElectronRepulsionPrimRecSLSD.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FF|1/|r-r'||PD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_ffpd(T& distributor,
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

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(9140, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(5120, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(14184, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(88860, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(6615, 1);

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

        cfactors.replicate_points(c_coords, ket_range, 0, 1);

        cfactors.replicate_points(d_coords, ket_range, 3, 1);

        t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        ckbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // loop over basis function pairs on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            // zero integral buffers

            cbuffer.zero();

            ckbuffer.zero();

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 195, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 198, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 201, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 204, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 207, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 210, 7, 8, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 213, 2, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 222, 3, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 231, 4, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 240, 5, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 249, 6, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 258, 7, 30, 33, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 267, 8, 33, 36, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 276, 15, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 294, 18, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 312, 21, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 330, 24, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 348, 27, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 366, 30, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 384, 33, 81, 87, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 402, 36, 87, 93, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 420, 51, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 450, 57, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 480, 63, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 510, 69, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 540, 75, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 570, 81, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 600, 87, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 630, 93, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 660, 2, 3, 195, 198, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 666, 3, 4, 198, 201, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 672, 4, 5, 201, 204, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 678, 5, 6, 204, 207, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 684, 6, 7, 207, 210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 690, 15, 18, 195, 213, 222, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 708, 18, 21, 198, 222, 231, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 726, 21, 24, 201, 231, 240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 744, 24, 27, 204, 240, 249, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 762, 27, 30, 207, 249, 258, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 780, 30, 33, 210, 258, 267, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 798, 45, 51, 213, 276, 294, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 834, 51, 57, 222, 294, 312, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 870, 57, 63, 231, 312, 330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 906, 63, 69, 240, 330, 348, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 942, 69, 75, 249, 348, 366, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 978, 75, 81, 258, 366, 384, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1014, 81, 87, 267, 384, 402, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1050, 105, 115, 294, 420, 450, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1110, 115, 125, 312, 450, 480, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1170, 125, 135, 330, 480, 510, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1230, 135, 145, 348, 510, 540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1290, 145, 155, 366, 540, 570, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1350, 155, 165, 384, 570, 600, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1410, 165, 175, 402, 600, 630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1470, 195, 198, 660, 666, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1480, 198, 201, 666, 672, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1490, 201, 204, 672, 678, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1500, 204, 207, 678, 684, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1510, 213, 222, 660, 690, 708, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1540, 222, 231, 666, 708, 726, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1570, 231, 240, 672, 726, 744, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1600, 240, 249, 678, 744, 762, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1630, 249, 258, 684, 762, 780, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1660, 276, 294, 690, 798, 834, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1720, 294, 312, 708, 834, 870, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1780, 312, 330, 726, 870, 906, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1840, 330, 348, 744, 906, 942, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1900, 348, 366, 762, 942, 978, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1960, 366, 384, 780, 978, 1014, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2020, 420, 450, 834, 1050, 1110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2120, 450, 480, 870, 1110, 1170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2220, 480, 510, 906, 1170, 1230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2320, 510, 540, 942, 1230, 1290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2420, 540, 570, 978, 1290, 1350, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2520, 570, 600, 1014, 1350, 1410, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2620, 660, 666, 1470, 1480, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2635, 666, 672, 1480, 1490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2650, 672, 678, 1490, 1500, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2665, 690, 708, 1470, 1510, 1540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2710, 708, 726, 1480, 1540, 1570, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2755, 726, 744, 1490, 1570, 1600, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2800, 744, 762, 1500, 1600, 1630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2845, 798, 834, 1510, 1660, 1720, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2935, 834, 870, 1540, 1720, 1780, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3025, 870, 906, 1570, 1780, 1840, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3115, 906, 942, 1600, 1840, 1900, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3205, 942, 978, 1630, 1900, 1960, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3295, 1050, 1110, 1720, 2020, 2120, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3445, 1110, 1170, 1780, 2120, 2220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3595, 1170, 1230, 1840, 2220, 2320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3745, 1230, 1290, 1900, 2320, 2420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3895, 1290, 1350, 1960, 2420, 2520, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 4045, 1470, 1480, 2620, 2635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 4066, 1480, 1490, 2635, 2650, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 4087, 1510, 1540, 2620, 2665, 2710, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 4150, 1540, 1570, 2635, 2710, 2755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 4213, 1570, 1600, 2650, 2755, 2800, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 4276, 1660, 1720, 2665, 2845, 2935, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 4402, 1720, 1780, 2710, 2935, 3025, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 4528, 1780, 1840, 2755, 3025, 3115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 4654, 1840, 1900, 2800, 3115, 3205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4780, 2020, 2120, 2935, 3295, 3445, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4990, 2120, 2220, 3025, 3445, 3595, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5200, 2220, 2320, 3115, 3595, 3745, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5410, 2320, 2420, 3205, 3745, 3895, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 5620, 2620, 2635, 4045, 4066, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 5648, 2665, 2710, 4045, 4087, 4150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 5732, 2710, 2755, 4066, 4150, 4213, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 5816, 2845, 2935, 4087, 4276, 4402, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 5984, 2935, 3025, 4150, 4402, 4528, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 6152, 3025, 3115, 4213, 4528, 4654, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 6320, 3295, 3445, 4402, 4780, 4990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 6600, 3445, 3595, 4528, 4990, 5200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 6880, 3595, 3745, 4654, 5200, 5410, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 7160, 4087, 4150, 5620, 5648, 5732, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 7268, 4276, 4402, 5648, 5816, 5984, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 7484, 4402, 4528, 5732, 5984, 6152, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 7700, 4780, 4990, 5984, 6320, 6600, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 8060, 4990, 5200, 6152, 6600, 6880, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsd(pbuffer, 8420, 5816, 5984, 7160, 7268, 7484, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsf(pbuffer, 8690, 6320, 6600, 7484, 7700, 8060, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 798, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 1050, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 96, pbuffer, 1660, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 156, pbuffer, 2020, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 256, pbuffer, 2845, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 346, pbuffer, 3295, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {2845, 2935});

                pbuffer.scale(2.0 * b_exp, {3295, 3445});

                pbuffer.scale(2.0 * b_exp, {4276, 4402});

                pbuffer.scale(2.0 * b_exp, {4780, 4990});

                pbuffer.scale(2.0 * b_exp, {5816, 5984});

                pbuffer.scale(2.0 * b_exp, {6320, 6600});

                t2cfunc::reduce(cbuffer, 496, pbuffer, 2845, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 586, pbuffer, 3295, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 736, pbuffer, 4276, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 862, pbuffer, 4780, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1072, pbuffer, 5816, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1240, pbuffer, 6320, 280, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {798, 834});

                pbuffer.scale(2.0 * a_exp, {1050, 1110});

                pbuffer.scale(2.0 * a_exp, {1660, 1720});

                pbuffer.scale(2.0 * a_exp, {2020, 2120});

                pbuffer.scale(a_exp / b_exp, {2845, 2935});

                pbuffer.scale(a_exp / b_exp, {3295, 3445});

                pbuffer.scale(a_exp / b_exp, {4276, 4402});

                pbuffer.scale(a_exp / b_exp, {4780, 4990});

                pbuffer.scale(a_exp / b_exp, {5816, 5984});

                pbuffer.scale(a_exp / b_exp, {6320, 6600});

                t2cfunc::reduce(cbuffer, 1520, pbuffer, 798, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1556, pbuffer, 1050, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1616, pbuffer, 1660, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1676, pbuffer, 2020, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1776, pbuffer, 2845, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1866, pbuffer, 3295, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2016, pbuffer, 4276, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2142, pbuffer, 4780, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2352, pbuffer, 5816, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2520, pbuffer, 6320, 280, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {2845, 2935});

                pbuffer.scale(2.0 * b_exp, {3295, 3445});

                pbuffer.scale(2.0 * b_exp, {4276, 4402});

                pbuffer.scale(2.0 * b_exp, {4780, 4990});

                pbuffer.scale(2.0 * b_exp, {5816, 5984});

                pbuffer.scale(2.0 * b_exp, {6320, 6600});

                pbuffer.scale(4.0 * a_exp * b_exp, {7268, 7484});

                pbuffer.scale(4.0 * a_exp * b_exp, {7700, 8060});

                pbuffer.scale(4.0 * a_exp * b_exp, {8420, 8690});

                pbuffer.scale(4.0 * a_exp * b_exp, {8690, 9140});

                t2cfunc::reduce(cbuffer, 2800, pbuffer, 2845, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2890, pbuffer, 3295, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3040, pbuffer, 4276, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3166, pbuffer, 4780, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3376, pbuffer, 5816, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3544, pbuffer, 6320, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3824, pbuffer, 7268, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4040, pbuffer, 7700, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4400, pbuffer, 8420, 270, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4670, pbuffer, 8690, 450, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 36, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 108, cbuffer, 96, 156, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 828, cbuffer, 256, 346, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3042, cbuffer, 496, 586, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3312, cbuffer, 736, 862, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3690, cbuffer, 1072, 1240, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4194, cbuffer, 1520, 1556, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4302, cbuffer, 1616, 1676, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5022, cbuffer, 1776, 1866, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 6102, cbuffer, 2016, 2142, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 7614, cbuffer, 2352, 2520, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 11574, cbuffer, 2800, 2890, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 11844, cbuffer, 3040, 3166, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 12222, cbuffer, 3376, 3544, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 12726, cbuffer, 3824, 4040, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 13374, cbuffer, 4400, 4670, cfactors, 6, 0, 8);

            t4cfunc::ket_transform<1, 2>(skbuffer, 0, ckbuffer, 0, 0, 2);

            t4cfunc::ket_transform<1, 2>(skbuffer, 90, ckbuffer, 108, 0, 3);

            t4cfunc::ket_transform<1, 2>(skbuffer, 2040, ckbuffer, 828, 0, 4);

            t4cfunc::ket_transform<1, 2>(skbuffer, 77505, ckbuffer, 3042, 0, 4);

            t4cfunc::ket_transform<1, 2>(skbuffer, 77730, ckbuffer, 3312, 0, 5);

            t4cfunc::ket_transform<1, 2>(skbuffer, 78045, ckbuffer, 3690, 0, 6);

            t4cfunc::ket_transform<1, 2>(skbuffer, 78465, ckbuffer, 4194, 0, 2);

            t4cfunc::ket_transform<1, 2>(skbuffer, 78555, ckbuffer, 4302, 0, 3);

            t4cfunc::ket_transform<1, 2>(skbuffer, 79155, ckbuffer, 5022, 0, 4);

            t4cfunc::ket_transform<1, 2>(skbuffer, 80055, ckbuffer, 6102, 0, 5);

            t4cfunc::ket_transform<1, 2>(skbuffer, 81315, ckbuffer, 7614, 0, 6);

            t4cfunc::ket_transform<1, 2>(skbuffer, 86685, ckbuffer, 11574, 0, 4);

            t4cfunc::ket_transform<1, 2>(skbuffer, 86910, ckbuffer, 11844, 0, 5);

            t4cfunc::ket_transform<1, 2>(skbuffer, 87225, ckbuffer, 12222, 0, 6);

            t4cfunc::ket_transform<1, 2>(skbuffer, 87645, ckbuffer, 12726, 0, 7);

            t4cfunc::ket_transform<1, 2>(skbuffer, 88185, ckbuffer, 13374, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 12525, 90, 2040, r_ab, 1, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 84615, 78555, 79155, r_ab, 1, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 85065, 79155, 80055, r_ab, 1, 2);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 85740, 80055, 81315, r_ab, 1, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 14325, 90, 84615, 85065, r_ab, 1, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 21750, 2040, 85065, 85740, r_ab, 1, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 41055, 12525, 14325, 21750, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 240, 0, 77505, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 2265, 90, 77730, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 4965, 2040, 78045, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 12975, 90, 240, 2265, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pgxx(skbuffer, 19725, 2040, 2265, 4965, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dfxx(skbuffer, 38355, 12525, 12975, 19725, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 78705, 78465, 86685, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 79380, 78555, 86910, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 80370, 79155, 87225, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 81735, 80055, 87645, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_skxx(skbuffer, 82995, 81315, 88185, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 690, 78555, 78705, 79380, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 2940, 79155, 79380, 80370, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 5910, 80055, 80370, 81735, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sixx(skbuffer, 8745, 81315, 81735, 82995, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 15675, 240, 84615, 690, 2940, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 23775, 2265, 85065, 2940, 5910, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_phxx(skbuffer, 29850, 4965, 85740, 5910, 8745, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 43755, 12975, 14325, 15675, 23775, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dgxx(skbuffer, 51855, 19725, 21750, 23775, 29850, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ffxx(skbuffer, 64005, 38355, 41055, 43755, 51855, r_ab, 1, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 64005, 1, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 735, skbuffer, 65505, 1, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1470, skbuffer, 67005, 1, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2205, skbuffer, 68505, 1, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2940, skbuffer, 70005, 1, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 3675, skbuffer, 71505, 1, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 4410, skbuffer, 73005, 1, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 5145, skbuffer, 74505, 1, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 5880, skbuffer, 76005, 1, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 1, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFFPD_hpp */
