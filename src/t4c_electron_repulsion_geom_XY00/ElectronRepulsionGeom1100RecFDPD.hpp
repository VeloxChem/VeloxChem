#ifndef ElectronRepulsionGeom1100RecFDPD_hpp
#define ElectronRepulsionGeom1100RecFDPD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0100ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSIXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecFDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSHXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FD|1/|r-r'||PD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_fdpd(T& distributor,
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

    CSimdArray<double> pbuffer(6056, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(3680, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(10134, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(58305, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(4725, 1);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 11, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 11);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 29, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 32, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 35, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 38, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 41, 0, 1, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 47, 1, 2, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 2, 3, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 3, 4, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 65, 4, 5, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 71, 5, 6, 26, 29, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 77, 6, 7, 29, 32, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 83, 7, 8, 32, 35, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 89, 8, 9, 35, 38, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 11, 14, 41, 47, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 14, 17, 47, 53, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 17, 20, 53, 59, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 20, 23, 59, 65, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 23, 26, 65, 71, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 26, 29, 71, 77, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 29, 32, 77, 83, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 32, 35, 83, 89, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 175, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 178, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 181, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 184, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 187, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 190, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 199, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 208, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 217, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 226, 6, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 235, 7, 29, 32, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 244, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 262, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 280, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 298, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 316, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 334, 29, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 352, 32, 77, 83, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 370, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 400, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 430, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 460, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 490, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 520, 77, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 550, 83, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 580, 2, 3, 175, 178, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 586, 3, 4, 178, 181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 592, 4, 5, 181, 184, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 598, 5, 6, 184, 187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 604, 14, 17, 175, 190, 199, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 622, 17, 20, 178, 199, 208, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 640, 20, 23, 181, 208, 217, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 658, 23, 26, 184, 217, 226, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 676, 26, 29, 187, 226, 235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 694, 41, 47, 190, 244, 262, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 730, 47, 53, 199, 262, 280, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 766, 53, 59, 208, 280, 298, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 802, 59, 65, 217, 298, 316, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 838, 65, 71, 226, 316, 334, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 874, 71, 77, 235, 334, 352, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 910, 95, 105, 262, 370, 400, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 970, 105, 115, 280, 400, 430, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1030, 115, 125, 298, 430, 460, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1090, 125, 135, 316, 460, 490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1150, 135, 145, 334, 490, 520, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1210, 145, 155, 352, 520, 550, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1270, 175, 178, 580, 586, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1280, 178, 181, 586, 592, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1290, 181, 184, 592, 598, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1300, 190, 199, 580, 604, 622, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1330, 199, 208, 586, 622, 640, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1360, 208, 217, 592, 640, 658, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1390, 217, 226, 598, 658, 676, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1420, 244, 262, 604, 694, 730, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1480, 262, 280, 622, 730, 766, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1540, 280, 298, 640, 766, 802, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1600, 298, 316, 658, 802, 838, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1660, 316, 334, 676, 838, 874, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1720, 370, 400, 730, 910, 970, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1820, 400, 430, 766, 970, 1030, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1920, 430, 460, 802, 1030, 1090, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2020, 460, 490, 838, 1090, 1150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2120, 490, 520, 874, 1150, 1210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2220, 580, 586, 1270, 1280, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2235, 586, 592, 1280, 1290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2250, 604, 622, 1270, 1300, 1330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2295, 622, 640, 1280, 1330, 1360, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2340, 640, 658, 1290, 1360, 1390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2385, 694, 730, 1300, 1420, 1480, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2475, 730, 766, 1330, 1480, 1540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2565, 766, 802, 1360, 1540, 1600, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2655, 802, 838, 1390, 1600, 1660, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2745, 910, 970, 1480, 1720, 1820, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2895, 970, 1030, 1540, 1820, 1920, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3045, 1030, 1090, 1600, 1920, 2020, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3195, 1090, 1150, 1660, 2020, 2120, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 3345, 1270, 1280, 2220, 2235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3366, 1300, 1330, 2220, 2250, 2295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3429, 1330, 1360, 2235, 2295, 2340, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3492, 1420, 1480, 2250, 2385, 2475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3618, 1480, 1540, 2295, 2475, 2565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3744, 1540, 1600, 2340, 2565, 2655, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3870, 1720, 1820, 2475, 2745, 2895, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4080, 1820, 1920, 2565, 2895, 3045, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4290, 1920, 2020, 2655, 3045, 3195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 4500, 2250, 2295, 3345, 3366, 3429, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 4584, 2385, 2475, 3366, 3492, 3618, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 4752, 2475, 2565, 3429, 3618, 3744, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 4920, 2745, 2895, 3618, 3870, 4080, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 5200, 2895, 3045, 3744, 4080, 4290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 5480, 3492, 3618, 4500, 4584, 4752, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 5696, 3870, 4080, 4752, 4920, 5200, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 244, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 370, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 48, pbuffer, 694, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 84, pbuffer, 910, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 144, pbuffer, 1420, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 204, pbuffer, 1720, 100, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {1420, 1480});

                pbuffer.scale(2.0 * b_exp, {1720, 1820});

                pbuffer.scale(2.0 * b_exp, {2385, 2475});

                pbuffer.scale(2.0 * b_exp, {2745, 2895});

                pbuffer.scale(2.0 * b_exp, {3492, 3618});

                pbuffer.scale(2.0 * b_exp, {3870, 4080});

                t2cfunc::reduce(cbuffer, 304, pbuffer, 1420, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 364, pbuffer, 1720, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 464, pbuffer, 2385, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 554, pbuffer, 2745, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 704, pbuffer, 3492, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 830, pbuffer, 3870, 210, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {244, 262});

                pbuffer.scale(2.0 * a_exp, {370, 400});

                pbuffer.scale(2.0 * a_exp, {694, 730});

                pbuffer.scale(2.0 * a_exp, {910, 970});

                pbuffer.scale(a_exp / b_exp, {1420, 1480});

                pbuffer.scale(a_exp / b_exp, {1720, 1820});

                pbuffer.scale(a_exp / b_exp, {2385, 2475});

                pbuffer.scale(a_exp / b_exp, {2745, 2895});

                pbuffer.scale(a_exp / b_exp, {3492, 3618});

                pbuffer.scale(a_exp / b_exp, {3870, 4080});

                t2cfunc::reduce(cbuffer, 1040, pbuffer, 244, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1058, pbuffer, 370, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1088, pbuffer, 694, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1124, pbuffer, 910, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1184, pbuffer, 1420, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1244, pbuffer, 1720, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1344, pbuffer, 2385, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1434, pbuffer, 2745, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1584, pbuffer, 3492, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1710, pbuffer, 3870, 210, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {1420, 1480});

                pbuffer.scale(2.0 * b_exp, {1720, 1820});

                pbuffer.scale(2.0 * b_exp, {2385, 2475});

                pbuffer.scale(2.0 * b_exp, {2745, 2895});

                pbuffer.scale(2.0 * b_exp, {3492, 3618});

                pbuffer.scale(2.0 * b_exp, {3870, 4080});

                pbuffer.scale(4.0 * a_exp * b_exp, {4584, 4752});

                pbuffer.scale(4.0 * a_exp * b_exp, {4920, 5200});

                pbuffer.scale(4.0 * a_exp * b_exp, {5480, 5696});

                pbuffer.scale(4.0 * a_exp * b_exp, {5696, 6056});

                t2cfunc::reduce(cbuffer, 1920, pbuffer, 1420, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1980, pbuffer, 1720, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2080, pbuffer, 2385, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2170, pbuffer, 2745, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2320, pbuffer, 3492, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2446, pbuffer, 3870, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2656, pbuffer, 4584, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2824, pbuffer, 4920, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3104, pbuffer, 5480, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3320, pbuffer, 5696, 360, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 18, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 54, cbuffer, 48, 84, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 486, cbuffer, 144, 204, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2016, cbuffer, 304, 364, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2196, cbuffer, 464, 554, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2466, cbuffer, 704, 830, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2844, cbuffer, 1040, 1058, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2898, cbuffer, 1088, 1124, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3330, cbuffer, 1184, 1244, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4050, cbuffer, 1344, 1434, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5130, cbuffer, 1584, 1710, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8154, cbuffer, 1920, 1980, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8334, cbuffer, 2080, 2170, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8604, cbuffer, 2320, 2446, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8982, cbuffer, 2656, 2824, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 9486, cbuffer, 3104, 3320, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<1, 2>(skbuffer, 0, ckbuffer, 0, 0, 1);

            t4cfunc::ket_transform<1, 2>(skbuffer, 45, ckbuffer, 54, 0, 2);

            t4cfunc::ket_transform<1, 2>(skbuffer, 1215, ckbuffer, 486, 0, 3);

            t4cfunc::ket_transform<1, 2>(skbuffer, 50145, ckbuffer, 2016, 0, 3);

            t4cfunc::ket_transform<1, 2>(skbuffer, 50295, ckbuffer, 2196, 0, 4);

            t4cfunc::ket_transform<1, 2>(skbuffer, 50520, ckbuffer, 2466, 0, 5);

            t4cfunc::ket_transform<1, 2>(skbuffer, 50835, ckbuffer, 2844, 0, 1);

            t4cfunc::ket_transform<1, 2>(skbuffer, 50880, ckbuffer, 2898, 0, 2);

            t4cfunc::ket_transform<1, 2>(skbuffer, 51240, ckbuffer, 3330, 0, 3);

            t4cfunc::ket_transform<1, 2>(skbuffer, 51840, ckbuffer, 4050, 0, 4);

            t4cfunc::ket_transform<1, 2>(skbuffer, 52740, ckbuffer, 5130, 0, 5);

            t4cfunc::ket_transform<1, 2>(skbuffer, 56655, ckbuffer, 8154, 0, 3);

            t4cfunc::ket_transform<1, 2>(skbuffer, 56805, ckbuffer, 8334, 0, 4);

            t4cfunc::ket_transform<1, 2>(skbuffer, 57030, ckbuffer, 8604, 0, 5);

            t4cfunc::ket_transform<1, 2>(skbuffer, 57345, ckbuffer, 8982, 0, 6);

            t4cfunc::ket_transform<1, 2>(skbuffer, 57765, ckbuffer, 9486, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 8700, 45, 1215, r_ab, 1, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 55260, 50880, 51240, r_ab, 1, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 55530, 51240, 51840, r_ab, 1, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 55980, 51840, 52740, r_ab, 1, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 9780, 45, 55260, 55530, r_ab, 1, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 14370, 1215, 55530, 55980, r_ab, 1, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ddxx(skbuffer, 27465, 8700, 9780, 14370, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 135, 0, 50145, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 1365, 45, 50295, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 3165, 1215, 50520, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 8970, 45, 135, 1365, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 13020, 1215, 1365, 3165, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ddxx(skbuffer, 25845, 8700, 8970, 13020, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 50970, 50835, 56655, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 51390, 50880, 56805, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 52065, 51240, 57030, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 53055, 51840, 57345, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 54000, 52740, 57765, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 405, 50880, 50970, 51390, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 1815, 51240, 51390, 52065, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 3840, 51840, 52065, 53055, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 5865, 52740, 53055, 54000, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 10590, 135, 55260, 405, 1815, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 15720, 1365, 55530, 1815, 3840, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 19770, 3165, 55980, 3840, 5865, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 29085, 8970, 9780, 10590, 15720, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 33945, 13020, 14370, 15720, 19770, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_fdxx(skbuffer, 42045, 25845, 27465, 29085, 33945, r_ab, 1, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 42045, 1, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 525, skbuffer, 42945, 1, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1050, skbuffer, 43845, 1, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1575, skbuffer, 44745, 1, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2100, skbuffer, 45645, 1, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2625, skbuffer, 46545, 1, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 3150, skbuffer, 47445, 1, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 3675, skbuffer, 48345, 1, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 4200, skbuffer, 49245, 1, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 1, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFDPD_hpp */
