#ifndef ElectronRepulsionGeom1010RecFPDS_hpp
#define ElectronRepulsionGeom1010RecFPDS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
#include "ElectronRepulsionGeom1010ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FP|1/|r-r'||DS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fpds(T& distributor,
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

    CSimdArray<double> pbuffer(2535, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1776, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(6438, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(9660, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(945, 1);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 9, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 9);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 9, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 12, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 21, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 24, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 27, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 30, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 33, 0, 1, 9, 12, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 39, 1, 2, 12, 15, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 45, 2, 3, 15, 18, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 51, 3, 4, 18, 21, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 57, 4, 5, 21, 24, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 63, 5, 6, 24, 27, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 69, 6, 7, 27, 30, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 75, 9, 12, 33, 39, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 12, 15, 39, 45, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 15, 18, 45, 51, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 18, 21, 51, 57, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 21, 24, 57, 63, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 24, 27, 63, 69, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 135, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 138, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 141, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 144, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 147, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 150, 1, 9, 12, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 159, 2, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 168, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 177, 4, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 186, 5, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 195, 12, 33, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 213, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 231, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 249, 21, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 267, 24, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 285, 39, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 315, 45, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 345, 51, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 375, 57, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 405, 63, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 435, 0, 1, 135, 138, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 441, 1, 2, 138, 141, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 447, 2, 3, 141, 144, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 453, 3, 4, 144, 147, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 459, 9, 12, 138, 150, 159, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 477, 12, 15, 141, 159, 168, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 495, 15, 18, 144, 168, 177, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 513, 18, 21, 147, 177, 186, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 531, 33, 39, 159, 195, 213, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 567, 39, 45, 168, 213, 231, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 603, 45, 51, 177, 231, 249, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 639, 51, 57, 186, 249, 267, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 675, 75, 85, 213, 285, 315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 735, 85, 95, 231, 315, 345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 795, 95, 105, 249, 345, 375, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 855, 105, 115, 267, 375, 405, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 915, 135, 138, 435, 441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 925, 138, 141, 441, 447, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 935, 141, 144, 447, 453, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 945, 150, 159, 441, 459, 477, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 975, 159, 168, 447, 477, 495, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1005, 168, 177, 453, 495, 513, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1035, 195, 213, 477, 531, 567, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1095, 213, 231, 495, 567, 603, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1155, 231, 249, 513, 603, 639, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1215, 285, 315, 567, 675, 735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1315, 315, 345, 603, 735, 795, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1415, 345, 375, 639, 795, 855, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1515, 435, 441, 915, 925, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1530, 441, 447, 925, 935, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1545, 459, 477, 925, 945, 975, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1590, 477, 495, 935, 975, 1005, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1635, 531, 567, 975, 1035, 1095, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1725, 567, 603, 1005, 1095, 1155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 1815, 675, 735, 1095, 1215, 1315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 1965, 735, 795, 1155, 1315, 1415, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 2115, 915, 925, 1515, 1530, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2136, 945, 975, 1530, 1545, 1590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2199, 1035, 1095, 1590, 1635, 1725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 2325, 1215, 1315, 1725, 1815, 1965, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 135, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 150, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12, pbuffer, 435, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 459, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 915, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 46, pbuffer, 945, 30, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {135, 138});

                pbuffer.scale(2.0 * a_exp, {150, 159});

                pbuffer.scale(2.0 * a_exp, {435, 441});

                pbuffer.scale(2.0 * a_exp, {459, 477});

                pbuffer.scale(2.0 * a_exp, {915, 925});

                pbuffer.scale(2.0 * a_exp, {945, 975});

                pbuffer.scale(2.0 * a_exp, {1515, 1530});

                pbuffer.scale(2.0 * a_exp, {1545, 1590});

                pbuffer.scale(2.0 * a_exp, {2115, 2136});

                pbuffer.scale(2.0 * a_exp, {2136, 2199});

                t2cfunc::reduce(cbuffer, 456, pbuffer, 135, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 459, pbuffer, 150, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 468, pbuffer, 435, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 474, pbuffer, 459, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 492, pbuffer, 915, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 502, pbuffer, 945, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 532, pbuffer, 1515, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 547, pbuffer, 1545, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 592, pbuffer, 2115, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 613, pbuffer, 2136, 63, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {135, 138});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {150, 159});

                pbuffer.scale(pfactors, 0, 2.0, {195, 213});

                pbuffer.scale(pfactors, 0, 2.0, {285, 315});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {435, 441});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {459, 477});

                pbuffer.scale(pfactors, 0, 2.0, {531, 567});

                pbuffer.scale(pfactors, 0, 2.0, {675, 735});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {915, 925});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {945, 975});

                pbuffer.scale(pfactors, 0, 2.0, {1035, 1095});

                pbuffer.scale(pfactors, 0, 2.0, {1215, 1315});

                t2cfunc::reduce(cbuffer, 76, pbuffer, 135, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 79, pbuffer, 150, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 88, pbuffer, 195, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 106, pbuffer, 285, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 136, pbuffer, 435, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 142, pbuffer, 459, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 531, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 196, pbuffer, 675, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 256, pbuffer, 915, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 266, pbuffer, 945, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 296, pbuffer, 1035, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 356, pbuffer, 1215, 100, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {135, 138});

                pbuffer.scale(2.0 * a_exp, {150, 159});

                pbuffer.scale(2.0 * a_exp, {195, 213});

                pbuffer.scale(2.0 * a_exp, {285, 315});

                pbuffer.scale(2.0 * a_exp, {435, 441});

                pbuffer.scale(2.0 * a_exp, {459, 477});

                pbuffer.scale(2.0 * a_exp, {531, 567});

                pbuffer.scale(2.0 * a_exp, {675, 735});

                pbuffer.scale(2.0 * a_exp, {915, 925});

                pbuffer.scale(2.0 * a_exp, {945, 975});

                pbuffer.scale(2.0 * a_exp, {1035, 1095});

                pbuffer.scale(2.0 * a_exp, {1215, 1315});

                pbuffer.scale(pfactors, 0, 2.0, {1515, 1530});

                pbuffer.scale(pfactors, 0, 2.0, {1545, 1590});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1635, 1725});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1815, 1965});

                pbuffer.scale(pfactors, 0, 2.0, {2115, 2136});

                pbuffer.scale(pfactors, 0, 2.0, {2136, 2199});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2199, 2325});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2325, 2535});

                t2cfunc::reduce(cbuffer, 676, pbuffer, 135, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 679, pbuffer, 150, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 688, pbuffer, 195, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 706, pbuffer, 285, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 736, pbuffer, 435, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 742, pbuffer, 459, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 760, pbuffer, 531, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 796, pbuffer, 675, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 856, pbuffer, 915, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 866, pbuffer, 945, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 896, pbuffer, 1035, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 956, pbuffer, 1215, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1056, pbuffer, 1515, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1071, pbuffer, 1545, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1116, pbuffer, 1635, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1206, pbuffer, 1815, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1356, pbuffer, 2115, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1377, pbuffer, 2136, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1440, pbuffer, 2199, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1566, pbuffer, 2325, 210, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 90, cbuffer, 0, 3, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 441, cbuffer, 12, 18, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1083, cbuffer, 36, 46, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1743, cbuffer, 456, 459, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2094, cbuffer, 468, 474, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2736, cbuffer, 492, 502, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 3756, cbuffer, 532, 547, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 5241, cbuffer, 592, 613, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 76, 79, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9, cbuffer, 79, 88, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 36, cbuffer, 88, 106, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 99, cbuffer, 0, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 126, cbuffer, 3, 9, 36, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 207, 90, 99, 126, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 261, cbuffer, 136, 142, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 279, cbuffer, 142, 160, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 333, cbuffer, 160, 196, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 459, cbuffer, 12, 261, 279, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 513, cbuffer, 18, 279, 333, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 675, 441, 459, 513, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 783, cbuffer, 256, 266, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 813, cbuffer, 266, 296, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 903, cbuffer, 296, 356, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1113, cbuffer, 36, 783, 813, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1203, cbuffer, 46, 813, 903, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1473, 1083, 1113, 1203, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1653, cbuffer, 676, 679, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1662, cbuffer, 679, 688, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1689, cbuffer, 688, 706, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1752, cbuffer, 456, 1653, 1662, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1779, cbuffer, 459, 1662, 1689, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1860, 1743, 1752, 1779, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1914, cbuffer, 736, 742, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1932, cbuffer, 742, 760, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1986, cbuffer, 760, 796, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2112, cbuffer, 468, 1914, 1932, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2166, cbuffer, 474, 1932, 1986, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 2328, 2094, 2112, 2166, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2436, cbuffer, 856, 866, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2466, cbuffer, 866, 896, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2556, cbuffer, 896, 956, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2766, cbuffer, 492, 2436, 2466, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2856, cbuffer, 502, 2466, 2556, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 3126, 2736, 2766, 2856, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 3306, cbuffer, 1056, 1071, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3351, cbuffer, 1071, 1116, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3486, cbuffer, 1116, 1206, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 3801, cbuffer, 532, 3306, 3351, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3936, cbuffer, 547, 3351, 3486, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 4341, 3756, 3801, 3936, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 4611, cbuffer, 1356, 1377, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4674, cbuffer, 1377, 1440, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4863, cbuffer, 1440, 1566, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 5304, cbuffer, 592, 4611, 4674, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 5493, cbuffer, 613, 4674, 4863, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 6060, 5241, 5304, 5493, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 0, ckbuffer, 207, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15, ckbuffer, 225, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 30, ckbuffer, 243, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 180, ckbuffer, 675, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 210, ckbuffer, 711, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 240, ckbuffer, 747, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 540, ckbuffer, 1473, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 590, ckbuffer, 1533, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 640, ckbuffer, 1593, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 8835, ckbuffer, 1860, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 8850, ckbuffer, 1878, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 8865, ckbuffer, 1896, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 8880, ckbuffer, 2328, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 8910, ckbuffer, 2364, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 8940, ckbuffer, 2400, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 8970, ckbuffer, 3126, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9020, ckbuffer, 3186, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9070, ckbuffer, 3246, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9120, ckbuffer, 4341, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9195, ckbuffer, 4431, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9270, ckbuffer, 4521, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9345, ckbuffer, 6060, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9450, ckbuffer, 6186, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9555, ckbuffer, 6312, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1815, 0, 180, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1860, 15, 210, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1905, 30, 240, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2355, 180, 540, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2445, 210, 590, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2535, 240, 640, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 4785, 1815, 2355, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 4875, 1860, 2445, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 4965, 1905, 2535, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 45, 8835, 8880, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 270, 8880, 8970, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 690, 8970, 9120, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 1140, 9120, 9345, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 1950, 0, 45, 270, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 2625, 180, 270, 690, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 3435, 540, 690, 1140, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 5055, 1815, 1950, 2625, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 5865, 2355, 2625, 3435, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fpxx(skbuffer, 7485, 4785, 5055, 5865, r_ab, 2, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 7485, 2, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 105, skbuffer, 7635, 2, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 210, skbuffer, 7785, 2, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 315, skbuffer, 7935, 2, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 420, skbuffer, 8085, 2, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 525, skbuffer, 8235, 2, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 630, skbuffer, 8385, 2, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 735, skbuffer, 8535, 2, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 840, skbuffer, 8685, 2, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 2, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFPDS_hpp */
