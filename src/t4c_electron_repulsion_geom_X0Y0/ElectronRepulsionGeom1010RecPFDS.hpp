#ifndef ElectronRepulsionGeom1010RecPFDS_hpp
#define ElectronRepulsionGeom1010RecPFDS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(PF|1/|r-r'||DS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_pfds(T& distributor,
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

    CSimdArray<double> cbuffer(1344, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(4872, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(3315, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 915, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 945, 30, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {915, 925});

                pbuffer.scale(2.0 * a_exp, {945, 975});

                pbuffer.scale(2.0 * a_exp, {1515, 1530});

                pbuffer.scale(2.0 * a_exp, {1545, 1590});

                pbuffer.scale(2.0 * a_exp, {2115, 2136});

                pbuffer.scale(2.0 * a_exp, {2136, 2199});

                t2cfunc::reduce(cbuffer, 240, pbuffer, 915, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 250, pbuffer, 945, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 280, pbuffer, 1515, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 295, pbuffer, 1545, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 340, pbuffer, 2115, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 361, pbuffer, 2136, 63, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {915, 925});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {945, 975});

                pbuffer.scale(pfactors, 0, 2.0, {1035, 1095});

                pbuffer.scale(pfactors, 0, 2.0, {1215, 1315});

                t2cfunc::reduce(cbuffer, 40, pbuffer, 915, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 50, pbuffer, 945, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 80, pbuffer, 1035, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 140, pbuffer, 1215, 100, ket_width, ket_npgtos);

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

                t2cfunc::reduce(cbuffer, 424, pbuffer, 915, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 434, pbuffer, 945, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 464, pbuffer, 1035, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 524, pbuffer, 1215, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 624, pbuffer, 1515, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 639, pbuffer, 1545, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 684, pbuffer, 1635, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 774, pbuffer, 1815, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 924, pbuffer, 2115, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 945, pbuffer, 2136, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1008, pbuffer, 2199, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1134, pbuffer, 2325, 210, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 300, cbuffer, 0, 10, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1170, cbuffer, 240, 250, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2190, cbuffer, 280, 295, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 3675, cbuffer, 340, 361, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 40, 50, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 30, cbuffer, 50, 80, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 120, cbuffer, 80, 140, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 330, cbuffer, 0, 0, 30, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 420, cbuffer, 10, 30, 120, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 690, 300, 330, 420, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 870, cbuffer, 424, 434, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 900, cbuffer, 434, 464, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 990, cbuffer, 464, 524, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1200, cbuffer, 240, 870, 900, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1290, cbuffer, 250, 900, 990, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1560, 1170, 1200, 1290, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1740, cbuffer, 624, 639, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1785, cbuffer, 639, 684, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1920, cbuffer, 684, 774, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2235, cbuffer, 280, 1740, 1785, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2370, cbuffer, 295, 1785, 1920, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 2775, 2190, 2235, 2370, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 3045, cbuffer, 924, 945, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3108, cbuffer, 945, 1008, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3297, cbuffer, 1008, 1134, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 3738, cbuffer, 340, 3045, 3108, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3927, cbuffer, 361, 3108, 3297, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 4494, 3675, 3738, 3927, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 0, ckbuffer, 690, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 50, ckbuffer, 750, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 100, ckbuffer, 810, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 2625, ckbuffer, 1560, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 2675, ckbuffer, 1620, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 2725, ckbuffer, 1680, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 2775, ckbuffer, 2775, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 2850, ckbuffer, 2865, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 2925, ckbuffer, 2955, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3000, ckbuffer, 4494, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3105, ckbuffer, 4620, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3210, ckbuffer, 4746, 0, 5);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 150, 2625, 2775, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 600, 2775, 3000, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 1275, 0, 150, 600, r_ab, 2, 0);

            t4cfunc::bra_transform<1, 3>(sbuffer, 0, skbuffer, 1275, 2, 0);

            t4cfunc::bra_transform<1, 3>(sbuffer, 105, skbuffer, 1425, 2, 0);

            t4cfunc::bra_transform<1, 3>(sbuffer, 210, skbuffer, 1575, 2, 0);

            t4cfunc::bra_transform<1, 3>(sbuffer, 315, skbuffer, 1725, 2, 0);

            t4cfunc::bra_transform<1, 3>(sbuffer, 420, skbuffer, 1875, 2, 0);

            t4cfunc::bra_transform<1, 3>(sbuffer, 525, skbuffer, 2025, 2, 0);

            t4cfunc::bra_transform<1, 3>(sbuffer, 630, skbuffer, 2175, 2, 0);

            t4cfunc::bra_transform<1, 3>(sbuffer, 735, skbuffer, 2325, 2, 0);

            t4cfunc::bra_transform<1, 3>(sbuffer, 840, skbuffer, 2475, 2, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 1, 3, 2, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecPFDS_hpp */
