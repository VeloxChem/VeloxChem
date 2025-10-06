#ifndef ElectronRepulsionGeom1100RecFPPP_hpp
#define ElectronRepulsionGeom1100RecFPPP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0100ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecFPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FP|1/|r-r'||PP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_fppp(T& distributor,
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

    CSimdArray<double> pbuffer(2022, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1404, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(3402, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(20601, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1701, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 75, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 78, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 81, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 84, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 87, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 90, 1, 9, 12, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 99, 2, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 108, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 117, 4, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 126, 5, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 135, 6, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 144, 12, 33, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 162, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 180, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 198, 21, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 216, 24, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 234, 27, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 252, 1, 2, 75, 78, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 258, 2, 3, 78, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 264, 3, 4, 81, 84, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 270, 4, 5, 84, 87, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 276, 9, 12, 75, 90, 99, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 294, 12, 15, 78, 99, 108, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 312, 15, 18, 81, 108, 117, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 330, 18, 21, 84, 117, 126, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 348, 21, 24, 87, 126, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 366, 33, 39, 99, 144, 162, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 402, 39, 45, 108, 162, 180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 438, 45, 51, 117, 180, 198, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 474, 51, 57, 126, 198, 216, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 510, 57, 63, 135, 216, 234, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 546, 75, 78, 252, 258, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 556, 78, 81, 258, 264, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 566, 81, 84, 264, 270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 576, 90, 99, 252, 276, 294, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 606, 99, 108, 258, 294, 312, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 636, 108, 117, 264, 312, 330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 666, 117, 126, 270, 330, 348, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 696, 144, 162, 294, 366, 402, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 756, 162, 180, 312, 402, 438, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 816, 180, 198, 330, 438, 474, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 876, 198, 216, 348, 474, 510, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 936, 252, 258, 546, 556, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 951, 258, 264, 556, 566, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 966, 276, 294, 546, 576, 606, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1011, 294, 312, 556, 606, 636, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1056, 312, 330, 566, 636, 666, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1101, 366, 402, 606, 696, 756, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1191, 402, 438, 636, 756, 816, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1281, 438, 474, 666, 816, 876, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 1371, 546, 556, 936, 951, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1392, 576, 606, 936, 966, 1011, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1455, 606, 636, 951, 1011, 1056, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 1518, 696, 756, 1011, 1101, 1191, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 1644, 756, 816, 1056, 1191, 1281, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 1770, 966, 1011, 1371, 1392, 1455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 1854, 1101, 1191, 1455, 1518, 1644, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 9, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 90, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 144, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 276, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 54, pbuffer, 366, 36, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {276, 294});

                pbuffer.scale(2.0 * b_exp, {366, 402});

                pbuffer.scale(2.0 * b_exp, {576, 606});

                pbuffer.scale(2.0 * b_exp, {696, 756});

                pbuffer.scale(2.0 * b_exp, {966, 1011});

                pbuffer.scale(2.0 * b_exp, {1101, 1191});

                t2cfunc::reduce(cbuffer, 90, pbuffer, 276, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 108, pbuffer, 366, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 144, pbuffer, 576, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 174, pbuffer, 696, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 234, pbuffer, 966, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 279, pbuffer, 1101, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {9, 12});

                pbuffer.scale(2.0 * a_exp, {33, 39});

                pbuffer.scale(2.0 * a_exp, {90, 99});

                pbuffer.scale(2.0 * a_exp, {144, 162});

                pbuffer.scale(a_exp / b_exp, {276, 294});

                pbuffer.scale(a_exp / b_exp, {366, 402});

                pbuffer.scale(a_exp / b_exp, {576, 606});

                pbuffer.scale(a_exp / b_exp, {696, 756});

                pbuffer.scale(a_exp / b_exp, {966, 1011});

                pbuffer.scale(a_exp / b_exp, {1101, 1191});

                t2cfunc::reduce(cbuffer, 369, pbuffer, 9, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 372, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 378, pbuffer, 90, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 387, pbuffer, 144, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 405, pbuffer, 276, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 423, pbuffer, 366, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 459, pbuffer, 576, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 489, pbuffer, 696, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 549, pbuffer, 966, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 594, pbuffer, 1101, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {276, 294});

                pbuffer.scale(2.0 * b_exp, {366, 402});

                pbuffer.scale(2.0 * b_exp, {576, 606});

                pbuffer.scale(2.0 * b_exp, {696, 756});

                pbuffer.scale(2.0 * b_exp, {966, 1011});

                pbuffer.scale(2.0 * b_exp, {1101, 1191});

                pbuffer.scale(4.0 * a_exp * b_exp, {1392, 1455});

                pbuffer.scale(4.0 * a_exp * b_exp, {1518, 1644});

                pbuffer.scale(4.0 * a_exp * b_exp, {1770, 1854});

                pbuffer.scale(4.0 * a_exp * b_exp, {1854, 2022});

                t2cfunc::reduce(cbuffer, 684, pbuffer, 276, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 702, pbuffer, 366, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 738, pbuffer, 576, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 768, pbuffer, 696, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 828, pbuffer, 966, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 873, pbuffer, 1101, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 963, pbuffer, 1392, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1026, pbuffer, 1518, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1152, pbuffer, 1770, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1236, pbuffer, 1854, 168, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 0, cbuffer, 0, 3, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 9, cbuffer, 9, 18, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 117, cbuffer, 36, 54, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 603, cbuffer, 90, 108, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 657, cbuffer, 144, 174, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 747, cbuffer, 234, 279, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 882, cbuffer, 369, 372, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 891, cbuffer, 378, 387, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 999, cbuffer, 405, 423, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1215, cbuffer, 459, 489, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1575, cbuffer, 549, 594, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2682, cbuffer, 684, 702, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2736, cbuffer, 738, 768, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2826, cbuffer, 828, 873, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2961, cbuffer, 963, 1026, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3150, cbuffer, 1152, 1236, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 0, ckbuffer, 0, 0, 0);

            t4cfunc::ket_transform<1, 1>(skbuffer, 9, ckbuffer, 9, 0, 1);

            t4cfunc::ket_transform<1, 1>(skbuffer, 360, ckbuffer, 117, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 17289, ckbuffer, 603, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 17343, ckbuffer, 657, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 17433, ckbuffer, 747, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 17568, ckbuffer, 882, 0, 0);

            t4cfunc::ket_transform<1, 1>(skbuffer, 17577, ckbuffer, 891, 0, 1);

            t4cfunc::ket_transform<1, 1>(skbuffer, 17685, ckbuffer, 999, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 17901, ckbuffer, 1215, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 18261, ckbuffer, 1575, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 19881, ckbuffer, 2682, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 19935, ckbuffer, 2736, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 20025, ckbuffer, 2826, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 20160, ckbuffer, 2961, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 20349, ckbuffer, 3150, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 3357, 9, 360, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 19368, 17577, 17685, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 19449, 17685, 17901, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 19611, 17901, 18261, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 3681, 9, 19368, 19449, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 5139, 360, 19449, 19611, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dpxx(skbuffer, 9999, 3357, 3681, 5139, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 36, 0, 17289, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 414, 9, 17343, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 1062, 360, 17433, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 3438, 9, 36, 414, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 4653, 360, 414, 1062, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dpxx(skbuffer, 9513, 3357, 3438, 4653, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 17604, 17568, 19881, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 17739, 17577, 19935, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 17991, 17685, 20025, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 18396, 17901, 20160, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 18801, 18261, 20349, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 117, 17577, 17604, 17739, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 576, 17685, 17739, 17991, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 1332, 17901, 17991, 18396, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 2142, 18261, 18396, 18801, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 3924, 36, 19368, 117, 576, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 5625, 414, 19449, 576, 1332, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 7083, 1062, 19611, 1332, 2142, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 10485, 3438, 3681, 3924, 5625, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 11943, 4653, 5139, 5625, 7083, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_fpxx(skbuffer, 14859, 9513, 9999, 10485, 11943, r_ab, 1, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 14859, 1, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 189, skbuffer, 15129, 1, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 378, skbuffer, 15399, 1, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 567, skbuffer, 15669, 1, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 756, skbuffer, 15939, 1, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 945, skbuffer, 16209, 1, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1134, skbuffer, 16479, 1, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1323, skbuffer, 16749, 1, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1512, skbuffer, 17019, 1, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 1, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFPPP_hpp */
