#ifndef ElectronRepulsionGeom1100RecFFSP_hpp
#define ElectronRepulsionGeom1100RecFFSP_hpp

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
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSKSP.hpp"
#include "ElectronRepulsionPrimRecSKSS.hpp"
#include "ElectronRepulsionPrimRecSLSP.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FF|1/|r-r'||SP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_ffsp(T& distributor,
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

    CSimdArray<double> pbuffer(1817, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(960, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(17772, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1323, 1);

    // setup Boys fuction data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_data(11, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 10, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 10);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 10, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 25, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 28, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 31, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 34, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 37, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 40, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 43, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 46, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 49, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 52, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 55, 7, 8, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 58, 1, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 67, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 76, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 85, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 94, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 103, 6, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 112, 7, 28, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 121, 8, 31, 34, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 130, 1, 2, 37, 40, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 136, 2, 3, 40, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 142, 3, 4, 43, 46, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 148, 4, 5, 46, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 154, 5, 6, 49, 52, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 160, 6, 7, 52, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 166, 10, 13, 37, 58, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 184, 13, 16, 40, 67, 76, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 202, 16, 19, 43, 76, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 220, 19, 22, 46, 85, 94, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 238, 22, 25, 49, 94, 103, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 256, 25, 28, 52, 103, 112, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 274, 28, 31, 55, 112, 121, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 292, 37, 40, 130, 136, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 302, 40, 43, 136, 142, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 312, 43, 46, 142, 148, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 322, 46, 49, 148, 154, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 332, 49, 52, 154, 160, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 342, 58, 67, 130, 166, 184, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 372, 67, 76, 136, 184, 202, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 402, 76, 85, 142, 202, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 432, 85, 94, 148, 220, 238, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 462, 94, 103, 154, 238, 256, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 492, 103, 112, 160, 256, 274, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 522, 130, 136, 292, 302, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 537, 136, 142, 302, 312, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 552, 142, 148, 312, 322, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 567, 148, 154, 322, 332, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 582, 166, 184, 292, 342, 372, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 627, 184, 202, 302, 372, 402, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 672, 202, 220, 312, 402, 432, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 717, 220, 238, 322, 432, 462, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 762, 238, 256, 332, 462, 492, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 807, 292, 302, 522, 537, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 828, 302, 312, 537, 552, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 849, 312, 322, 552, 567, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 870, 342, 372, 522, 582, 627, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 933, 372, 402, 537, 627, 672, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 996, 402, 432, 552, 672, 717, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1059, 432, 462, 567, 717, 762, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 1122, 522, 537, 807, 828, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 1150, 537, 552, 828, 849, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 1178, 582, 627, 807, 870, 933, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 1262, 627, 672, 828, 933, 996, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 1346, 672, 717, 849, 996, 1059, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_skss(pbuffer, 1430, 807, 828, 1122, 1150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 1466, 870, 933, 1122, 1178, 1262, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 1574, 933, 996, 1150, 1262, 1346, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsp(pbuffer, 1682, 1178, 1262, 1430, 1466, 1574, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 166, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 342, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 48, pbuffer, 582, 45, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {582, 627});

                pbuffer.scale(2.0 * b_exp, {870, 933});

                pbuffer.scale(2.0 * b_exp, {1178, 1262});

                t2cfunc::reduce(cbuffer, 93, pbuffer, 582, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 138, pbuffer, 870, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 201, pbuffer, 1178, 84, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {166, 184});

                pbuffer.scale(2.0 * a_exp, {342, 372});

                pbuffer.scale(a_exp / b_exp, {582, 627});

                pbuffer.scale(a_exp / b_exp, {870, 933});

                pbuffer.scale(a_exp / b_exp, {1178, 1262});

                t2cfunc::reduce(cbuffer, 285, pbuffer, 166, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 303, pbuffer, 342, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 333, pbuffer, 582, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 378, pbuffer, 870, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 441, pbuffer, 1178, 84, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {582, 627});

                pbuffer.scale(2.0 * b_exp, {870, 933});

                pbuffer.scale(2.0 * b_exp, {1178, 1262});

                pbuffer.scale(4.0 * a_exp * b_exp, {1466, 1574});

                pbuffer.scale(4.0 * a_exp * b_exp, {1682, 1817});

                t2cfunc::reduce(cbuffer, 525, pbuffer, 582, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 570, pbuffer, 870, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 633, pbuffer, 1178, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 717, pbuffer, 1466, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 825, pbuffer, 1682, 135, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 1>(skbuffer, 0, cbuffer, 0, 0, 2);

            t4cfunc::ket_transform<0, 1>(skbuffer, 18, cbuffer, 18, 0, 3);

            t4cfunc::ket_transform<0, 1>(skbuffer, 408, cbuffer, 48, 0, 4);

            t4cfunc::ket_transform<0, 1>(skbuffer, 15501, cbuffer, 93, 0, 4);

            t4cfunc::ket_transform<0, 1>(skbuffer, 15546, cbuffer, 138, 0, 5);

            t4cfunc::ket_transform<0, 1>(skbuffer, 15609, cbuffer, 201, 0, 6);

            t4cfunc::ket_transform<0, 1>(skbuffer, 15693, cbuffer, 285, 0, 2);

            t4cfunc::ket_transform<0, 1>(skbuffer, 15711, cbuffer, 303, 0, 3);

            t4cfunc::ket_transform<0, 1>(skbuffer, 15831, cbuffer, 333, 0, 4);

            t4cfunc::ket_transform<0, 1>(skbuffer, 16011, cbuffer, 378, 0, 5);

            t4cfunc::ket_transform<0, 1>(skbuffer, 16263, cbuffer, 441, 0, 6);

            t4cfunc::ket_transform<0, 1>(skbuffer, 17337, cbuffer, 525, 0, 4);

            t4cfunc::ket_transform<0, 1>(skbuffer, 17382, cbuffer, 570, 0, 5);

            t4cfunc::ket_transform<0, 1>(skbuffer, 17445, cbuffer, 633, 0, 6);

            t4cfunc::ket_transform<0, 1>(skbuffer, 17529, cbuffer, 717, 0, 7);

            t4cfunc::ket_transform<0, 1>(skbuffer, 17637, cbuffer, 825, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 2505, 18, 408, r_ab, 0, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 16923, 15711, 15831, r_ab, 0, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 17013, 15831, 16011, r_ab, 0, 1);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 17148, 16011, 16263, r_ab, 0, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 2865, 18, 16923, 17013, r_ab, 0, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 4350, 408, 17013, 17148, r_ab, 0, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 8211, 2505, 2865, 4350, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 48, 0, 15501, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 453, 18, 15546, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 993, 408, 15609, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 2595, 18, 48, 453, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pgxx(skbuffer, 3945, 408, 453, 993, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dfxx(skbuffer, 7671, 2505, 2595, 3945, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 15741, 15693, 17337, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 15876, 15711, 17382, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 16074, 15831, 17445, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 16347, 16011, 17529, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_skxx(skbuffer, 16599, 16263, 17637, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 138, 15711, 15741, 15876, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 588, 15831, 15876, 16074, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 1182, 16011, 16074, 16347, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sixx(skbuffer, 1749, 16263, 16347, 16599, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 3135, 48, 16923, 138, 588, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 4755, 453, 17013, 588, 1182, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_phxx(skbuffer, 5970, 993, 17148, 1182, 1749, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 8751, 2595, 2865, 3135, 4755, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dgxx(skbuffer, 10371, 3945, 4350, 4755, 5970, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ffxx(skbuffer, 12801, 7671, 8211, 8751, 10371, r_ab, 0, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 12801, 0, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 147, skbuffer, 13101, 0, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 294, skbuffer, 13401, 0, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 441, skbuffer, 13701, 0, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 588, skbuffer, 14001, 0, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 735, skbuffer, 14301, 0, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 882, skbuffer, 14601, 0, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1029, skbuffer, 14901, 0, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1176, skbuffer, 15201, 0, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 0, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFFSP_hpp */
