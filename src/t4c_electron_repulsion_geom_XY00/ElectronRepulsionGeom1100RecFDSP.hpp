#ifndef ElectronRepulsionGeom1100RecFDSP_hpp
#define ElectronRepulsionGeom1100RecFDSP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FD|1/|r-r'||SP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_fdsp(T& distributor,
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

    CSimdArray<double> pbuffer(1202, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(690, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(11661, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 33, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 36, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 39, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 42, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 45, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 48, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 51, 1, 9, 12, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 60, 2, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 69, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 78, 4, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 87, 5, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 96, 6, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 105, 7, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 114, 1, 2, 33, 36, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 120, 2, 3, 36, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 126, 3, 4, 39, 42, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 132, 4, 5, 42, 45, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 138, 5, 6, 45, 48, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 144, 9, 12, 33, 51, 60, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 162, 12, 15, 36, 60, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 180, 15, 18, 39, 69, 78, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 198, 18, 21, 42, 78, 87, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 216, 21, 24, 45, 87, 96, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 234, 24, 27, 48, 96, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 252, 33, 36, 114, 120, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 262, 36, 39, 120, 126, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 272, 39, 42, 126, 132, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 282, 42, 45, 132, 138, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 292, 51, 60, 114, 144, 162, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 322, 60, 69, 120, 162, 180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 352, 69, 78, 126, 180, 198, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 382, 78, 87, 132, 198, 216, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 412, 87, 96, 138, 216, 234, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 442, 114, 120, 252, 262, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 457, 120, 126, 262, 272, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 472, 126, 132, 272, 282, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 487, 144, 162, 252, 292, 322, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 532, 162, 180, 262, 322, 352, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 577, 180, 198, 272, 352, 382, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 622, 198, 216, 282, 382, 412, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 667, 252, 262, 442, 457, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 688, 262, 272, 457, 472, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 709, 292, 322, 442, 487, 532, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 772, 322, 352, 457, 532, 577, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 835, 352, 382, 472, 577, 622, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 898, 442, 457, 667, 688, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 926, 487, 532, 667, 709, 772, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 1010, 532, 577, 688, 772, 835, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 1094, 709, 772, 898, 926, 1010, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 51, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 144, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 27, pbuffer, 292, 30, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {292, 322});

                pbuffer.scale(2.0 * b_exp, {487, 532});

                pbuffer.scale(2.0 * b_exp, {709, 772});

                t2cfunc::reduce(cbuffer, 57, pbuffer, 292, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 87, pbuffer, 487, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 132, pbuffer, 709, 63, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {51, 60});

                pbuffer.scale(2.0 * a_exp, {144, 162});

                pbuffer.scale(a_exp / b_exp, {292, 322});

                pbuffer.scale(a_exp / b_exp, {487, 532});

                pbuffer.scale(a_exp / b_exp, {709, 772});

                t2cfunc::reduce(cbuffer, 195, pbuffer, 51, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 204, pbuffer, 144, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 222, pbuffer, 292, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 252, pbuffer, 487, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 297, pbuffer, 709, 63, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {292, 322});

                pbuffer.scale(2.0 * b_exp, {487, 532});

                pbuffer.scale(2.0 * b_exp, {709, 772});

                pbuffer.scale(4.0 * a_exp * b_exp, {926, 1010});

                pbuffer.scale(4.0 * a_exp * b_exp, {1094, 1202});

                t2cfunc::reduce(cbuffer, 360, pbuffer, 292, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 390, pbuffer, 487, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 435, pbuffer, 709, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 498, pbuffer, 926, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 582, pbuffer, 1094, 108, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 1>(skbuffer, 0, cbuffer, 0, 0, 1);

            t4cfunc::ket_transform<0, 1>(skbuffer, 9, cbuffer, 9, 0, 2);

            t4cfunc::ket_transform<0, 1>(skbuffer, 243, cbuffer, 27, 0, 3);

            t4cfunc::ket_transform<0, 1>(skbuffer, 10029, cbuffer, 57, 0, 3);

            t4cfunc::ket_transform<0, 1>(skbuffer, 10059, cbuffer, 87, 0, 4);

            t4cfunc::ket_transform<0, 1>(skbuffer, 10104, cbuffer, 132, 0, 5);

            t4cfunc::ket_transform<0, 1>(skbuffer, 10167, cbuffer, 195, 0, 1);

            t4cfunc::ket_transform<0, 1>(skbuffer, 10176, cbuffer, 204, 0, 2);

            t4cfunc::ket_transform<0, 1>(skbuffer, 10248, cbuffer, 222, 0, 3);

            t4cfunc::ket_transform<0, 1>(skbuffer, 10368, cbuffer, 252, 0, 4);

            t4cfunc::ket_transform<0, 1>(skbuffer, 10548, cbuffer, 297, 0, 5);

            t4cfunc::ket_transform<0, 1>(skbuffer, 11331, cbuffer, 360, 0, 3);

            t4cfunc::ket_transform<0, 1>(skbuffer, 11361, cbuffer, 390, 0, 4);

            t4cfunc::ket_transform<0, 1>(skbuffer, 11406, cbuffer, 435, 0, 5);

            t4cfunc::ket_transform<0, 1>(skbuffer, 11469, cbuffer, 498, 0, 6);

            t4cfunc::ket_transform<0, 1>(skbuffer, 11553, cbuffer, 582, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 1740, 9, 243, r_ab, 0, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 11052, 10176, 10248, r_ab, 0, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 11106, 10248, 10368, r_ab, 0, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 11196, 10368, 10548, r_ab, 0, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 1956, 9, 11052, 11106, r_ab, 0, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 2874, 243, 11106, 11196, r_ab, 0, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ddxx(skbuffer, 5493, 1740, 1956, 2874, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 27, 0, 10029, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 273, 9, 10059, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 633, 243, 10104, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 1794, 9, 27, 273, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 2604, 243, 273, 633, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ddxx(skbuffer, 5169, 1740, 1794, 2604, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 10194, 10167, 11331, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 10278, 10176, 11361, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 10413, 10248, 11406, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 10611, 10368, 11469, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 10800, 10548, 11553, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 81, 10176, 10194, 10278, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 363, 10248, 10278, 10413, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 768, 10368, 10413, 10611, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 1173, 10548, 10611, 10800, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 2118, 27, 11052, 81, 363, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 3144, 273, 11106, 363, 768, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 3954, 633, 11196, 768, 1173, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 5817, 1794, 1956, 2118, 3144, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 6789, 2604, 2874, 3144, 3954, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_fdxx(skbuffer, 8409, 5169, 5493, 5817, 6789, r_ab, 0, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 8409, 0, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 105, skbuffer, 8589, 0, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 210, skbuffer, 8769, 0, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 315, skbuffer, 8949, 0, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 420, skbuffer, 9129, 0, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 525, skbuffer, 9309, 0, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 630, skbuffer, 9489, 0, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 735, skbuffer, 9669, 0, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 840, skbuffer, 9849, 0, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 0, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFDSP_hpp */
