#ifndef ElectronRepulsionGeom1100RecDSDD_hpp
#define ElectronRepulsionGeom1100RecDSDD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPSXX.hpp"
#include "ElectronRepulsionContrRecXXDD.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0100ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DS|1/|r-r'||DD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dsdd(T& distributor,
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

    CSimdArray<double> pbuffer(2330, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1488, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(6408, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(9900, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1125, 1);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 9);

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

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 135, 33, 39, 75, 85, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 150, 39, 45, 85, 95, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 165, 45, 51, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 180, 51, 57, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 195, 57, 63, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 210, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 213, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 216, 2, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 225, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 234, 4, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 243, 12, 33, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 261, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 279, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 297, 21, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 315, 39, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 345, 45, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 375, 51, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 405, 57, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 435, 85, 135, 150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 480, 95, 150, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 525, 105, 165, 180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 570, 115, 180, 195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 615, 2, 3, 210, 213, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 621, 12, 15, 210, 216, 225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 639, 15, 18, 213, 225, 234, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 657, 33, 39, 216, 243, 261, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 693, 39, 45, 225, 261, 279, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 729, 45, 51, 234, 279, 297, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 765, 75, 85, 261, 315, 345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 825, 85, 95, 279, 345, 375, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 885, 95, 105, 297, 375, 405, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 945, 135, 150, 345, 435, 480, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1035, 150, 165, 375, 480, 525, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1125, 165, 180, 405, 525, 570, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1215, 216, 225, 615, 621, 639, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1245, 243, 261, 621, 657, 693, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1305, 261, 279, 639, 693, 729, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1365, 315, 345, 693, 765, 825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1465, 345, 375, 729, 825, 885, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1565, 435, 480, 825, 945, 1035, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1715, 480, 525, 885, 1035, 1125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1865, 657, 693, 1215, 1245, 1305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 1955, 765, 825, 1305, 1365, 1465, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 2105, 945, 1035, 1465, 1565, 1715, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 135, 15, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {243, 261});

                pbuffer.scale(2.0 * b_exp, {315, 345});

                pbuffer.scale(2.0 * b_exp, {435, 480});

                t2cfunc::reduce(cbuffer, 31, pbuffer, 243, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 49, pbuffer, 315, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 79, pbuffer, 435, 45, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {33, 39});

                pbuffer.scale(2.0 * a_exp, {75, 85});

                pbuffer.scale(2.0 * a_exp, {135, 150});

                pbuffer.scale(a_exp / b_exp, {243, 261});

                pbuffer.scale(a_exp / b_exp, {315, 345});

                pbuffer.scale(a_exp / b_exp, {435, 480});

                pbuffer.scale(2.0 * a_exp, {657, 693});

                pbuffer.scale(2.0 * a_exp, {765, 825});

                pbuffer.scale(2.0 * a_exp, {945, 1035});

                t2cfunc::reduce(cbuffer, 124, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 130, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 140, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 155, pbuffer, 243, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 173, pbuffer, 315, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 203, pbuffer, 435, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 248, pbuffer, 657, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 284, pbuffer, 765, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 344, pbuffer, 945, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {243, 261});

                pbuffer.scale(2.0 * b_exp, {315, 345});

                pbuffer.scale(2.0 * b_exp, {435, 480});

                pbuffer.scale(2.0 * b_exp, {657, 693});

                pbuffer.scale(2.0 * b_exp, {765, 825});

                pbuffer.scale(2.0 * b_exp, {945, 1035});

                pbuffer.scale(4.0 * a_exp * b_exp, {1245, 1305});

                pbuffer.scale(4.0 * a_exp * b_exp, {1365, 1465});

                pbuffer.scale(4.0 * a_exp * b_exp, {1565, 1715});

                pbuffer.scale(4.0 * a_exp * b_exp, {1865, 1955});

                pbuffer.scale(4.0 * a_exp * b_exp, {1955, 2105});

                pbuffer.scale(4.0 * a_exp * b_exp, {2105, 2330});

                t2cfunc::reduce(cbuffer, 434, pbuffer, 243, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 452, pbuffer, 315, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 482, pbuffer, 435, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 527, pbuffer, 657, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 563, pbuffer, 765, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 623, pbuffer, 945, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 713, pbuffer, 1245, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 773, pbuffer, 1365, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 873, pbuffer, 1565, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1023, pbuffer, 1865, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1113, pbuffer, 1955, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1263, pbuffer, 2105, 225, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 6, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 18, cbuffer, 6, 16, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 48, 0, 18, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 408, cbuffer, 31, 49, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 462, cbuffer, 49, 79, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 552, 408, 462, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 660, cbuffer, 124, 130, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 678, cbuffer, 130, 140, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 708, 660, 678, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 744, cbuffer, 155, 173, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 798, cbuffer, 173, 203, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 888, 744, 798, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1320, cbuffer, 248, 284, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1428, cbuffer, 284, 344, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 1608, 1320, 1428, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3552, cbuffer, 434, 452, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 3606, cbuffer, 452, 482, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 3696, 3552, 3606, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3804, cbuffer, 527, 563, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 3912, cbuffer, 563, 623, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 4092, 3804, 3912, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4308, cbuffer, 713, 773, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4488, cbuffer, 773, 873, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 4788, 4308, 4488, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5148, cbuffer, 1023, 1113, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5418, cbuffer, 1113, 1263, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 5868, 5148, 5418, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 0, ckbuffer, 48, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 6775, ckbuffer, 552, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 6850, ckbuffer, 708, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 6875, ckbuffer, 888, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 7175, ckbuffer, 1608, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 8825, ckbuffer, 3696, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 8900, ckbuffer, 4092, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 9050, ckbuffer, 4788, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 9300, ckbuffer, 5868, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 8525, 6850, 6875, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 8600, 6875, 7175, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 9675, 8825, 8900, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_psxx(skbuffer, 2500, 0, 8525, 8600, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 25, 6850, 8900, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_psxx(skbuffer, 2275, 0, 6775, 25, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 6950, 6850, 8900, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 7325, 6875, 9050, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 7775, 7175, 9300, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 250, 6875, 6950, 7325, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 925, 7175, 7325, 7775, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_psxx(skbuffer, 2725, 6775, 8525, 9675, 250, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 3400, 25, 8600, 250, 925, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dsxx(skbuffer, 5425, 2275, 2500, 2725, 3400, r_ab, 2, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 0, skbuffer, 5425, 2, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 125, skbuffer, 5575, 2, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 250, skbuffer, 5725, 2, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 375, skbuffer, 5875, 2, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 500, skbuffer, 6025, 2, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 625, skbuffer, 6175, 2, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 750, skbuffer, 6325, 2, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 875, skbuffer, 6475, 2, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 1000, skbuffer, 6625, 2, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 0, 2, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDSDD_hpp */
