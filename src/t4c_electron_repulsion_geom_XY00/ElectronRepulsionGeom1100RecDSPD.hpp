#ifndef ElectronRepulsionGeom1100RecDSPD_hpp
#define ElectronRepulsionGeom1100RecDSPD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPSXX.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
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
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DS|1/|r-r'||PD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dspd(T& distributor,
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

    CSimdArray<double> pbuffer(1260, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(864, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(2160, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(6030, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(675, 1);

    // setup Boys fuction data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 8);

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 6, pfactors, 16, bf_data, 6);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 7, pfactors, 16, bf_data, 7);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 8, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 29, 0, 1, 8, 11, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 35, 1, 2, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 41, 2, 3, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 47, 3, 4, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 4, 5, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 5, 6, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 65, 8, 11, 29, 35, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 75, 11, 14, 35, 41, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 14, 17, 41, 47, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 17, 20, 47, 53, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 20, 23, 53, 59, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 115, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 118, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 121, 2, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 130, 3, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 139, 4, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 148, 11, 29, 35, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 166, 14, 35, 41, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 184, 17, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 202, 20, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 220, 35, 65, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 250, 41, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 280, 47, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 310, 53, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 340, 2, 3, 115, 118, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 346, 11, 14, 115, 121, 130, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 364, 14, 17, 118, 130, 139, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 382, 29, 35, 121, 148, 166, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 418, 35, 41, 130, 166, 184, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 454, 41, 47, 139, 184, 202, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 490, 65, 75, 166, 220, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 550, 75, 85, 184, 250, 280, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 610, 85, 95, 202, 280, 310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 670, 121, 130, 340, 346, 364, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 700, 148, 166, 346, 382, 418, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 760, 166, 184, 364, 418, 454, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 820, 220, 250, 418, 490, 550, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 920, 250, 280, 454, 550, 610, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1020, 382, 418, 670, 700, 760, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 1110, 490, 550, 760, 820, 920, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 65, 10, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {148, 166});

                pbuffer.scale(2.0 * b_exp, {220, 250});

                pbuffer.scale(2.0 * b_exp, {382, 418});

                pbuffer.scale(2.0 * b_exp, {490, 550});

                t2cfunc::reduce(cbuffer, 16, pbuffer, 148, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 34, pbuffer, 220, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 64, pbuffer, 382, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 490, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {29, 35});

                pbuffer.scale(2.0 * a_exp, {65, 75});

                pbuffer.scale(a_exp / b_exp, {148, 166});

                pbuffer.scale(a_exp / b_exp, {220, 250});

                pbuffer.scale(a_exp / b_exp, {382, 418});

                pbuffer.scale(a_exp / b_exp, {490, 550});

                t2cfunc::reduce(cbuffer, 160, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 166, pbuffer, 65, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 176, pbuffer, 148, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 194, pbuffer, 220, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 224, pbuffer, 382, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 260, pbuffer, 490, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {148, 166});

                pbuffer.scale(2.0 * b_exp, {220, 250});

                pbuffer.scale(2.0 * b_exp, {382, 418});

                pbuffer.scale(2.0 * b_exp, {490, 550});

                pbuffer.scale(4.0 * a_exp * b_exp, {700, 760});

                pbuffer.scale(4.0 * a_exp * b_exp, {820, 920});

                pbuffer.scale(4.0 * a_exp * b_exp, {1020, 1110});

                pbuffer.scale(4.0 * a_exp * b_exp, {1110, 1260});

                t2cfunc::reduce(cbuffer, 320, pbuffer, 148, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 338, pbuffer, 220, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 368, pbuffer, 382, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 404, pbuffer, 490, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 464, pbuffer, 700, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 524, pbuffer, 820, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 624, pbuffer, 1020, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 714, pbuffer, 1110, 150, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 6, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 180, cbuffer, 16, 34, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 234, cbuffer, 64, 100, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 342, cbuffer, 160, 166, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 360, cbuffer, 176, 194, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 576, cbuffer, 224, 260, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1548, cbuffer, 320, 338, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1602, cbuffer, 368, 404, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1710, cbuffer, 464, 524, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1890, cbuffer, 624, 714, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<1, 2>(skbuffer, 0, ckbuffer, 0, 0, 0);

            t4cfunc::ket_transform<1, 2>(skbuffer, 4065, ckbuffer, 180, 0, 1);

            t4cfunc::ket_transform<1, 2>(skbuffer, 4110, ckbuffer, 234, 0, 2);

            t4cfunc::ket_transform<1, 2>(skbuffer, 4200, ckbuffer, 342, 0, 0);

            t4cfunc::ket_transform<1, 2>(skbuffer, 4215, ckbuffer, 360, 0, 1);

            t4cfunc::ket_transform<1, 2>(skbuffer, 4395, ckbuffer, 576, 0, 2);

            t4cfunc::ket_transform<1, 2>(skbuffer, 5385, ckbuffer, 1548, 0, 1);

            t4cfunc::ket_transform<1, 2>(skbuffer, 5430, ckbuffer, 1602, 0, 2);

            t4cfunc::ket_transform<1, 2>(skbuffer, 5520, ckbuffer, 1710, 0, 3);

            t4cfunc::ket_transform<1, 2>(skbuffer, 5670, ckbuffer, 1890, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 5205, 4200, 4215, r_ab, 1, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 5250, 4215, 4395, r_ab, 1, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 5895, 5385, 5430, r_ab, 1, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_psxx(skbuffer, 1500, 0, 5205, 5250, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 15, 0, 4110, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_psxx(skbuffer, 1365, 0, 4065, 15, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 4260, 4200, 5430, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 4485, 4215, 5520, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 4755, 4395, 5670, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 150, 4215, 4260, 4485, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 555, 4395, 4485, 4755, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_psxx(skbuffer, 1635, 4065, 5205, 5895, 150, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 2040, 15, 5250, 150, 555, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dsxx(skbuffer, 3255, 1365, 1500, 1635, 2040, r_ab, 1, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 0, skbuffer, 3255, 1, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 75, skbuffer, 3345, 1, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 150, skbuffer, 3435, 1, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 225, skbuffer, 3525, 1, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 300, skbuffer, 3615, 1, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 375, skbuffer, 3705, 1, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 450, skbuffer, 3795, 1, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 525, skbuffer, 3885, 1, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 600, skbuffer, 3975, 1, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 0, 1, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDSPD_hpp */
