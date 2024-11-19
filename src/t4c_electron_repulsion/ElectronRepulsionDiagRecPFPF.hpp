#ifndef ElectronRepulsionDiagRecPFPF_hpp
#define ElectronRepulsionDiagRecPFPF_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "BoysFunc.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
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
#include "GtoPairBlock.hpp"
#include "SimdArray.hpp"
#include "T2CUtils.hpp"
#include "T4CUtils.hpp"

namespace erirec {  // erirec namespace

/// @brief Computes (PF|1/|r-r'||PF)  integrals for GTOs pair block.
/// @param distributor The pointer to screening data distributor.
/// @param gto_pair_block The GTOs pair block.
/// @param gto_indices The range [gto_first, gto_last) of GTOs on bra and ket sides.
template <class T>
auto
comp_diag_electron_repulsion_pfpf(T& distributor, const CGtoPairBlock& gto_pair_block, const std::pair<size_t, size_t>& gto_indices) -> void
{
    // intialize GTOs pair data

    const auto a_coords = gto_pair_block.bra_coordinates();

    const auto b_coords = gto_pair_block.ket_coordinates();

    const auto a_vec_exps = gto_pair_block.bra_exponents();

    const auto b_vec_exps = gto_pair_block.ket_exponents();

    const auto ab_vec_norms = gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = gto_pair_block.overlap_factors();

    const auto a_indices = gto_pair_block.bra_orbital_indices();

    const auto b_indices = gto_pair_block.ket_orbital_indices();

    const auto ncgtos = gto_pair_block.number_of_contracted_pairs();

    const auto npgtos = gto_pair_block.number_of_primitive_pairs();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> pfactors(29, npgtos);

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(2060, npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(625, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(750, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(1155, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(441, 1);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, npgtos);

    // allocate aligned array to store max. integral values

    const auto gto_dims = gto_indices.second - gto_indices.first;

    std::vector<double> max_values(gto_dims, 0.0);

    // loop over contracted GTOs on bra and ket sides

    for (auto i = gto_indices.first; i < gto_indices.second; i++)
    {
        // set up indices on ket side

        auto ket_range = std::pair<size_t, size_t>{i, i + 1};

        pfactors.load(a_vec_exps, ket_range, 0, npgtos);

        pfactors.load(b_vec_exps, ket_range, 1, npgtos);

        pfactors.load(ab_vec_ovls, ket_range, 2, npgtos);

        pfactors.load(ab_vec_norms, ket_range, 3, npgtos);

        pfactors.replicate_points(a_coords, ket_range, 4, npgtos);

        pfactors.replicate_points(b_coords, ket_range, 7, npgtos);

        cfactors.replicate_points(a_coords, ket_range, 0, 1);

        cfactors.replicate_points(b_coords, ket_range, 3, 1);

        t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        ckbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // zero integral buffers

        cbuffer.zero();

        ckbuffer.zero();

        skbuffer.zero();

        sbuffer.zero();

        // set up coordinates on bra side

        const auto r_a = a_coords[i];

        const auto r_b = b_coords[i];

        const auto a_xyz = r_a.coordinates();

        const auto b_xyz = r_b.coordinates();

        const auto r_ab = TPoint<double>({a_xyz[0] - b_xyz[0], a_xyz[1] - b_xyz[1], a_xyz[2] - b_xyz[2]});

        for (int j = 0; j < npgtos; j++)
        {
            const auto a_exp = a_vec_exps[j * ncgtos + i];

            const auto b_exp = b_vec_exps[j * ncgtos + i];

            const auto ab_norm = ab_vec_norms[j * ncgtos + i];

            const auto ab_ovl = ab_vec_ovls[j * ncgtos + i];

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

            erirec::comp_prim_electron_repulsion_spss(pbuffer, 210, 3, 4, pfactors, 26, r_pb);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 213, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 222, 4, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 231, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 249, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 267, 21, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 285, 39, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 315, 45, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 345, 51, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 375, 57, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 405, 85, 135, 150, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 450, 95, 150, 165, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 495, 105, 165, 180, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 540, 115, 180, 195, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 585, 15, 18, 210, 213, 222, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 603, 39, 45, 213, 231, 249, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 639, 45, 51, 222, 249, 267, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 675, 75, 85, 231, 285, 315, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 735, 85, 95, 249, 315, 345, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 795, 95, 105, 267, 345, 375, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 855, 135, 150, 315, 405, 450, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 945, 150, 165, 345, 450, 495, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1035, 165, 180, 375, 495, 540, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1125, 231, 249, 585, 603, 639, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1185, 285, 315, 603, 675, 735, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1285, 315, 345, 639, 735, 795, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1385, 405, 450, 735, 855, 945, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1535, 450, 495, 795, 945, 1035, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 1685, 675, 735, 1125, 1185, 1285, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 1835, 855, 945, 1285, 1385, 1535, pfactors, 26, r_pb, a_exp, b_exp);

            t2cfunc::reduce(cbuffer, 0, pbuffer, 1185, 100, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 100, pbuffer, 1385, 150, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 250, pbuffer, 1685, 150, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 400, pbuffer, 1835, 225, ket_width, npgtos);
        }

        erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 100, cfactors, 6, 0, 3);

        erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 300, cbuffer, 250, 400, cfactors, 6, 0, 4);

        t4cfunc::ket_transform<1, 3>(skbuffer, 0, ckbuffer, 0, 0, 3);

        t4cfunc::ket_transform<1, 3>(skbuffer, 210, ckbuffer, 300, 0, 4);

        erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 525, 0, 210, r_ab, 1, 3);

        t4cfunc::bra_transform<1, 3>(sbuffer, 0, skbuffer, 525, 1, 3);

        t4cfunc::update_max_values(max_values, sbuffer, i - gto_indices.first);
    }

    distributor.distribute(max_values, gto_indices);
}

template <class T>
auto
comp_diag_electron_repulsion_fpfp(T& distributor, const CGtoPairBlock& gto_pair_block, const std::pair<size_t, size_t>& gto_indices) -> void
{
    // intialize GTOs pair data

    const auto a_coords = gto_pair_block.ket_coordinates();

    const auto b_coords = gto_pair_block.bra_coordinates();

    const auto a_vec_exps = gto_pair_block.ket_exponents();

    const auto b_vec_exps = gto_pair_block.bra_exponents();

    const auto ab_vec_norms = gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = gto_pair_block.overlap_factors();

    const auto a_indices = gto_pair_block.ket_orbital_indices();

    const auto b_indices = gto_pair_block.bra_orbital_indices();

    const auto ncgtos = gto_pair_block.number_of_contracted_pairs();

    const auto npgtos = gto_pair_block.number_of_primitive_pairs();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> pfactors(29, npgtos);

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(2060, npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(625, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(750, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(1155, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(441, 1);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, npgtos);

    // allocate aligned array to store max. integral values

    const auto gto_dims = gto_indices.second - gto_indices.first;

    std::vector<double> max_values(gto_dims, 0.0);

    // loop over contracted GTOs on bra and ket sides

    for (auto i = gto_indices.first; i < gto_indices.second; i++)
    {
        // set up indices on ket side

        auto ket_range = std::pair<size_t, size_t>{i, i + 1};

        pfactors.load(a_vec_exps, ket_range, 0, npgtos);

        pfactors.load(b_vec_exps, ket_range, 1, npgtos);

        pfactors.load(ab_vec_ovls, ket_range, 2, npgtos);

        pfactors.load(ab_vec_norms, ket_range, 3, npgtos);

        pfactors.replicate_points(a_coords, ket_range, 4, npgtos);

        pfactors.replicate_points(b_coords, ket_range, 7, npgtos);

        cfactors.replicate_points(a_coords, ket_range, 0, 1);

        cfactors.replicate_points(b_coords, ket_range, 3, 1);

        t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        ckbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // zero integral buffers

        cbuffer.zero();

        ckbuffer.zero();

        skbuffer.zero();

        sbuffer.zero();

        // set up coordinates on bra side

        const auto r_a = a_coords[i];

        const auto r_b = b_coords[i];

        const auto a_xyz = r_a.coordinates();

        const auto b_xyz = r_b.coordinates();

        const auto r_ab = TPoint<double>({a_xyz[0] - b_xyz[0], a_xyz[1] - b_xyz[1], a_xyz[2] - b_xyz[2]});

        for (int j = 0; j < npgtos; j++)
        {
            const auto a_exp = a_vec_exps[j * ncgtos + i];

            const auto b_exp = b_vec_exps[j * ncgtos + i];

            const auto ab_norm = ab_vec_norms[j * ncgtos + i];

            const auto ab_ovl = ab_vec_ovls[j * ncgtos + i];

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

            erirec::comp_prim_electron_repulsion_spss(pbuffer, 210, 3, 4, pfactors, 26, r_pb);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 213, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 222, 4, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 231, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 249, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 267, 21, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 285, 39, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 315, 45, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 345, 51, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 375, 57, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 405, 85, 135, 150, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 450, 95, 150, 165, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 495, 105, 165, 180, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 540, 115, 180, 195, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 585, 15, 18, 210, 213, 222, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 603, 39, 45, 213, 231, 249, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 639, 45, 51, 222, 249, 267, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 675, 75, 85, 231, 285, 315, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 735, 85, 95, 249, 315, 345, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 795, 95, 105, 267, 345, 375, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 855, 135, 150, 315, 405, 450, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 945, 150, 165, 345, 450, 495, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1035, 165, 180, 375, 495, 540, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1125, 231, 249, 585, 603, 639, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1185, 285, 315, 603, 675, 735, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1285, 315, 345, 639, 735, 795, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1385, 405, 450, 735, 855, 945, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1535, 450, 495, 795, 945, 1035, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 1685, 675, 735, 1125, 1185, 1285, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 1835, 855, 945, 1285, 1385, 1535, pfactors, 26, r_pb, a_exp, b_exp);

            t2cfunc::reduce(cbuffer, 0, pbuffer, 1185, 100, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 100, pbuffer, 1385, 150, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 250, pbuffer, 1685, 150, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 400, pbuffer, 1835, 225, ket_width, npgtos);
        }

        erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 100, cfactors, 6, 0, 3);

        erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 300, cbuffer, 250, 400, cfactors, 6, 0, 4);

        t4cfunc::ket_transform<1, 3>(skbuffer, 0, ckbuffer, 0, 0, 3);

        t4cfunc::ket_transform<1, 3>(skbuffer, 210, ckbuffer, 300, 0, 4);

        erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 525, 0, 210, r_ab, 1, 3);

        t4cfunc::bra_transform<1, 3>(sbuffer, 0, skbuffer, 525, 1, 3);

        t4cfunc::update_max_values(max_values, sbuffer, i - gto_indices.first);
    }

    distributor.distribute(max_values, gto_indices);
}

}  // namespace erirec

#endif /* ElectronRepulsionDiagRecPFPF_hpp */
