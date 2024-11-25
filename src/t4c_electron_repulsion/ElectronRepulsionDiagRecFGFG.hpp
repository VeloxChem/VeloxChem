#ifndef ElectronRepulsionDiagRecFGFG_hpp
#define ElectronRepulsionDiagRecFGFG_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "BoysFunc.hpp"
#include "ElectronRepulsionContrRecDGXX.hpp"
#include "ElectronRepulsionContrRecDHXX.hpp"
#include "ElectronRepulsionContrRecFGXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecPIXX.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXDH.hpp"
#include "ElectronRepulsionContrRecXXFG.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
#include "ElectronRepulsionContrRecXXPI.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSI.hpp"
#include "ElectronRepulsionPrimRecSDSK.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSI.hpp"
#include "ElectronRepulsionPrimRecSFSK.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSI.hpp"
#include "ElectronRepulsionPrimRecSGSK.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSI.hpp"
#include "ElectronRepulsionPrimRecSHSK.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
#include "ElectronRepulsionPrimRecSISI.hpp"
#include "ElectronRepulsionPrimRecSISK.hpp"
#include "ElectronRepulsionPrimRecSKSG.hpp"
#include "ElectronRepulsionPrimRecSKSH.hpp"
#include "ElectronRepulsionPrimRecSKSI.hpp"
#include "ElectronRepulsionPrimRecSKSK.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
#include "ElectronRepulsionPrimRecSPSI.hpp"
#include "ElectronRepulsionPrimRecSPSK.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSH.hpp"
#include "ElectronRepulsionPrimRecSSSI.hpp"
#include "ElectronRepulsionPrimRecSSSK.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "GtoPairBlock.hpp"
#include "SimdArray.hpp"
#include "T2CUtils.hpp"
#include "T4CUtils.hpp"

namespace erirec {  // erirec namespace

/// @brief Computes (FG|1/|r-r'||FG)  integrals for GTOs pair block.
/// @param distributor The pointer to screening data distributor.
/// @param gto_pair_block The GTOs pair block.
/// @param gto_indices The range [gto_first, gto_last) of GTOs on bra and ket sides.
template <class T>
auto
comp_diag_electron_repulsion_fgfg(T& distributor, const CGtoPairBlock& gto_pair_block, const std::pair<size_t, size_t>& gto_indices) -> void
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

    CSimdArray<double> pbuffer(36346, npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(10000, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(55800, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(41454, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3969, 1);

    // setup Boys fuction data

    const CBoysFunc<14> bf_table;

    CSimdArray<double> bf_data(16, npgtos);

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

            t4cfunc::comp_boys_args(bf_data, 15, pfactors, 13, a_exp, b_exp);

            bf_table.compute(bf_data, 0, 15);

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

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 12, pfactors, 16, bf_data, 12);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 13, pfactors, 16, bf_data, 13);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 14, pfactors, 16, bf_data, 14);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 0, 1, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 1, 2, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 21, 2, 3, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 24, 3, 4, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 27, 4, 5, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 30, 5, 6, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 33, 6, 7, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 36, 7, 8, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 39, 8, 9, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 42, 9, 10, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 45, 10, 11, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 48, 11, 12, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 51, 12, 13, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 54, 13, 14, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 57, 0, 1, 15, 18, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 63, 1, 2, 18, 21, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 69, 2, 3, 21, 24, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 75, 3, 4, 24, 27, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 81, 4, 5, 27, 30, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 87, 5, 6, 30, 33, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 93, 6, 7, 33, 36, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 99, 7, 8, 36, 39, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 105, 8, 9, 39, 42, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 111, 9, 10, 42, 45, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 117, 10, 11, 45, 48, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 123, 11, 12, 48, 51, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 129, 12, 13, 51, 54, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 15, 18, 57, 63, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 18, 21, 63, 69, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 21, 24, 69, 75, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 24, 27, 75, 81, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 175, 27, 30, 81, 87, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 185, 30, 33, 87, 93, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 195, 33, 36, 93, 99, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 205, 36, 39, 99, 105, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 215, 39, 42, 105, 111, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 225, 42, 45, 111, 117, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 235, 45, 48, 117, 123, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 245, 48, 51, 123, 129, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 255, 57, 63, 135, 145, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 270, 63, 69, 145, 155, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 285, 69, 75, 155, 165, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 300, 75, 81, 165, 175, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 315, 81, 87, 175, 185, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 330, 87, 93, 185, 195, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 345, 93, 99, 195, 205, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 360, 99, 105, 205, 215, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 375, 105, 111, 215, 225, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 390, 111, 117, 225, 235, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 405, 117, 123, 235, 245, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 420, 135, 145, 255, 270, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 441, 145, 155, 270, 285, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 462, 155, 165, 285, 300, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 483, 165, 175, 300, 315, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 504, 175, 185, 315, 330, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 525, 185, 195, 330, 345, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 546, 195, 205, 345, 360, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 567, 205, 215, 360, 375, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 588, 215, 225, 375, 390, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 609, 225, 235, 390, 405, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 630, 255, 270, 420, 441, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 658, 270, 285, 441, 462, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 686, 285, 300, 462, 483, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 714, 300, 315, 483, 504, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 742, 315, 330, 504, 525, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 770, 330, 345, 525, 546, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 798, 345, 360, 546, 567, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 826, 360, 375, 567, 588, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 854, 375, 390, 588, 609, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 882, 420, 441, 630, 658, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 918, 441, 462, 658, 686, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 954, 462, 483, 686, 714, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 990, 483, 504, 714, 742, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1026, 504, 525, 742, 770, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1062, 525, 546, 770, 798, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1098, 546, 567, 798, 826, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1134, 567, 588, 826, 854, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spss(pbuffer, 1170, 4, 5, pfactors, 26, r_pb);

            erirec::comp_prim_electron_repulsion_spss(pbuffer, 1173, 5, 6, pfactors, 26, r_pb);

            erirec::comp_prim_electron_repulsion_spss(pbuffer, 1176, 6, 7, pfactors, 26, r_pb);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1179, 4, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1188, 5, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1197, 6, 30, 33, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1206, 7, 33, 36, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1215, 24, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1233, 27, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1251, 30, 81, 87, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1269, 33, 87, 93, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1287, 36, 93, 99, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1305, 69, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1335, 75, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1365, 81, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1395, 87, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1425, 93, 185, 195, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1455, 99, 195, 205, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1485, 145, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1530, 155, 270, 285, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1575, 165, 285, 300, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1620, 175, 300, 315, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1665, 185, 315, 330, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1710, 195, 330, 345, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1755, 205, 345, 360, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1800, 270, 420, 441, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1863, 285, 441, 462, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1926, 300, 462, 483, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1989, 315, 483, 504, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2052, 330, 504, 525, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2115, 345, 525, 546, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2178, 360, 546, 567, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2241, 441, 630, 658, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2325, 462, 658, 686, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2409, 483, 686, 714, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2493, 504, 714, 742, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2577, 525, 742, 770, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2661, 546, 770, 798, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2745, 567, 798, 826, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2829, 658, 882, 918, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2937, 686, 918, 954, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3045, 714, 954, 990, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3153, 742, 990, 1026, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3261, 770, 1026, 1062, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3369, 798, 1062, 1098, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3477, 826, 1098, 1134, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdss(pbuffer, 3585, 4, 5, 1170, 1173, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdss(pbuffer, 3591, 5, 6, 1173, 1176, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3597, 24, 27, 1170, 1179, 1188, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3615, 27, 30, 1173, 1188, 1197, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3633, 30, 33, 1176, 1197, 1206, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3651, 69, 75, 1179, 1215, 1233, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3687, 75, 81, 1188, 1233, 1251, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3723, 81, 87, 1197, 1251, 1269, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3759, 87, 93, 1206, 1269, 1287, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3795, 145, 155, 1215, 1305, 1335, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3855, 155, 165, 1233, 1335, 1365, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3915, 165, 175, 1251, 1365, 1395, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3975, 175, 185, 1269, 1395, 1425, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 4035, 185, 195, 1287, 1425, 1455, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4095, 255, 270, 1305, 1485, 1530, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4185, 270, 285, 1335, 1530, 1575, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4275, 285, 300, 1365, 1575, 1620, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4365, 300, 315, 1395, 1620, 1665, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4455, 315, 330, 1425, 1665, 1710, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4545, 330, 345, 1455, 1710, 1755, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4635, 420, 441, 1530, 1800, 1863, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4761, 441, 462, 1575, 1863, 1926, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4887, 462, 483, 1620, 1926, 1989, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5013, 483, 504, 1665, 1989, 2052, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5139, 504, 525, 1710, 2052, 2115, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5265, 525, 546, 1755, 2115, 2178, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5391, 630, 658, 1863, 2241, 2325, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5559, 658, 686, 1926, 2325, 2409, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5727, 686, 714, 1989, 2409, 2493, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5895, 714, 742, 2052, 2493, 2577, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6063, 742, 770, 2115, 2577, 2661, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6231, 770, 798, 2178, 2661, 2745, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6399, 882, 918, 2325, 2829, 2937, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6615, 918, 954, 2409, 2937, 3045, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6831, 954, 990, 2493, 3045, 3153, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7047, 990, 1026, 2577, 3153, 3261, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7263, 1026, 1062, 2661, 3261, 3369, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7479, 1062, 1098, 2745, 3369, 3477, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfss(pbuffer, 7695, 1170, 1173, 3585, 3591, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 7705, 1179, 1188, 3585, 3597, 3615, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 7735, 1188, 1197, 3591, 3615, 3633, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 7765, 1215, 1233, 3597, 3651, 3687, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 7825, 1233, 1251, 3615, 3687, 3723, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 7885, 1251, 1269, 3633, 3723, 3759, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 7945, 1305, 1335, 3651, 3795, 3855, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 8045, 1335, 1365, 3687, 3855, 3915, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 8145, 1365, 1395, 3723, 3915, 3975, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 8245, 1395, 1425, 3759, 3975, 4035, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8345, 1485, 1530, 3795, 4095, 4185, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8495, 1530, 1575, 3855, 4185, 4275, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8645, 1575, 1620, 3915, 4275, 4365, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8795, 1620, 1665, 3975, 4365, 4455, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8945, 1665, 1710, 4035, 4455, 4545, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9095, 1800, 1863, 4185, 4635, 4761, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9305, 1863, 1926, 4275, 4761, 4887, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9515, 1926, 1989, 4365, 4887, 5013, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9725, 1989, 2052, 4455, 5013, 5139, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9935, 2052, 2115, 4545, 5139, 5265, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10145, 2241, 2325, 4761, 5391, 5559, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10425, 2325, 2409, 4887, 5559, 5727, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10705, 2409, 2493, 5013, 5727, 5895, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10985, 2493, 2577, 5139, 5895, 6063, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 11265, 2577, 2661, 5265, 6063, 6231, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 11545, 2829, 2937, 5559, 6399, 6615, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 11905, 2937, 3045, 5727, 6615, 6831, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 12265, 3045, 3153, 5895, 6831, 7047, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 12625, 3153, 3261, 6063, 7047, 7263, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 12985, 3261, 3369, 6231, 7263, 7479, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 13345, 3597, 3615, 7695, 7705, 7735, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 13390, 3651, 3687, 7705, 7765, 7825, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 13480, 3687, 3723, 7735, 7825, 7885, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 13570, 3795, 3855, 7765, 7945, 8045, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 13720, 3855, 3915, 7825, 8045, 8145, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 13870, 3915, 3975, 7885, 8145, 8245, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 14020, 4095, 4185, 7945, 8345, 8495, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 14245, 4185, 4275, 8045, 8495, 8645, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 14470, 4275, 4365, 8145, 8645, 8795, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 14695, 4365, 4455, 8245, 8795, 8945, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 14920, 4635, 4761, 8495, 9095, 9305, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 15235, 4761, 4887, 8645, 9305, 9515, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 15550, 4887, 5013, 8795, 9515, 9725, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 15865, 5013, 5139, 8945, 9725, 9935, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 16180, 5391, 5559, 9305, 10145, 10425, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 16600, 5559, 5727, 9515, 10425, 10705, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 17020, 5727, 5895, 9725, 10705, 10985, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 17440, 5895, 6063, 9935, 10985, 11265, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 17860, 6399, 6615, 10425, 11545, 11905, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 18400, 6615, 6831, 10705, 11905, 12265, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 18940, 6831, 7047, 10985, 12265, 12625, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 19480, 7047, 7263, 11265, 12625, 12985, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsd(pbuffer, 20020, 7765, 7825, 13345, 13390, 13480, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsf(pbuffer, 20146, 7945, 8045, 13390, 13570, 13720, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsf(pbuffer, 20356, 8045, 8145, 13480, 13720, 13870, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsg(pbuffer, 20566, 8345, 8495, 13570, 14020, 14245, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsg(pbuffer, 20881, 8495, 8645, 13720, 14245, 14470, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsg(pbuffer, 21196, 8645, 8795, 13870, 14470, 14695, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsh(pbuffer, 21511, 9095, 9305, 14245, 14920, 15235, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsh(pbuffer, 21952, 9305, 9515, 14470, 15235, 15550, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsh(pbuffer, 22393, 9515, 9725, 14695, 15550, 15865, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsi(pbuffer, 22834, 10145, 10425, 15235, 16180, 16600, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsi(pbuffer, 23422, 10425, 10705, 15550, 16600, 17020, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsi(pbuffer, 24010, 10705, 10985, 15865, 17020, 17440, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsk(pbuffer, 24598, 11545, 11905, 16600, 17860, 18400, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsk(pbuffer, 25354, 11905, 12265, 17020, 18400, 18940, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsk(pbuffer, 26110, 12265, 12625, 17440, 18940, 19480, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisf(pbuffer, 26866, 13570, 13720, 20020, 20146, 20356, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisg(pbuffer, 27146, 14020, 14245, 20146, 20566, 20881, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisg(pbuffer, 27566, 14245, 14470, 20356, 20881, 21196, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sish(pbuffer, 27986, 14920, 15235, 20881, 21511, 21952, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sish(pbuffer, 28574, 15235, 15550, 21196, 21952, 22393, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisi(pbuffer, 29162, 16180, 16600, 21952, 22834, 23422, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisi(pbuffer, 29946, 16600, 17020, 22393, 23422, 24010, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisk(pbuffer, 30730, 17860, 18400, 23422, 24598, 25354, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisk(pbuffer, 31738, 18400, 18940, 24010, 25354, 26110, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sksg(pbuffer, 32746, 20566, 20881, 26866, 27146, 27566, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sksh(pbuffer, 33286, 21511, 21952, 27566, 27986, 28574, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sksi(pbuffer, 34042, 22834, 23422, 28574, 29162, 29946, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sksk(pbuffer, 35050, 24598, 25354, 29946, 30730, 31738, pfactors, 26, r_pb, a_exp, b_exp);

            t2cfunc::reduce(cbuffer, 0, pbuffer, 14020, 225, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 225, pbuffer, 14920, 315, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 540, pbuffer, 16180, 420, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 960, pbuffer, 17860, 540, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 1500, pbuffer, 20566, 315, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 1815, pbuffer, 21511, 441, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 2256, pbuffer, 22834, 588, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 2844, pbuffer, 24598, 756, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 3600, pbuffer, 27146, 420, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 4020, pbuffer, 27986, 588, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 4608, pbuffer, 29162, 784, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 5392, pbuffer, 30730, 1008, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 6400, pbuffer, 32746, 540, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 6940, pbuffer, 33286, 756, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 7696, pbuffer, 34042, 1008, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 8704, pbuffer, 35050, 1296, ket_width, npgtos);
        }

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 0, cbuffer, 0, 225, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 675, cbuffer, 225, 540, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 1620, cbuffer, 540, 960, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 2880, 0, 675, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 4230, 675, 1620, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 6120, 2880, 4230, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 8370, cbuffer, 1500, 1815, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 9315, cbuffer, 1815, 2256, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 10638, cbuffer, 2256, 2844, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 12402, 8370, 9315, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 14292, 9315, 10638, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 16938, 12402, 14292, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 20088, cbuffer, 3600, 4020, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 21348, cbuffer, 4020, 4608, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 23112, cbuffer, 4608, 5392, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 25464, 20088, 21348, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 27984, 21348, 23112, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 31512, 25464, 27984, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 35712, cbuffer, 6400, 6940, cfactors, 6, 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 37332, cbuffer, 6940, 7696, cfactors, 6, 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 39600, cbuffer, 7696, 8704, cfactors, 6, 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 42624, 35712, 37332, cfactors, 6, 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 45864, 37332, 39600, cfactors, 6, 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 50400, 42624, 45864, cfactors, 6, 0, 7);

        t4cfunc::ket_transform<3, 4>(skbuffer, 0, ckbuffer, 6120, 0, 4);

        t4cfunc::ket_transform<3, 4>(skbuffer, 945, ckbuffer, 16938, 0, 5);

        t4cfunc::ket_transform<3, 4>(skbuffer, 2268, ckbuffer, 31512, 0, 6);

        t4cfunc::ket_transform<3, 4>(skbuffer, 4032, ckbuffer, 50400, 0, 7);

        erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 6300, 0, 945, r_ab, 3, 4);

        erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 9135, 945, 2268, r_ab, 3, 4);

        erirec::comp_bra_hrr_electron_repulsion_pixx(skbuffer, 13104, 2268, 4032, r_ab, 3, 4);

        erirec::comp_bra_hrr_electron_repulsion_dgxx(skbuffer, 18396, 6300, 9135, r_ab, 3, 4);

        erirec::comp_bra_hrr_electron_repulsion_dhxx(skbuffer, 24066, 9135, 13104, r_ab, 3, 4);

        erirec::comp_bra_hrr_electron_repulsion_fgxx(skbuffer, 32004, 18396, 24066, r_ab, 3, 4);

        t4cfunc::bra_transform<3, 4>(sbuffer, 0, skbuffer, 32004, 3, 4);

        t4cfunc::update_max_values(max_values, sbuffer, i - gto_indices.first);
    }

    distributor.distribute(max_values, gto_indices);
}

template <class T>
auto
comp_diag_electron_repulsion_gfgf(T& distributor, const CGtoPairBlock& gto_pair_block, const std::pair<size_t, size_t>& gto_indices) -> void
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

    CSimdArray<double> pbuffer(36346, npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(10000, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(55800, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(41454, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3969, 1);

    // setup Boys fuction data

    const CBoysFunc<14> bf_table;

    CSimdArray<double> bf_data(16, npgtos);

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

            t4cfunc::comp_boys_args(bf_data, 15, pfactors, 13, a_exp, b_exp);

            bf_table.compute(bf_data, 0, 15);

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

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 12, pfactors, 16, bf_data, 12);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 13, pfactors, 16, bf_data, 13);

            erirec::comp_prim_electron_repulsion_ssss(pbuffer, 14, pfactors, 16, bf_data, 14);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 0, 1, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 1, 2, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 21, 2, 3, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 24, 3, 4, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 27, 4, 5, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 30, 5, 6, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 33, 6, 7, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 36, 7, 8, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 39, 8, 9, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 42, 9, 10, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 45, 10, 11, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 48, 11, 12, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 51, 12, 13, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssp(pbuffer, 54, 13, 14, pfactors, 20, 23);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 57, 0, 1, 15, 18, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 63, 1, 2, 18, 21, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 69, 2, 3, 21, 24, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 75, 3, 4, 24, 27, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 81, 4, 5, 27, 30, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 87, 5, 6, 30, 33, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 93, 6, 7, 33, 36, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 99, 7, 8, 36, 39, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 105, 8, 9, 39, 42, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 111, 9, 10, 42, 45, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 117, 10, 11, 45, 48, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 123, 11, 12, 48, 51, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssd(pbuffer, 129, 12, 13, 51, 54, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 15, 18, 57, 63, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 18, 21, 63, 69, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 21, 24, 69, 75, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 24, 27, 75, 81, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 175, 27, 30, 81, 87, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 185, 30, 33, 87, 93, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 195, 33, 36, 93, 99, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 205, 36, 39, 99, 105, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 215, 39, 42, 105, 111, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 225, 42, 45, 111, 117, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 235, 45, 48, 117, 123, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssf(pbuffer, 245, 48, 51, 123, 129, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 255, 57, 63, 135, 145, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 270, 63, 69, 145, 155, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 285, 69, 75, 155, 165, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 300, 75, 81, 165, 175, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 315, 81, 87, 175, 185, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 330, 87, 93, 185, 195, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 345, 93, 99, 195, 205, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 360, 99, 105, 205, 215, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 375, 105, 111, 215, 225, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 390, 111, 117, 225, 235, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssg(pbuffer, 405, 117, 123, 235, 245, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 420, 135, 145, 255, 270, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 441, 145, 155, 270, 285, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 462, 155, 165, 285, 300, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 483, 165, 175, 300, 315, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 504, 175, 185, 315, 330, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 525, 185, 195, 330, 345, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 546, 195, 205, 345, 360, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 567, 205, 215, 360, 375, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 588, 215, 225, 375, 390, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssh(pbuffer, 609, 225, 235, 390, 405, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 630, 255, 270, 420, 441, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 658, 270, 285, 441, 462, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 686, 285, 300, 462, 483, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 714, 300, 315, 483, 504, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 742, 315, 330, 504, 525, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 770, 330, 345, 525, 546, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 798, 345, 360, 546, 567, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 826, 360, 375, 567, 588, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssi(pbuffer, 854, 375, 390, 588, 609, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 882, 420, 441, 630, 658, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 918, 441, 462, 658, 686, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 954, 462, 483, 686, 714, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 990, 483, 504, 714, 742, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1026, 504, 525, 742, 770, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1062, 525, 546, 770, 798, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1098, 546, 567, 798, 826, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1134, 567, 588, 826, 854, pfactors, 20, 23, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spss(pbuffer, 1170, 4, 5, pfactors, 26, r_pb);

            erirec::comp_prim_electron_repulsion_spss(pbuffer, 1173, 5, 6, pfactors, 26, r_pb);

            erirec::comp_prim_electron_repulsion_spss(pbuffer, 1176, 6, 7, pfactors, 26, r_pb);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1179, 4, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1188, 5, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1197, 6, 30, 33, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1206, 7, 33, 36, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1215, 24, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1233, 27, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1251, 30, 81, 87, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1269, 33, 87, 93, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1287, 36, 93, 99, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1305, 69, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1335, 75, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1365, 81, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1395, 87, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1425, 93, 185, 195, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1455, 99, 195, 205, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1485, 145, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1530, 155, 270, 285, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1575, 165, 285, 300, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1620, 175, 300, 315, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1665, 185, 315, 330, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1710, 195, 330, 345, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1755, 205, 345, 360, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1800, 270, 420, 441, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1863, 285, 441, 462, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1926, 300, 462, 483, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1989, 315, 483, 504, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2052, 330, 504, 525, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2115, 345, 525, 546, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2178, 360, 546, 567, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2241, 441, 630, 658, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2325, 462, 658, 686, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2409, 483, 686, 714, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2493, 504, 714, 742, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2577, 525, 742, 770, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2661, 546, 770, 798, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2745, 567, 798, 826, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2829, 658, 882, 918, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2937, 686, 918, 954, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3045, 714, 954, 990, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3153, 742, 990, 1026, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3261, 770, 1026, 1062, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3369, 798, 1062, 1098, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3477, 826, 1098, 1134, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdss(pbuffer, 3585, 4, 5, 1170, 1173, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdss(pbuffer, 3591, 5, 6, 1173, 1176, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3597, 24, 27, 1170, 1179, 1188, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3615, 27, 30, 1173, 1188, 1197, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3633, 30, 33, 1176, 1197, 1206, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3651, 69, 75, 1179, 1215, 1233, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3687, 75, 81, 1188, 1233, 1251, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3723, 81, 87, 1197, 1251, 1269, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3759, 87, 93, 1206, 1269, 1287, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3795, 145, 155, 1215, 1305, 1335, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3855, 155, 165, 1233, 1335, 1365, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3915, 165, 175, 1251, 1365, 1395, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3975, 175, 185, 1269, 1395, 1425, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 4035, 185, 195, 1287, 1425, 1455, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4095, 255, 270, 1305, 1485, 1530, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4185, 270, 285, 1335, 1530, 1575, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4275, 285, 300, 1365, 1575, 1620, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4365, 300, 315, 1395, 1620, 1665, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4455, 315, 330, 1425, 1665, 1710, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4545, 330, 345, 1455, 1710, 1755, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4635, 420, 441, 1530, 1800, 1863, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4761, 441, 462, 1575, 1863, 1926, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4887, 462, 483, 1620, 1926, 1989, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5013, 483, 504, 1665, 1989, 2052, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5139, 504, 525, 1710, 2052, 2115, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5265, 525, 546, 1755, 2115, 2178, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5391, 630, 658, 1863, 2241, 2325, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5559, 658, 686, 1926, 2325, 2409, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5727, 686, 714, 1989, 2409, 2493, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5895, 714, 742, 2052, 2493, 2577, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6063, 742, 770, 2115, 2577, 2661, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6231, 770, 798, 2178, 2661, 2745, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6399, 882, 918, 2325, 2829, 2937, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6615, 918, 954, 2409, 2937, 3045, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6831, 954, 990, 2493, 3045, 3153, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7047, 990, 1026, 2577, 3153, 3261, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7263, 1026, 1062, 2661, 3261, 3369, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7479, 1062, 1098, 2745, 3369, 3477, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfss(pbuffer, 7695, 1170, 1173, 3585, 3591, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 7705, 1179, 1188, 3585, 3597, 3615, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 7735, 1188, 1197, 3591, 3615, 3633, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 7765, 1215, 1233, 3597, 3651, 3687, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 7825, 1233, 1251, 3615, 3687, 3723, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 7885, 1251, 1269, 3633, 3723, 3759, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 7945, 1305, 1335, 3651, 3795, 3855, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 8045, 1335, 1365, 3687, 3855, 3915, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 8145, 1365, 1395, 3723, 3915, 3975, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 8245, 1395, 1425, 3759, 3975, 4035, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8345, 1485, 1530, 3795, 4095, 4185, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8495, 1530, 1575, 3855, 4185, 4275, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8645, 1575, 1620, 3915, 4275, 4365, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8795, 1620, 1665, 3975, 4365, 4455, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8945, 1665, 1710, 4035, 4455, 4545, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9095, 1800, 1863, 4185, 4635, 4761, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9305, 1863, 1926, 4275, 4761, 4887, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9515, 1926, 1989, 4365, 4887, 5013, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9725, 1989, 2052, 4455, 5013, 5139, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9935, 2052, 2115, 4545, 5139, 5265, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10145, 2241, 2325, 4761, 5391, 5559, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10425, 2325, 2409, 4887, 5559, 5727, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10705, 2409, 2493, 5013, 5727, 5895, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10985, 2493, 2577, 5139, 5895, 6063, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 11265, 2577, 2661, 5265, 6063, 6231, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 11545, 2829, 2937, 5559, 6399, 6615, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 11905, 2937, 3045, 5727, 6615, 6831, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 12265, 3045, 3153, 5895, 6831, 7047, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 12625, 3153, 3261, 6063, 7047, 7263, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 12985, 3261, 3369, 6231, 7263, 7479, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 13345, 3597, 3615, 7695, 7705, 7735, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 13390, 3651, 3687, 7705, 7765, 7825, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 13480, 3687, 3723, 7735, 7825, 7885, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 13570, 3795, 3855, 7765, 7945, 8045, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 13720, 3855, 3915, 7825, 8045, 8145, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 13870, 3915, 3975, 7885, 8145, 8245, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 14020, 4095, 4185, 7945, 8345, 8495, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 14245, 4185, 4275, 8045, 8495, 8645, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 14470, 4275, 4365, 8145, 8645, 8795, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 14695, 4365, 4455, 8245, 8795, 8945, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 14920, 4635, 4761, 8495, 9095, 9305, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 15235, 4761, 4887, 8645, 9305, 9515, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 15550, 4887, 5013, 8795, 9515, 9725, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 15865, 5013, 5139, 8945, 9725, 9935, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 16180, 5391, 5559, 9305, 10145, 10425, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 16600, 5559, 5727, 9515, 10425, 10705, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 17020, 5727, 5895, 9725, 10705, 10985, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 17440, 5895, 6063, 9935, 10985, 11265, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 17860, 6399, 6615, 10425, 11545, 11905, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 18400, 6615, 6831, 10705, 11905, 12265, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 18940, 6831, 7047, 10985, 12265, 12625, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 19480, 7047, 7263, 11265, 12625, 12985, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsd(pbuffer, 20020, 7765, 7825, 13345, 13390, 13480, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsf(pbuffer, 20146, 7945, 8045, 13390, 13570, 13720, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsf(pbuffer, 20356, 8045, 8145, 13480, 13720, 13870, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsg(pbuffer, 20566, 8345, 8495, 13570, 14020, 14245, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsg(pbuffer, 20881, 8495, 8645, 13720, 14245, 14470, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsg(pbuffer, 21196, 8645, 8795, 13870, 14470, 14695, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsh(pbuffer, 21511, 9095, 9305, 14245, 14920, 15235, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsh(pbuffer, 21952, 9305, 9515, 14470, 15235, 15550, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsh(pbuffer, 22393, 9515, 9725, 14695, 15550, 15865, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsi(pbuffer, 22834, 10145, 10425, 15235, 16180, 16600, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsi(pbuffer, 23422, 10425, 10705, 15550, 16600, 17020, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsi(pbuffer, 24010, 10705, 10985, 15865, 17020, 17440, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsk(pbuffer, 24598, 11545, 11905, 16600, 17860, 18400, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsk(pbuffer, 25354, 11905, 12265, 17020, 18400, 18940, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsk(pbuffer, 26110, 12265, 12625, 17440, 18940, 19480, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisf(pbuffer, 26866, 13570, 13720, 20020, 20146, 20356, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisg(pbuffer, 27146, 14020, 14245, 20146, 20566, 20881, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisg(pbuffer, 27566, 14245, 14470, 20356, 20881, 21196, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sish(pbuffer, 27986, 14920, 15235, 20881, 21511, 21952, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sish(pbuffer, 28574, 15235, 15550, 21196, 21952, 22393, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisi(pbuffer, 29162, 16180, 16600, 21952, 22834, 23422, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisi(pbuffer, 29946, 16600, 17020, 22393, 23422, 24010, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisk(pbuffer, 30730, 17860, 18400, 23422, 24598, 25354, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sisk(pbuffer, 31738, 18400, 18940, 24010, 25354, 26110, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sksg(pbuffer, 32746, 20566, 20881, 26866, 27146, 27566, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sksh(pbuffer, 33286, 21511, 21952, 27566, 27986, 28574, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sksi(pbuffer, 34042, 22834, 23422, 28574, 29162, 29946, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sksk(pbuffer, 35050, 24598, 25354, 29946, 30730, 31738, pfactors, 26, r_pb, a_exp, b_exp);

            t2cfunc::reduce(cbuffer, 0, pbuffer, 14020, 225, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 225, pbuffer, 14920, 315, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 540, pbuffer, 16180, 420, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 960, pbuffer, 17860, 540, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 1500, pbuffer, 20566, 315, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 1815, pbuffer, 21511, 441, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 2256, pbuffer, 22834, 588, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 2844, pbuffer, 24598, 756, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 3600, pbuffer, 27146, 420, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 4020, pbuffer, 27986, 588, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 4608, pbuffer, 29162, 784, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 5392, pbuffer, 30730, 1008, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 6400, pbuffer, 32746, 540, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 6940, pbuffer, 33286, 756, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 7696, pbuffer, 34042, 1008, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 8704, pbuffer, 35050, 1296, ket_width, npgtos);
        }

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 0, cbuffer, 0, 225, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 675, cbuffer, 225, 540, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 1620, cbuffer, 540, 960, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 2880, 0, 675, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 4230, 675, 1620, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 6120, 2880, 4230, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 8370, cbuffer, 1500, 1815, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 9315, cbuffer, 1815, 2256, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 10638, cbuffer, 2256, 2844, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 12402, 8370, 9315, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 14292, 9315, 10638, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 16938, 12402, 14292, cfactors, 6, 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 20088, cbuffer, 3600, 4020, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 21348, cbuffer, 4020, 4608, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 23112, cbuffer, 4608, 5392, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 25464, 20088, 21348, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 27984, 21348, 23112, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 31512, 25464, 27984, cfactors, 6, 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 35712, cbuffer, 6400, 6940, cfactors, 6, 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 37332, cbuffer, 6940, 7696, cfactors, 6, 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 39600, cbuffer, 7696, 8704, cfactors, 6, 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 42624, 35712, 37332, cfactors, 6, 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 45864, 37332, 39600, cfactors, 6, 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 50400, 42624, 45864, cfactors, 6, 0, 7);

        t4cfunc::ket_transform<3, 4>(skbuffer, 0, ckbuffer, 6120, 0, 4);

        t4cfunc::ket_transform<3, 4>(skbuffer, 945, ckbuffer, 16938, 0, 5);

        t4cfunc::ket_transform<3, 4>(skbuffer, 2268, ckbuffer, 31512, 0, 6);

        t4cfunc::ket_transform<3, 4>(skbuffer, 4032, ckbuffer, 50400, 0, 7);

        erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 6300, 0, 945, r_ab, 3, 4);

        erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 9135, 945, 2268, r_ab, 3, 4);

        erirec::comp_bra_hrr_electron_repulsion_pixx(skbuffer, 13104, 2268, 4032, r_ab, 3, 4);

        erirec::comp_bra_hrr_electron_repulsion_dgxx(skbuffer, 18396, 6300, 9135, r_ab, 3, 4);

        erirec::comp_bra_hrr_electron_repulsion_dhxx(skbuffer, 24066, 9135, 13104, r_ab, 3, 4);

        erirec::comp_bra_hrr_electron_repulsion_fgxx(skbuffer, 32004, 18396, 24066, r_ab, 3, 4);

        t4cfunc::bra_transform<3, 4>(sbuffer, 0, skbuffer, 32004, 3, 4);

        t4cfunc::update_max_values(max_values, sbuffer, i - gto_indices.first);
    }

    distributor.distribute(max_values, gto_indices);
}

}  // namespace erirec

#endif /* ElectronRepulsionDiagRecFGFG_hpp */
