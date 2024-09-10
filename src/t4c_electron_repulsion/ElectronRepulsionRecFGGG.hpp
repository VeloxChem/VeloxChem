#ifndef ElectronRepulsionRecFGGG_hpp
#define ElectronRepulsionRecFGGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDGXX.hpp"
#include "ElectronRepulsionContrRecDHXX.hpp"
#include "ElectronRepulsionContrRecFGXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecPIXX.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXDH.hpp"
#include "ElectronRepulsionContrRecXXDI.hpp"
#include "ElectronRepulsionContrRecXXFG.hpp"
#include "ElectronRepulsionContrRecXXFH.hpp"
#include "ElectronRepulsionContrRecXXGG.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
#include "ElectronRepulsionContrRecXXPI.hpp"
#include "ElectronRepulsionContrRecXXPK.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSI.hpp"
#include "ElectronRepulsionPrimRecSDSK.hpp"
#include "ElectronRepulsionPrimRecSDSL.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSI.hpp"
#include "ElectronRepulsionPrimRecSFSK.hpp"
#include "ElectronRepulsionPrimRecSFSL.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSI.hpp"
#include "ElectronRepulsionPrimRecSGSK.hpp"
#include "ElectronRepulsionPrimRecSGSL.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSI.hpp"
#include "ElectronRepulsionPrimRecSHSK.hpp"
#include "ElectronRepulsionPrimRecSHSL.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
#include "ElectronRepulsionPrimRecSISI.hpp"
#include "ElectronRepulsionPrimRecSISK.hpp"
#include "ElectronRepulsionPrimRecSISL.hpp"
#include "ElectronRepulsionPrimRecSKSG.hpp"
#include "ElectronRepulsionPrimRecSKSH.hpp"
#include "ElectronRepulsionPrimRecSKSI.hpp"
#include "ElectronRepulsionPrimRecSKSK.hpp"
#include "ElectronRepulsionPrimRecSKSL.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
#include "ElectronRepulsionPrimRecSPSI.hpp"
#include "ElectronRepulsionPrimRecSPSK.hpp"
#include "ElectronRepulsionPrimRecSPSL.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSH.hpp"
#include "ElectronRepulsionPrimRecSSSI.hpp"
#include "ElectronRepulsionPrimRecSSSK.hpp"
#include "ElectronRepulsionPrimRecSSSL.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "GtoPairBlock.hpp"
#include "BatchFunc.hpp"

namespace erirec { // erirec namespace

/// @brief Computes (FG|1/|r-r'||GG)  integrals for two basis function pairs blocks.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
/// @param bra_eq_ket True if basis function pairs blocks on bra and ket are the same, False otherwise.
template <class T>
inline auto
comp_electron_repulsion_fggg(T& distributor,
                             const CGtoPairBlock& bra_gto_pair_block,
                             const CGtoPairBlock& ket_gto_pair_block,
                             const std::pair<size_t, size_t>& bra_indices,
                             const std::pair<size_t, size_t>& ket_indices,
                             const bool bra_eq_ket) -> void
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

    CSimdArray<double> pbuffer(51316, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(14500, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(126900, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(53298, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(5103, 1);

    // setup Boys fuction data

    const CBoysFunc<15> bf_table;

    CSimdArray<double> bf_data(17, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 16, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 16);

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

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 15, pfactors, 16, bf_data, 15);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 25, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 28, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 31, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 34, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 37, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 40, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 43, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 46, 10, 11, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 49, 11, 12, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 52, 12, 13, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 55, 13, 14, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 58, 14, 15, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 61, 0, 1, 16, 19, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 67, 1, 2, 19, 22, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 73, 2, 3, 22, 25, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 79, 3, 4, 25, 28, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 85, 4, 5, 28, 31, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 91, 5, 6, 31, 34, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 97, 6, 7, 34, 37, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 103, 7, 8, 37, 40, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 109, 8, 9, 40, 43, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 115, 9, 10, 43, 46, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 121, 10, 11, 46, 49, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 127, 11, 12, 49, 52, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 133, 12, 13, 52, 55, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 139, 13, 14, 55, 58, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 16, 19, 61, 67, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 19, 22, 67, 73, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 22, 25, 73, 79, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 175, 25, 28, 79, 85, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 185, 28, 31, 85, 91, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 195, 31, 34, 91, 97, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 205, 34, 37, 97, 103, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 215, 37, 40, 103, 109, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 225, 40, 43, 109, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 235, 43, 46, 115, 121, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 245, 46, 49, 121, 127, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 255, 49, 52, 127, 133, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 265, 52, 55, 133, 139, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 275, 61, 67, 145, 155, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 290, 67, 73, 155, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 305, 73, 79, 165, 175, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 320, 79, 85, 175, 185, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 335, 85, 91, 185, 195, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 350, 91, 97, 195, 205, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 365, 97, 103, 205, 215, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 380, 103, 109, 215, 225, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 395, 109, 115, 225, 235, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 410, 115, 121, 235, 245, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 425, 121, 127, 245, 255, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 440, 127, 133, 255, 265, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 455, 145, 155, 275, 290, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 476, 155, 165, 290, 305, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 497, 165, 175, 305, 320, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 518, 175, 185, 320, 335, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 539, 185, 195, 335, 350, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 560, 195, 205, 350, 365, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 581, 205, 215, 365, 380, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 602, 215, 225, 380, 395, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 623, 225, 235, 395, 410, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 644, 235, 245, 410, 425, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 665, 245, 255, 425, 440, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 686, 275, 290, 455, 476, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 714, 290, 305, 476, 497, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 742, 305, 320, 497, 518, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 770, 320, 335, 518, 539, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 798, 335, 350, 539, 560, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 826, 350, 365, 560, 581, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 854, 365, 380, 581, 602, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 882, 380, 395, 602, 623, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 910, 395, 410, 623, 644, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 938, 410, 425, 644, 665, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 966, 455, 476, 686, 714, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1002, 476, 497, 714, 742, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1038, 497, 518, 742, 770, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1074, 518, 539, 770, 798, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1110, 539, 560, 798, 826, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1146, 560, 581, 826, 854, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1182, 581, 602, 854, 882, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1218, 602, 623, 882, 910, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1254, 623, 644, 910, 938, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1290, 686, 714, 966, 1002, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1335, 714, 742, 1002, 1038, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1380, 742, 770, 1038, 1074, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1425, 770, 798, 1074, 1110, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1470, 798, 826, 1110, 1146, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1515, 826, 854, 1146, 1182, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1560, 854, 882, 1182, 1218, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1605, 882, 910, 1218, 1254, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1650, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1653, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1656, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1659, 4, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1668, 5, 28, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1677, 6, 31, 34, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1686, 7, 34, 37, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1695, 25, 73, 79, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1713, 28, 79, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1731, 31, 85, 91, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1749, 34, 91, 97, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1767, 37, 97, 103, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1785, 73, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1815, 79, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1845, 85, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1875, 91, 185, 195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1905, 97, 195, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1935, 103, 205, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1965, 155, 275, 290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 2010, 165, 290, 305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 2055, 175, 305, 320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 2100, 185, 320, 335, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 2145, 195, 335, 350, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 2190, 205, 350, 365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 2235, 215, 365, 380, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2280, 290, 455, 476, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2343, 305, 476, 497, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2406, 320, 497, 518, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2469, 335, 518, 539, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2532, 350, 539, 560, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2595, 365, 560, 581, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2658, 380, 581, 602, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2721, 476, 686, 714, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2805, 497, 714, 742, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2889, 518, 742, 770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2973, 539, 770, 798, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 3057, 560, 798, 826, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 3141, 581, 826, 854, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 3225, 602, 854, 882, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3309, 714, 966, 1002, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3417, 742, 1002, 1038, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3525, 770, 1038, 1074, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3633, 798, 1074, 1110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3741, 826, 1110, 1146, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3849, 854, 1146, 1182, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3957, 882, 1182, 1218, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 4065, 1002, 1290, 1335, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 4200, 1038, 1335, 1380, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 4335, 1074, 1380, 1425, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 4470, 1110, 1425, 1470, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 4605, 1146, 1470, 1515, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 4740, 1182, 1515, 1560, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 4875, 1218, 1560, 1605, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 5010, 4, 5, 1650, 1653, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 5016, 5, 6, 1653, 1656, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 5022, 25, 28, 1650, 1659, 1668, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 5040, 28, 31, 1653, 1668, 1677, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 5058, 31, 34, 1656, 1677, 1686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 5076, 73, 79, 1659, 1695, 1713, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 5112, 79, 85, 1668, 1713, 1731, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 5148, 85, 91, 1677, 1731, 1749, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 5184, 91, 97, 1686, 1749, 1767, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 5220, 155, 165, 1695, 1785, 1815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 5280, 165, 175, 1713, 1815, 1845, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 5340, 175, 185, 1731, 1845, 1875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 5400, 185, 195, 1749, 1875, 1905, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 5460, 195, 205, 1767, 1905, 1935, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 5520, 275, 290, 1785, 1965, 2010, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 5610, 290, 305, 1815, 2010, 2055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 5700, 305, 320, 1845, 2055, 2100, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 5790, 320, 335, 1875, 2100, 2145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 5880, 335, 350, 1905, 2145, 2190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 5970, 350, 365, 1935, 2190, 2235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 6060, 455, 476, 2010, 2280, 2343, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 6186, 476, 497, 2055, 2343, 2406, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 6312, 497, 518, 2100, 2406, 2469, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 6438, 518, 539, 2145, 2469, 2532, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 6564, 539, 560, 2190, 2532, 2595, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 6690, 560, 581, 2235, 2595, 2658, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6816, 686, 714, 2343, 2721, 2805, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6984, 714, 742, 2406, 2805, 2889, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 7152, 742, 770, 2469, 2889, 2973, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 7320, 770, 798, 2532, 2973, 3057, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 7488, 798, 826, 2595, 3057, 3141, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 7656, 826, 854, 2658, 3141, 3225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7824, 966, 1002, 2805, 3309, 3417, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 8040, 1002, 1038, 2889, 3417, 3525, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 8256, 1038, 1074, 2973, 3525, 3633, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 8472, 1074, 1110, 3057, 3633, 3741, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 8688, 1110, 1146, 3141, 3741, 3849, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 8904, 1146, 1182, 3225, 3849, 3957, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsl(pbuffer, 9120, 1290, 1335, 3417, 4065, 4200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsl(pbuffer, 9390, 1335, 1380, 3525, 4200, 4335, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsl(pbuffer, 9660, 1380, 1425, 3633, 4335, 4470, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsl(pbuffer, 9930, 1425, 1470, 3741, 4470, 4605, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsl(pbuffer, 10200, 1470, 1515, 3849, 4605, 4740, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsl(pbuffer, 10470, 1515, 1560, 3957, 4740, 4875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 10740, 1650, 1653, 5010, 5016, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 10750, 1659, 1668, 5010, 5022, 5040, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 10780, 1668, 1677, 5016, 5040, 5058, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 10810, 1695, 1713, 5022, 5076, 5112, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 10870, 1713, 1731, 5040, 5112, 5148, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 10930, 1731, 1749, 5058, 5148, 5184, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 10990, 1785, 1815, 5076, 5220, 5280, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 11090, 1815, 1845, 5112, 5280, 5340, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 11190, 1845, 1875, 5148, 5340, 5400, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 11290, 1875, 1905, 5184, 5400, 5460, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 11390, 1965, 2010, 5220, 5520, 5610, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 11540, 2010, 2055, 5280, 5610, 5700, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 11690, 2055, 2100, 5340, 5700, 5790, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 11840, 2100, 2145, 5400, 5790, 5880, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 11990, 2145, 2190, 5460, 5880, 5970, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 12140, 2280, 2343, 5610, 6060, 6186, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 12350, 2343, 2406, 5700, 6186, 6312, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 12560, 2406, 2469, 5790, 6312, 6438, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 12770, 2469, 2532, 5880, 6438, 6564, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 12980, 2532, 2595, 5970, 6564, 6690, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 13190, 2721, 2805, 6186, 6816, 6984, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 13470, 2805, 2889, 6312, 6984, 7152, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 13750, 2889, 2973, 6438, 7152, 7320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 14030, 2973, 3057, 6564, 7320, 7488, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 14310, 3057, 3141, 6690, 7488, 7656, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 14590, 3309, 3417, 6984, 7824, 8040, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 14950, 3417, 3525, 7152, 8040, 8256, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 15310, 3525, 3633, 7320, 8256, 8472, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 15670, 3633, 3741, 7488, 8472, 8688, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 16030, 3741, 3849, 7656, 8688, 8904, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsl(pbuffer, 16390, 4065, 4200, 8040, 9120, 9390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsl(pbuffer, 16840, 4200, 4335, 8256, 9390, 9660, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsl(pbuffer, 17290, 4335, 4470, 8472, 9660, 9930, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsl(pbuffer, 17740, 4470, 4605, 8688, 9930, 10200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsl(pbuffer, 18190, 4605, 4740, 8904, 10200, 10470, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 18640, 5022, 5040, 10740, 10750, 10780, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 18685, 5076, 5112, 10750, 10810, 10870, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 18775, 5112, 5148, 10780, 10870, 10930, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 18865, 5220, 5280, 10810, 10990, 11090, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 19015, 5280, 5340, 10870, 11090, 11190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 19165, 5340, 5400, 10930, 11190, 11290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 19315, 5520, 5610, 10990, 11390, 11540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 19540, 5610, 5700, 11090, 11540, 11690, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 19765, 5700, 5790, 11190, 11690, 11840, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 19990, 5790, 5880, 11290, 11840, 11990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 20215, 6060, 6186, 11540, 12140, 12350, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 20530, 6186, 6312, 11690, 12350, 12560, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 20845, 6312, 6438, 11840, 12560, 12770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 21160, 6438, 6564, 11990, 12770, 12980, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 21475, 6816, 6984, 12350, 13190, 13470, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 21895, 6984, 7152, 12560, 13470, 13750, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 22315, 7152, 7320, 12770, 13750, 14030, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 22735, 7320, 7488, 12980, 14030, 14310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 23155, 7824, 8040, 13470, 14590, 14950, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 23695, 8040, 8256, 13750, 14950, 15310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 24235, 8256, 8472, 14030, 15310, 15670, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 24775, 8472, 8688, 14310, 15670, 16030, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsl(pbuffer, 25315, 9120, 9390, 14950, 16390, 16840, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsl(pbuffer, 25990, 9390, 9660, 15310, 16840, 17290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsl(pbuffer, 26665, 9660, 9930, 15670, 17290, 17740, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsl(pbuffer, 27340, 9930, 10200, 16030, 17740, 18190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 28015, 10810, 10870, 18640, 18685, 18775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 28141, 10990, 11090, 18685, 18865, 19015, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 28351, 11090, 11190, 18775, 19015, 19165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 28561, 11390, 11540, 18865, 19315, 19540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 28876, 11540, 11690, 19015, 19540, 19765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 29191, 11690, 11840, 19165, 19765, 19990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 29506, 12140, 12350, 19540, 20215, 20530, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 29947, 12350, 12560, 19765, 20530, 20845, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 30388, 12560, 12770, 19990, 20845, 21160, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 30829, 13190, 13470, 20530, 21475, 21895, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 31417, 13470, 13750, 20845, 21895, 22315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 32005, 13750, 14030, 21160, 22315, 22735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsk(pbuffer, 32593, 14590, 14950, 21895, 23155, 23695, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsk(pbuffer, 33349, 14950, 15310, 22315, 23695, 24235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsk(pbuffer, 34105, 15310, 15670, 22735, 24235, 24775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsl(pbuffer, 34861, 16390, 16840, 23695, 25315, 25990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsl(pbuffer, 35806, 16840, 17290, 24235, 25990, 26665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsl(pbuffer, 36751, 17290, 17740, 24775, 26665, 27340, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 37696, 18865, 19015, 28015, 28141, 28351, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 37976, 19315, 19540, 28141, 28561, 28876, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 38396, 19540, 19765, 28351, 28876, 29191, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 38816, 20215, 20530, 28876, 29506, 29947, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 39404, 20530, 20845, 29191, 29947, 30388, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 39992, 21475, 21895, 29947, 30829, 31417, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 40776, 21895, 22315, 30388, 31417, 32005, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisk(pbuffer, 41560, 23155, 23695, 31417, 32593, 33349, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisk(pbuffer, 42568, 23695, 24235, 32005, 33349, 34105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisl(pbuffer, 43576, 25315, 25990, 33349, 34861, 35806, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisl(pbuffer, 44836, 25990, 26665, 34105, 35806, 36751, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 46096, 28561, 28876, 37696, 37976, 38396, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksh(pbuffer, 46636, 29506, 29947, 38396, 38816, 39404, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksi(pbuffer, 47392, 30829, 31417, 39404, 39992, 40776, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksk(pbuffer, 48400, 32593, 33349, 40776, 41560, 42568, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksl(pbuffer, 49696, 34861, 35806, 42568, 43576, 44836, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 19315, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 225, pbuffer, 20215, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 540, pbuffer, 21475, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 960, pbuffer, 23155, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1500, pbuffer, 25315, 675, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2175, pbuffer, 28561, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2490, pbuffer, 29506, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2931, pbuffer, 30829, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3519, pbuffer, 32593, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4275, pbuffer, 34861, 945, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5220, pbuffer, 37976, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5640, pbuffer, 38816, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6228, pbuffer, 39992, 784, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7012, pbuffer, 41560, 1008, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8020, pbuffer, 43576, 1260, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9280, pbuffer, 46096, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9820, pbuffer, 46636, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10576, pbuffer, 47392, 1008, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11584, pbuffer, 48400, 1296, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12880, pbuffer, 49696, 1620, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 0, cbuffer, 0, 225, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 675, cbuffer, 225, 540, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 1620, cbuffer, 540, 960, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpk(ckbuffer, 2880, cbuffer, 960, 1500, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 4500, 0, 675, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 5850, 675, 1620, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdi(ckbuffer, 7740, 1620, 2880, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 10260, 4500, 5850, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxfh(ckbuffer, 12510, 5850, 7740, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxgg(ckbuffer, 15660, 10260, 12510, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 19035, cbuffer, 2175, 2490, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 19980, cbuffer, 2490, 2931, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 21303, cbuffer, 2931, 3519, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpk(ckbuffer, 23067, cbuffer, 3519, 4275, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 25335, 19035, 19980, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 27225, 19980, 21303, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdi(ckbuffer, 29871, 21303, 23067, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 33399, 25335, 27225, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxfh(ckbuffer, 36549, 27225, 29871, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxgg(ckbuffer, 40959, 33399, 36549, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 45684, cbuffer, 5220, 5640, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 46944, cbuffer, 5640, 6228, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 48708, cbuffer, 6228, 7012, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpk(ckbuffer, 51060, cbuffer, 7012, 8020, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 54084, 45684, 46944, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 56604, 46944, 48708, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdi(ckbuffer, 60132, 48708, 51060, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 64836, 54084, 56604, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxfh(ckbuffer, 69036, 56604, 60132, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxgg(ckbuffer, 74916, 64836, 69036, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 81216, cbuffer, 9280, 9820, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 82836, cbuffer, 9820, 10576, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 85104, cbuffer, 10576, 11584, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpk(ckbuffer, 88128, cbuffer, 11584, 12880, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 92016, 81216, 82836, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 95256, 82836, 85104, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdi(ckbuffer, 99792, 85104, 88128, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 105840, 92016, 95256, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxfh(ckbuffer, 111240, 95256, 99792, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxgg(ckbuffer, 118800, 105840, 111240, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<4, 4>(skbuffer, 0, ckbuffer, 15660, 0, 4);

            t4cfunc::ket_transform<4, 4>(skbuffer, 1215, ckbuffer, 40959, 0, 5);

            t4cfunc::ket_transform<4, 4>(skbuffer, 2916, ckbuffer, 74916, 0, 6);

            t4cfunc::ket_transform<4, 4>(skbuffer, 5184, ckbuffer, 118800, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 8100, 0, 1215, r_ab, 4, 4);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 11745, 1215, 2916, r_ab, 4, 4);

            erirec::comp_bra_hrr_electron_repulsion_pixx(skbuffer, 16848, 2916, 5184, r_ab, 4, 4);

            erirec::comp_bra_hrr_electron_repulsion_dgxx(skbuffer, 23652, 8100, 11745, r_ab, 4, 4);

            erirec::comp_bra_hrr_electron_repulsion_dhxx(skbuffer, 30942, 11745, 16848, r_ab, 4, 4);

            erirec::comp_bra_hrr_electron_repulsion_fgxx(skbuffer, 41148, 23652, 30942, r_ab, 4, 4);

            t4cfunc::bra_transform<3, 4>(sbuffer, 0, skbuffer, 41148, 4, 4);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 4, 4, 4, j, ket_range, bra_eq_ket);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionRecFGGG_hpp */
