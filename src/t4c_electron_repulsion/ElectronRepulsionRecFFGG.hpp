#ifndef ElectronRepulsionRecFFGG_hpp
#define ElectronRepulsionRecFFGG_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "ElectronRepulsionContrRecDFXX.hpp"
#include "ElectronRepulsionContrRecDGXX.hpp"
#include "ElectronRepulsionContrRecFFXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
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
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSI.hpp"
#include "ElectronRepulsionPrimRecSGSK.hpp"
#include "ElectronRepulsionPrimRecSGSL.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSI.hpp"
#include "ElectronRepulsionPrimRecSHSK.hpp"
#include "ElectronRepulsionPrimRecSHSL.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
#include "ElectronRepulsionPrimRecSISI.hpp"
#include "ElectronRepulsionPrimRecSISK.hpp"
#include "ElectronRepulsionPrimRecSISL.hpp"
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
#include "GtoPairBlock.hpp"
#include "SimdArray.hpp"
#include "T2CUtils.hpp"
#include "T4CUtils.hpp"

namespace erirec {  // erirec namespace

/// @brief Computes (FF|1/|r-r'||GG)  integrals for two basis function pairs blocks.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
/// @param bra_eq_ket True if basis function pairs blocks on bra and ket are the same, False otherwise.
template <class T>
inline auto
comp_electron_repulsion_ffgg(T&                               distributor,
                             const CGtoPairBlock&             bra_gto_pair_block,
                             const CGtoPairBlock&             ket_gto_pair_block,
                             const std::pair<size_t, size_t>& bra_indices,
                             const std::pair<size_t, size_t>& ket_indices,
                             const bool                       bra_eq_ket) -> void
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

    CSimdArray<double> pbuffer(32615, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(10730, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(93906, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(37422, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3969, 1);

    // setup Boys fuction data

    const CBoysFunc<14> bf_table;

    CSimdArray<double> bf_data(16, ket_npgtos);

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

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1170, 630, 658, 882, 918, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1215, 658, 686, 918, 954, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1260, 686, 714, 954, 990, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1305, 714, 742, 990, 1026, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1350, 742, 770, 1026, 1062, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1395, 770, 798, 1062, 1098, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssl(pbuffer, 1440, 798, 826, 1098, 1134, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1485, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1488, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1491, 4, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1500, 5, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1509, 6, 30, 33, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1518, 24, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1536, 27, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1554, 30, 81, 87, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1572, 33, 87, 93, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1590, 69, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1620, 75, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1650, 81, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1680, 87, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1710, 93, 185, 195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1740, 145, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1785, 155, 270, 285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1830, 165, 285, 300, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1875, 175, 300, 315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1920, 185, 315, 330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1965, 195, 330, 345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2010, 270, 420, 441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2073, 285, 441, 462, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2136, 300, 462, 483, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2199, 315, 483, 504, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2262, 330, 504, 525, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2325, 345, 525, 546, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2388, 441, 630, 658, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2472, 462, 658, 686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2556, 483, 686, 714, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2640, 504, 714, 742, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2724, 525, 742, 770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2808, 546, 770, 798, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2892, 658, 882, 918, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3000, 686, 918, 954, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3108, 714, 954, 990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3216, 742, 990, 1026, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3324, 770, 1026, 1062, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3432, 798, 1062, 1098, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 3540, 918, 1170, 1215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 3675, 954, 1215, 1260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 3810, 990, 1260, 1305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 3945, 1026, 1305, 1350, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 4080, 1062, 1350, 1395, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsl(pbuffer, 4215, 1098, 1395, 1440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 4350, 4, 5, 1485, 1488, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 4356, 24, 27, 1485, 1491, 1500, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 4374, 27, 30, 1488, 1500, 1509, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 4392, 69, 75, 1491, 1518, 1536, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 4428, 75, 81, 1500, 1536, 1554, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 4464, 81, 87, 1509, 1554, 1572, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 4500, 145, 155, 1518, 1590, 1620, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 4560, 155, 165, 1536, 1620, 1650, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 4620, 165, 175, 1554, 1650, 1680, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 4680, 175, 185, 1572, 1680, 1710, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4740, 255, 270, 1590, 1740, 1785, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4830, 270, 285, 1620, 1785, 1830, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4920, 285, 300, 1650, 1830, 1875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 5010, 300, 315, 1680, 1875, 1920, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 5100, 315, 330, 1710, 1920, 1965, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5190, 420, 441, 1785, 2010, 2073, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5316, 441, 462, 1830, 2073, 2136, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5442, 462, 483, 1875, 2136, 2199, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5568, 483, 504, 1920, 2199, 2262, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5694, 504, 525, 1965, 2262, 2325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5820, 630, 658, 2073, 2388, 2472, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5988, 658, 686, 2136, 2472, 2556, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6156, 686, 714, 2199, 2556, 2640, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6324, 714, 742, 2262, 2640, 2724, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6492, 742, 770, 2325, 2724, 2808, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6660, 882, 918, 2472, 2892, 3000, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6876, 918, 954, 2556, 3000, 3108, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7092, 954, 990, 2640, 3108, 3216, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7308, 990, 1026, 2724, 3216, 3324, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7524, 1026, 1062, 2808, 3324, 3432, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsl(pbuffer, 7740, 1170, 1215, 3000, 3540, 3675, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsl(pbuffer, 8010, 1215, 1260, 3108, 3675, 3810, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsl(pbuffer, 8280, 1260, 1305, 3216, 3810, 3945, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsl(pbuffer, 8550, 1305, 1350, 3324, 3945, 4080, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsl(pbuffer, 8820, 1350, 1395, 3432, 4080, 4215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 9090, 1491, 1500, 4350, 4356, 4374, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 9120, 1518, 1536, 4356, 4392, 4428, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 9180, 1536, 1554, 4374, 4428, 4464, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 9240, 1590, 1620, 4392, 4500, 4560, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 9340, 1620, 1650, 4428, 4560, 4620, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 9440, 1650, 1680, 4464, 4620, 4680, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 9540, 1740, 1785, 4500, 4740, 4830, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 9690, 1785, 1830, 4560, 4830, 4920, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 9840, 1830, 1875, 4620, 4920, 5010, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 9990, 1875, 1920, 4680, 5010, 5100, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 10140, 2010, 2073, 4830, 5190, 5316, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 10350, 2073, 2136, 4920, 5316, 5442, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 10560, 2136, 2199, 5010, 5442, 5568, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 10770, 2199, 2262, 5100, 5568, 5694, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10980, 2388, 2472, 5316, 5820, 5988, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 11260, 2472, 2556, 5442, 5988, 6156, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 11540, 2556, 2640, 5568, 6156, 6324, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 11820, 2640, 2724, 5694, 6324, 6492, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 12100, 2892, 3000, 5988, 6660, 6876, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 12460, 3000, 3108, 6156, 6876, 7092, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 12820, 3108, 3216, 6324, 7092, 7308, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 13180, 3216, 3324, 6492, 7308, 7524, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsl(pbuffer, 13540, 3540, 3675, 6876, 7740, 8010, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsl(pbuffer, 13990, 3675, 3810, 7092, 8010, 8280, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsl(pbuffer, 14440, 3810, 3945, 7308, 8280, 8550, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsl(pbuffer, 14890, 3945, 4080, 7524, 8550, 8820, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 15340, 4392, 4428, 9090, 9120, 9180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 15430, 4500, 4560, 9120, 9240, 9340, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 15580, 4560, 4620, 9180, 9340, 9440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 15730, 4740, 4830, 9240, 9540, 9690, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 15955, 4830, 4920, 9340, 9690, 9840, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 16180, 4920, 5010, 9440, 9840, 9990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 16405, 5190, 5316, 9690, 10140, 10350, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 16720, 5316, 5442, 9840, 10350, 10560, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 17035, 5442, 5568, 9990, 10560, 10770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 17350, 5820, 5988, 10350, 10980, 11260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 17770, 5988, 6156, 10560, 11260, 11540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 18190, 6156, 6324, 10770, 11540, 11820, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 18610, 6660, 6876, 11260, 12100, 12460, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 19150, 6876, 7092, 11540, 12460, 12820, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 19690, 7092, 7308, 11820, 12820, 13180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsl(pbuffer, 20230, 7740, 8010, 12460, 13540, 13990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsl(pbuffer, 20905, 8010, 8280, 12820, 13990, 14440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsl(pbuffer, 21580, 8280, 8550, 13180, 14440, 14890, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 22255, 9240, 9340, 15340, 15430, 15580, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 22465, 9540, 9690, 15430, 15730, 15955, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 22780, 9690, 9840, 15580, 15955, 16180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 23095, 10140, 10350, 15955, 16405, 16720, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 23536, 10350, 10560, 16180, 16720, 17035, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 23977, 10980, 11260, 16720, 17350, 17770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 24565, 11260, 11540, 17035, 17770, 18190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsk(pbuffer, 25153, 12100, 12460, 17770, 18610, 19150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsk(pbuffer, 25909, 12460, 12820, 18190, 19150, 19690, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsl(pbuffer, 26665, 13540, 13990, 19150, 20230, 20905, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsl(pbuffer, 27610, 13990, 14440, 19690, 20905, 21580, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 28555, 15730, 15955, 22255, 22465, 22780, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 28975, 16405, 16720, 22780, 23095, 23536, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 29563, 17350, 17770, 23536, 23977, 24565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisk(pbuffer, 30347, 18610, 19150, 24565, 25153, 25909, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisl(pbuffer, 31355, 20230, 20905, 25909, 26665, 27610, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 9540, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 150, pbuffer, 10140, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 360, pbuffer, 10980, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 640, pbuffer, 12100, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1000, pbuffer, 13540, 450, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1450, pbuffer, 15730, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1675, pbuffer, 16405, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1990, pbuffer, 17350, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2410, pbuffer, 18610, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2950, pbuffer, 20230, 675, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3625, pbuffer, 22465, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3940, pbuffer, 23095, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4381, pbuffer, 23977, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4969, pbuffer, 25153, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5725, pbuffer, 26665, 945, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6670, pbuffer, 28555, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7090, pbuffer, 28975, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7678, pbuffer, 29563, 784, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8462, pbuffer, 30347, 1008, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9470, pbuffer, 31355, 1260, ket_width, ket_npgtos);
            }

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 0, cbuffer, 0, 150, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 450, cbuffer, 150, 360, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 1080, cbuffer, 360, 640, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpk(ckbuffer, 1920, cbuffer, 640, 1000, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 3000, 0, 450, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 3900, 450, 1080, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdi(ckbuffer, 5160, 1080, 1920, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 6840, 3000, 3900, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxfh(ckbuffer, 8340, 3900, 5160, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxgg(ckbuffer, 10440, 6840, 8340, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 12690, cbuffer, 1450, 1675, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 13365, cbuffer, 1675, 1990, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 14310, cbuffer, 1990, 2410, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpk(ckbuffer, 15570, cbuffer, 2410, 2950, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 17190, 12690, 13365, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 18540, 13365, 14310, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdi(ckbuffer, 20430, 14310, 15570, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 22950, 17190, 18540, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxfh(ckbuffer, 25200, 18540, 20430, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxgg(ckbuffer, 28350, 22950, 25200, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 31725, cbuffer, 3625, 3940, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 32670, cbuffer, 3940, 4381, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 33993, cbuffer, 4381, 4969, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpk(ckbuffer, 35757, cbuffer, 4969, 5725, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 38025, 31725, 32670, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 39915, 32670, 33993, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdi(ckbuffer, 42561, 33993, 35757, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 46089, 38025, 39915, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxfh(ckbuffer, 49239, 39915, 42561, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxgg(ckbuffer, 53649, 46089, 49239, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 58374, cbuffer, 6670, 7090, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 59634, cbuffer, 7090, 7678, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpi(ckbuffer, 61398, cbuffer, 7678, 8462, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpk(ckbuffer, 63750, cbuffer, 8462, 9470, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 66774, 58374, 59634, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdh(ckbuffer, 69294, 59634, 61398, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdi(ckbuffer, 72822, 61398, 63750, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxfg(ckbuffer, 77526, 66774, 69294, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxfh(ckbuffer, 81726, 69294, 72822, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxgg(ckbuffer, 87606, 77526, 81726, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<4, 4>(skbuffer, 0, ckbuffer, 10440, 0, 3);

            t4cfunc::ket_transform<4, 4>(skbuffer, 810, ckbuffer, 28350, 0, 4);

            t4cfunc::ket_transform<4, 4>(skbuffer, 2025, ckbuffer, 53649, 0, 5);

            t4cfunc::ket_transform<4, 4>(skbuffer, 3726, ckbuffer, 87606, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 5994, 0, 810, r_ab, 4, 4);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 8424, 810, 2025, r_ab, 4, 4);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 12069, 2025, 3726, r_ab, 4, 4);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 17172, 5994, 8424, r_ab, 4, 4);

            erirec::comp_bra_hrr_electron_repulsion_dgxx(skbuffer, 22032, 8424, 12069, r_ab, 4, 4);

            erirec::comp_bra_hrr_electron_repulsion_ffxx(skbuffer, 29322, 17172, 22032, r_ab, 4, 4);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 29322, 4, 4);

            const bool diagonal = bra_eq_ket && (j >= ket_range.first) && (j < ket_range.second);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 4, 4, j, ket_range, diagonal);
        }
    }
}

}  // namespace erirec

#endif /* ElectronRepulsionRecFFGG_hpp */
