#ifndef ElectronRepulsionGeom1100RecFFFF_hpp
#define ElectronRepulsionGeom1100RecFFFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecXXDF.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXFF.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
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
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSI.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSI.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSI.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSI.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
#include "ElectronRepulsionPrimRecSISI.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSF.hpp"
#include "ElectronRepulsionPrimRecSKSG.hpp"
#include "ElectronRepulsionPrimRecSKSH.hpp"
#include "ElectronRepulsionPrimRecSKSI.hpp"
#include "ElectronRepulsionPrimRecSLSF.hpp"
#include "ElectronRepulsionPrimRecSLSG.hpp"
#include "ElectronRepulsionPrimRecSLSH.hpp"
#include "ElectronRepulsionPrimRecSLSI.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
#include "ElectronRepulsionPrimRecSPSI.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSH.hpp"
#include "ElectronRepulsionPrimRecSSSI.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FF|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_ffff(T& distributor,
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

    CSimdArray<double> pbuffer(39507, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(23680, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(170960, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(290276, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(21609, 1);

    // setup Boys fuction data

    const CBoysFunc<14> bf_table;

    CSimdArray<double> bf_data(16, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 15, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 15, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 15, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 15);
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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 882, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 885, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 888, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 891, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 894, 7, 8, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 897, 3, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 906, 4, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 915, 5, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 924, 6, 30, 33, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 933, 7, 33, 36, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 942, 8, 36, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 951, 21, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 969, 24, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 987, 27, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1005, 30, 81, 87, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1023, 33, 87, 93, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1041, 36, 93, 99, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1059, 39, 99, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1077, 63, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1107, 69, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1137, 75, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1167, 81, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1197, 87, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1227, 93, 185, 195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1257, 99, 195, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1287, 105, 205, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1317, 145, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1362, 155, 270, 285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1407, 165, 285, 300, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1452, 175, 300, 315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1497, 185, 315, 330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1542, 195, 330, 345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1587, 205, 345, 360, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1632, 215, 360, 375, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1677, 270, 420, 441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1740, 285, 441, 462, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1803, 300, 462, 483, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1866, 315, 483, 504, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1929, 330, 504, 525, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1992, 345, 525, 546, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2055, 360, 546, 567, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2118, 375, 567, 588, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2181, 441, 630, 658, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2265, 462, 658, 686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2349, 483, 686, 714, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2433, 504, 714, 742, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2517, 525, 742, 770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2601, 546, 770, 798, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2685, 567, 798, 826, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2769, 588, 826, 854, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2853, 3, 4, 882, 885, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2859, 4, 5, 885, 888, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2865, 5, 6, 888, 891, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2871, 6, 7, 891, 894, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2877, 21, 24, 882, 897, 906, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2895, 24, 27, 885, 906, 915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2913, 27, 30, 888, 915, 924, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2931, 30, 33, 891, 924, 933, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2949, 33, 36, 894, 933, 942, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2967, 63, 69, 897, 951, 969, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3003, 69, 75, 906, 969, 987, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3039, 75, 81, 915, 987, 1005, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3075, 81, 87, 924, 1005, 1023, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3111, 87, 93, 933, 1023, 1041, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3147, 93, 99, 942, 1041, 1059, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3183, 135, 145, 951, 1077, 1107, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3243, 145, 155, 969, 1107, 1137, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3303, 155, 165, 987, 1137, 1167, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3363, 165, 175, 1005, 1167, 1197, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3423, 175, 185, 1023, 1197, 1227, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3483, 185, 195, 1041, 1227, 1257, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3543, 195, 205, 1059, 1257, 1287, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3603, 255, 270, 1107, 1317, 1362, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3693, 270, 285, 1137, 1362, 1407, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3783, 285, 300, 1167, 1407, 1452, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3873, 300, 315, 1197, 1452, 1497, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3963, 315, 330, 1227, 1497, 1542, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4053, 330, 345, 1257, 1542, 1587, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4143, 345, 360, 1287, 1587, 1632, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4233, 420, 441, 1362, 1677, 1740, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4359, 441, 462, 1407, 1740, 1803, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4485, 462, 483, 1452, 1803, 1866, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4611, 483, 504, 1497, 1866, 1929, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4737, 504, 525, 1542, 1929, 1992, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4863, 525, 546, 1587, 1992, 2055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4989, 546, 567, 1632, 2055, 2118, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5115, 630, 658, 1740, 2181, 2265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5283, 658, 686, 1803, 2265, 2349, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5451, 686, 714, 1866, 2349, 2433, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5619, 714, 742, 1929, 2433, 2517, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5787, 742, 770, 1992, 2517, 2601, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5955, 770, 798, 2055, 2601, 2685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6123, 798, 826, 2118, 2685, 2769, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 6291, 882, 885, 2853, 2859, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 6301, 885, 888, 2859, 2865, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 6311, 888, 891, 2865, 2871, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 6321, 897, 906, 2853, 2877, 2895, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 6351, 906, 915, 2859, 2895, 2913, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 6381, 915, 924, 2865, 2913, 2931, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 6411, 924, 933, 2871, 2931, 2949, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 6441, 951, 969, 2877, 2967, 3003, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 6501, 969, 987, 2895, 3003, 3039, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 6561, 987, 1005, 2913, 3039, 3075, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 6621, 1005, 1023, 2931, 3075, 3111, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 6681, 1023, 1041, 2949, 3111, 3147, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 6741, 1077, 1107, 2967, 3183, 3243, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 6841, 1107, 1137, 3003, 3243, 3303, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 6941, 1137, 1167, 3039, 3303, 3363, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 7041, 1167, 1197, 3075, 3363, 3423, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 7141, 1197, 1227, 3111, 3423, 3483, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 7241, 1227, 1257, 3147, 3483, 3543, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 7341, 1317, 1362, 3243, 3603, 3693, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 7491, 1362, 1407, 3303, 3693, 3783, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 7641, 1407, 1452, 3363, 3783, 3873, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 7791, 1452, 1497, 3423, 3873, 3963, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 7941, 1497, 1542, 3483, 3963, 4053, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8091, 1542, 1587, 3543, 4053, 4143, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 8241, 1677, 1740, 3693, 4233, 4359, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 8451, 1740, 1803, 3783, 4359, 4485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 8661, 1803, 1866, 3873, 4485, 4611, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 8871, 1866, 1929, 3963, 4611, 4737, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9081, 1929, 1992, 4053, 4737, 4863, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9291, 1992, 2055, 4143, 4863, 4989, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 9501, 2181, 2265, 4359, 5115, 5283, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 9781, 2265, 2349, 4485, 5283, 5451, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10061, 2349, 2433, 4611, 5451, 5619, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10341, 2433, 2517, 4737, 5619, 5787, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10621, 2517, 2601, 4863, 5787, 5955, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10901, 2601, 2685, 4989, 5955, 6123, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 11181, 2853, 2859, 6291, 6301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 11196, 2859, 2865, 6301, 6311, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 11211, 2877, 2895, 6291, 6321, 6351, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 11256, 2895, 2913, 6301, 6351, 6381, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 11301, 2913, 2931, 6311, 6381, 6411, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 11346, 2967, 3003, 6321, 6441, 6501, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 11436, 3003, 3039, 6351, 6501, 6561, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 11526, 3039, 3075, 6381, 6561, 6621, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 11616, 3075, 3111, 6411, 6621, 6681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 11706, 3183, 3243, 6441, 6741, 6841, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 11856, 3243, 3303, 6501, 6841, 6941, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 12006, 3303, 3363, 6561, 6941, 7041, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 12156, 3363, 3423, 6621, 7041, 7141, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 12306, 3423, 3483, 6681, 7141, 7241, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 12456, 3603, 3693, 6841, 7341, 7491, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 12681, 3693, 3783, 6941, 7491, 7641, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 12906, 3783, 3873, 7041, 7641, 7791, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 13131, 3873, 3963, 7141, 7791, 7941, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 13356, 3963, 4053, 7241, 7941, 8091, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 13581, 4233, 4359, 7491, 8241, 8451, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 13896, 4359, 4485, 7641, 8451, 8661, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 14211, 4485, 4611, 7791, 8661, 8871, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 14526, 4611, 4737, 7941, 8871, 9081, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 14841, 4737, 4863, 8091, 9081, 9291, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 15156, 5115, 5283, 8451, 9501, 9781, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 15576, 5283, 5451, 8661, 9781, 10061, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 15996, 5451, 5619, 8871, 10061, 10341, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 16416, 5619, 5787, 9081, 10341, 10621, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 16836, 5787, 5955, 9291, 10621, 10901, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 17256, 6291, 6301, 11181, 11196, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 17277, 6321, 6351, 11181, 11211, 11256, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 17340, 6351, 6381, 11196, 11256, 11301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 17403, 6441, 6501, 11211, 11346, 11436, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 17529, 6501, 6561, 11256, 11436, 11526, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 17655, 6561, 6621, 11301, 11526, 11616, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 17781, 6741, 6841, 11346, 11706, 11856, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 17991, 6841, 6941, 11436, 11856, 12006, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 18201, 6941, 7041, 11526, 12006, 12156, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 18411, 7041, 7141, 11616, 12156, 12306, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 18621, 7341, 7491, 11856, 12456, 12681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 18936, 7491, 7641, 12006, 12681, 12906, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 19251, 7641, 7791, 12156, 12906, 13131, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 19566, 7791, 7941, 12306, 13131, 13356, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 19881, 8241, 8451, 12681, 13581, 13896, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 20322, 8451, 8661, 12906, 13896, 14211, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 20763, 8661, 8871, 13131, 14211, 14526, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 21204, 8871, 9081, 13356, 14526, 14841, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 21645, 9501, 9781, 13896, 15156, 15576, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 22233, 9781, 10061, 14211, 15576, 15996, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 22821, 10061, 10341, 14526, 15996, 16416, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 23409, 10341, 10621, 14841, 16416, 16836, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 23997, 11211, 11256, 17256, 17277, 17340, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 24081, 11346, 11436, 17277, 17403, 17529, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 24249, 11436, 11526, 17340, 17529, 17655, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 24417, 11706, 11856, 17403, 17781, 17991, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 24697, 11856, 12006, 17529, 17991, 18201, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 24977, 12006, 12156, 17655, 18201, 18411, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 25257, 12456, 12681, 17991, 18621, 18936, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 25677, 12681, 12906, 18201, 18936, 19251, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 26097, 12906, 13131, 18411, 19251, 19566, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 26517, 13581, 13896, 18936, 19881, 20322, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 27105, 13896, 14211, 19251, 20322, 20763, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 27693, 14211, 14526, 19566, 20763, 21204, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 28281, 15156, 15576, 20322, 21645, 22233, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 29065, 15576, 15996, 20763, 22233, 22821, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 29849, 15996, 16416, 21204, 22821, 23409, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 30633, 17403, 17529, 23997, 24081, 24249, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 30849, 17781, 17991, 24081, 24417, 24697, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 31209, 17991, 18201, 24249, 24697, 24977, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 31569, 18621, 18936, 24697, 25257, 25677, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 32109, 18936, 19251, 24977, 25677, 26097, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksh(pbuffer, 32649, 19881, 20322, 25677, 26517, 27105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksh(pbuffer, 33405, 20322, 20763, 26097, 27105, 27693, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksi(pbuffer, 34161, 21645, 22233, 27105, 28281, 29065, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksi(pbuffer, 35169, 22233, 22821, 27693, 29065, 29849, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsf(pbuffer, 36177, 24417, 24697, 30633, 30849, 31209, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsg(pbuffer, 36627, 25257, 25677, 31209, 31569, 32109, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsh(pbuffer, 37302, 26517, 27105, 32109, 32649, 33405, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsi(pbuffer, 38247, 28281, 29065, 33405, 34161, 35169, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 3183, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 3603, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 150, pbuffer, 4233, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 276, pbuffer, 5115, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 444, pbuffer, 6741, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 544, pbuffer, 7341, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 694, pbuffer, 8241, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 904, pbuffer, 9501, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1184, pbuffer, 11706, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1334, pbuffer, 12456, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1559, pbuffer, 13581, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1874, pbuffer, 15156, 420, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {11706, 11856});

                pbuffer.scale(2.0 * b_exp, {12456, 12681});

                pbuffer.scale(2.0 * b_exp, {13581, 13896});

                pbuffer.scale(2.0 * b_exp, {15156, 15576});

                pbuffer.scale(2.0 * b_exp, {17781, 17991});

                pbuffer.scale(2.0 * b_exp, {18621, 18936});

                pbuffer.scale(2.0 * b_exp, {19881, 20322});

                pbuffer.scale(2.0 * b_exp, {21645, 22233});

                pbuffer.scale(2.0 * b_exp, {24417, 24697});

                pbuffer.scale(2.0 * b_exp, {25257, 25677});

                pbuffer.scale(2.0 * b_exp, {26517, 27105});

                pbuffer.scale(2.0 * b_exp, {28281, 29065});

                t2cfunc::reduce(cbuffer, 2294, pbuffer, 11706, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2444, pbuffer, 12456, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2669, pbuffer, 13581, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2984, pbuffer, 15156, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3404, pbuffer, 17781, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3614, pbuffer, 18621, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3929, pbuffer, 19881, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4370, pbuffer, 21645, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4958, pbuffer, 24417, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5238, pbuffer, 25257, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5658, pbuffer, 26517, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6246, pbuffer, 28281, 784, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {3183, 3243});

                pbuffer.scale(2.0 * a_exp, {3603, 3693});

                pbuffer.scale(2.0 * a_exp, {4233, 4359});

                pbuffer.scale(2.0 * a_exp, {5115, 5283});

                pbuffer.scale(2.0 * a_exp, {6741, 6841});

                pbuffer.scale(2.0 * a_exp, {7341, 7491});

                pbuffer.scale(2.0 * a_exp, {8241, 8451});

                pbuffer.scale(2.0 * a_exp, {9501, 9781});

                pbuffer.scale(a_exp / b_exp, {11706, 11856});

                pbuffer.scale(a_exp / b_exp, {12456, 12681});

                pbuffer.scale(a_exp / b_exp, {13581, 13896});

                pbuffer.scale(a_exp / b_exp, {15156, 15576});

                pbuffer.scale(a_exp / b_exp, {17781, 17991});

                pbuffer.scale(a_exp / b_exp, {18621, 18936});

                pbuffer.scale(a_exp / b_exp, {19881, 20322});

                pbuffer.scale(a_exp / b_exp, {21645, 22233});

                pbuffer.scale(a_exp / b_exp, {24417, 24697});

                pbuffer.scale(a_exp / b_exp, {25257, 25677});

                pbuffer.scale(a_exp / b_exp, {26517, 27105});

                pbuffer.scale(a_exp / b_exp, {28281, 29065});

                t2cfunc::reduce(cbuffer, 7030, pbuffer, 3183, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7090, pbuffer, 3603, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7180, pbuffer, 4233, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7306, pbuffer, 5115, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7474, pbuffer, 6741, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7574, pbuffer, 7341, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7724, pbuffer, 8241, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7934, pbuffer, 9501, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8214, pbuffer, 11706, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8364, pbuffer, 12456, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8589, pbuffer, 13581, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8904, pbuffer, 15156, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9324, pbuffer, 17781, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9534, pbuffer, 18621, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9849, pbuffer, 19881, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10290, pbuffer, 21645, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10878, pbuffer, 24417, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11158, pbuffer, 25257, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11578, pbuffer, 26517, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12166, pbuffer, 28281, 784, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {11706, 11856});

                pbuffer.scale(2.0 * b_exp, {12456, 12681});

                pbuffer.scale(2.0 * b_exp, {13581, 13896});

                pbuffer.scale(2.0 * b_exp, {15156, 15576});

                pbuffer.scale(2.0 * b_exp, {17781, 17991});

                pbuffer.scale(2.0 * b_exp, {18621, 18936});

                pbuffer.scale(2.0 * b_exp, {19881, 20322});

                pbuffer.scale(2.0 * b_exp, {21645, 22233});

                pbuffer.scale(2.0 * b_exp, {24417, 24697});

                pbuffer.scale(2.0 * b_exp, {25257, 25677});

                pbuffer.scale(2.0 * b_exp, {26517, 27105});

                pbuffer.scale(2.0 * b_exp, {28281, 29065});

                pbuffer.scale(4.0 * a_exp * b_exp, {30849, 31209});

                pbuffer.scale(4.0 * a_exp * b_exp, {31569, 32109});

                pbuffer.scale(4.0 * a_exp * b_exp, {32649, 33405});

                pbuffer.scale(4.0 * a_exp * b_exp, {34161, 35169});

                pbuffer.scale(4.0 * a_exp * b_exp, {36177, 36627});

                pbuffer.scale(4.0 * a_exp * b_exp, {36627, 37302});

                pbuffer.scale(4.0 * a_exp * b_exp, {37302, 38247});

                pbuffer.scale(4.0 * a_exp * b_exp, {38247, 39507});

                t2cfunc::reduce(cbuffer, 12950, pbuffer, 11706, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13100, pbuffer, 12456, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13325, pbuffer, 13581, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13640, pbuffer, 15156, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14060, pbuffer, 17781, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14270, pbuffer, 18621, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14585, pbuffer, 19881, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15026, pbuffer, 21645, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15614, pbuffer, 24417, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15894, pbuffer, 25257, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16314, pbuffer, 26517, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16902, pbuffer, 28281, 784, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 17686, pbuffer, 30849, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18046, pbuffer, 31569, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18586, pbuffer, 32649, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 19342, pbuffer, 34161, 1008, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 20350, pbuffer, 36177, 450, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 20800, pbuffer, 36627, 675, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 21475, pbuffer, 37302, 945, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 22420, pbuffer, 38247, 1260, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 60, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 180, cbuffer, 60, 150, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 450, cbuffer, 150, 276, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 828, 0, 180, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 1188, 180, 450, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 1728, 828, 1188, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2328, cbuffer, 444, 544, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 2628, cbuffer, 544, 694, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 3078, cbuffer, 694, 904, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 3708, 2328, 2628, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 4308, 2628, 3078, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 5208, 3708, 4308, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9208, cbuffer, 1184, 1334, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 9658, cbuffer, 1334, 1559, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 10333, cbuffer, 1559, 1874, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 11278, 9208, 9658, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 12178, 9658, 10333, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 13528, 11278, 12178, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 25828, cbuffer, 2294, 2444, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 26278, cbuffer, 2444, 2669, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 26953, cbuffer, 2669, 2984, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 27898, 25828, 26278, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 28798, 26278, 26953, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 30148, 27898, 28798, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 31648, cbuffer, 3404, 3614, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 32278, cbuffer, 3614, 3929, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 33223, cbuffer, 3929, 4370, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 34546, 31648, 32278, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 35806, 32278, 33223, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 37696, 34546, 35806, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 39796, cbuffer, 4958, 5238, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 40636, cbuffer, 5238, 5658, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 41896, cbuffer, 5658, 6246, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 43660, 39796, 40636, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 45340, 40636, 41896, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 47860, 43660, 45340, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 50660, cbuffer, 7030, 7090, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 50840, cbuffer, 7090, 7180, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 51110, cbuffer, 7180, 7306, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 51488, 50660, 50840, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 51848, 50840, 51110, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 52388, 51488, 51848, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 52988, cbuffer, 7474, 7574, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 53288, cbuffer, 7574, 7724, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 53738, cbuffer, 7724, 7934, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 54368, 52988, 53288, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 54968, 53288, 53738, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 55868, 54368, 54968, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 59868, cbuffer, 8214, 8364, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 60318, cbuffer, 8364, 8589, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 60993, cbuffer, 8589, 8904, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 61938, 59868, 60318, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 62838, 60318, 60993, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 64188, 61938, 62838, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 70188, cbuffer, 9324, 9534, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 70818, cbuffer, 9534, 9849, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 71763, cbuffer, 9849, 10290, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 73086, 70188, 70818, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 74346, 70818, 71763, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 76236, 73086, 74346, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 84636, cbuffer, 10878, 11158, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 85476, cbuffer, 11158, 11578, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 86736, cbuffer, 11578, 12166, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 88500, 84636, 85476, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 90180, 85476, 86736, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 92700, 88500, 90180, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 114700, cbuffer, 12950, 13100, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 115150, cbuffer, 13100, 13325, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 115825, cbuffer, 13325, 13640, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 116770, 114700, 115150, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 117670, 115150, 115825, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 119020, 116770, 117670, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 120520, cbuffer, 14060, 14270, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 121150, cbuffer, 14270, 14585, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 122095, cbuffer, 14585, 15026, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 123418, 120520, 121150, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 124678, 121150, 122095, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 126568, 123418, 124678, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 128668, cbuffer, 15614, 15894, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 129508, cbuffer, 15894, 16314, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 130768, cbuffer, 16314, 16902, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 132532, 128668, 129508, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 134212, 129508, 130768, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 136732, 132532, 134212, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 139532, cbuffer, 17686, 18046, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 140612, cbuffer, 18046, 18586, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 142232, cbuffer, 18586, 19342, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 144500, 139532, 140612, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 146660, 140612, 142232, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 149900, 144500, 146660, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 153500, cbuffer, 20350, 20800, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 154850, cbuffer, 20800, 21475, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 156875, cbuffer, 21475, 22420, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 159710, 153500, 154850, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 162410, 154850, 156875, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 166460, 159710, 162410, cfactors, 6, 0, 8);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 1728, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 294, ckbuffer, 5208, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 6664, ckbuffer, 13528, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 253183, ckbuffer, 30148, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 253918, ckbuffer, 37696, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 254947, ckbuffer, 47860, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 256319, ckbuffer, 52388, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 256613, ckbuffer, 55868, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 258573, ckbuffer, 64188, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 261513, ckbuffer, 76236, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 265629, ckbuffer, 92700, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 283171, ckbuffer, 119020, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 283906, ckbuffer, 126568, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 284935, ckbuffer, 136732, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 286307, ckbuffer, 149900, 0, 7);

            t4cfunc::ket_transform<3, 3>(skbuffer, 288071, ckbuffer, 166460, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 40915, 294, 6664, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 276409, 256613, 258573, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 277879, 258573, 261513, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 280084, 261513, 265629, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 46795, 294, 276409, 277879, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 71050, 6664, 277879, 280084, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 134113, 40915, 46795, 71050, r_ab, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 784, 0, 253183, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 7399, 294, 253918, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 16219, 6664, 254947, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 42385, 294, 784, 7399, r_ab, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pgxx(skbuffer, 64435, 6664, 7399, 16219, r_ab, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dfxx(skbuffer, 125293, 40915, 42385, 64435, r_ab, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 257103, 256319, 283171, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 259308, 256613, 283906, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 262542, 258573, 284935, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 267001, 261513, 286307, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_skxx(skbuffer, 271117, 265629, 288071, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 2254, 256613, 257103, 259308, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 9604, 258573, 259308, 262542, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 19306, 261513, 262542, 267001, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sixx(skbuffer, 28567, 265629, 267001, 271117, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 51205, 784, 276409, 2254, 9604, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 77665, 7399, 277879, 9604, 19306, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_phxx(skbuffer, 97510, 16219, 280084, 19306, 28567, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 142933, 42385, 46795, 51205, 77665, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dgxx(skbuffer, 169393, 64435, 71050, 77665, 97510, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ffxx(skbuffer, 209083, 125293, 134113, 142933, 169393, r_ab, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 209083, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2401, skbuffer, 213983, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 4802, skbuffer, 218883, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 7203, skbuffer, 223783, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 9604, skbuffer, 228683, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 12005, skbuffer, 233583, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 14406, skbuffer, 238483, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 16807, skbuffer, 243383, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 19208, skbuffer, 248283, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFFFF_hpp */
