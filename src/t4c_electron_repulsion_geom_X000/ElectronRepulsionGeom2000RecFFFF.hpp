#ifndef ElectronRepulsionGeom2000RecFFFF_hpp
#define ElectronRepulsionGeom2000RecFFFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDFXX.hpp"
#include "ElectronRepulsionContrRecDGXX.hpp"
#include "ElectronRepulsionContrRecDHXX.hpp"
#include "ElectronRepulsionContrRecDIXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecPIXX.hpp"
#include "ElectronRepulsionContrRecPKXX.hpp"
#include "ElectronRepulsionContrRecXXDF.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXFF.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
#include "ElectronRepulsionGeom1000ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDGXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecFFXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPHXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"
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

/// @brief Computes d^(2)/dA^(2)(FF|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_ffff(T& distributor,
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

    CSimdArray<double> cbuffer(18796, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(98552, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(214277, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(14406, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 6741, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 7341, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 250, pbuffer, 8241, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 460, pbuffer, 9501, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 740, pbuffer, 11706, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 890, pbuffer, 12456, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1115, pbuffer, 13581, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1430, pbuffer, 15156, 420, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {6741, 6841});

                pbuffer.scale(2.0 * a_exp, {7341, 7491});

                pbuffer.scale(2.0 * a_exp, {8241, 8451});

                pbuffer.scale(2.0 * a_exp, {9501, 9781});

                pbuffer.scale(2.0 * a_exp, {11706, 11856});

                pbuffer.scale(2.0 * a_exp, {12456, 12681});

                pbuffer.scale(2.0 * a_exp, {13581, 13896});

                pbuffer.scale(2.0 * a_exp, {15156, 15576});

                pbuffer.scale(2.0 * a_exp, {17781, 17991});

                pbuffer.scale(2.0 * a_exp, {18621, 18936});

                pbuffer.scale(2.0 * a_exp, {19881, 20322});

                pbuffer.scale(2.0 * a_exp, {21645, 22233});

                pbuffer.scale(2.0 * a_exp, {24417, 24697});

                pbuffer.scale(2.0 * a_exp, {25257, 25677});

                pbuffer.scale(2.0 * a_exp, {26517, 27105});

                pbuffer.scale(2.0 * a_exp, {28281, 29065});

                t2cfunc::reduce(cbuffer, 1850, pbuffer, 6741, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1950, pbuffer, 7341, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2100, pbuffer, 8241, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2310, pbuffer, 9501, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2590, pbuffer, 11706, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2740, pbuffer, 12456, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2965, pbuffer, 13581, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3280, pbuffer, 15156, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3700, pbuffer, 17781, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3910, pbuffer, 18621, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4225, pbuffer, 19881, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4666, pbuffer, 21645, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5254, pbuffer, 24417, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5534, pbuffer, 25257, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5954, pbuffer, 26517, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6542, pbuffer, 28281, 784, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {6741, 6841});

                pbuffer.scale(2.0 * a_exp, {7341, 7491});

                pbuffer.scale(2.0 * a_exp, {8241, 8451});

                pbuffer.scale(2.0 * a_exp, {9501, 9781});

                pbuffer.scale(2.0 * a_exp, {11706, 11856});

                pbuffer.scale(2.0 * a_exp, {12456, 12681});

                pbuffer.scale(2.0 * a_exp, {13581, 13896});

                pbuffer.scale(2.0 * a_exp, {15156, 15576});

                pbuffer.scale(2.0 * a_exp, {17781, 17991});

                pbuffer.scale(2.0 * a_exp, {18621, 18936});

                pbuffer.scale(2.0 * a_exp, {19881, 20322});

                pbuffer.scale(2.0 * a_exp, {21645, 22233});

                pbuffer.scale(2.0 * a_exp, {24417, 24697});

                pbuffer.scale(2.0 * a_exp, {25257, 25677});

                pbuffer.scale(2.0 * a_exp, {26517, 27105});

                pbuffer.scale(2.0 * a_exp, {28281, 29065});

                pbuffer.scale(4.0 * a_exp * a_exp, {30849, 31209});

                pbuffer.scale(4.0 * a_exp * a_exp, {31569, 32109});

                pbuffer.scale(4.0 * a_exp * a_exp, {32649, 33405});

                pbuffer.scale(4.0 * a_exp * a_exp, {34161, 35169});

                pbuffer.scale(4.0 * a_exp * a_exp, {36177, 36627});

                pbuffer.scale(4.0 * a_exp * a_exp, {36627, 37302});

                pbuffer.scale(4.0 * a_exp * a_exp, {37302, 38247});

                pbuffer.scale(4.0 * a_exp * a_exp, {38247, 39507});

                t2cfunc::reduce(cbuffer, 7326, pbuffer, 6741, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7426, pbuffer, 7341, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7576, pbuffer, 8241, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7786, pbuffer, 9501, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8066, pbuffer, 11706, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8216, pbuffer, 12456, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8441, pbuffer, 13581, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8756, pbuffer, 15156, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9176, pbuffer, 17781, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9386, pbuffer, 18621, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9701, pbuffer, 19881, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10142, pbuffer, 21645, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10730, pbuffer, 24417, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11010, pbuffer, 25257, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11430, pbuffer, 26517, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12018, pbuffer, 28281, 784, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12802, pbuffer, 30849, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13162, pbuffer, 31569, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13702, pbuffer, 32649, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14458, pbuffer, 34161, 1008, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15466, pbuffer, 36177, 450, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15916, pbuffer, 36627, 675, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16591, pbuffer, 37302, 945, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 17536, pbuffer, 38247, 1260, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 100, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 300, cbuffer, 100, 250, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 750, cbuffer, 250, 460, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 1380, 0, 300, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 1980, 300, 750, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 2880, 1380, 1980, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 3880, cbuffer, 740, 890, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 4330, cbuffer, 890, 1115, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 5005, cbuffer, 1115, 1430, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 5950, 3880, 4330, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 6850, 4330, 5005, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 8200, 5950, 6850, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9700, cbuffer, 1850, 1950, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 10000, cbuffer, 1950, 2100, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 10450, cbuffer, 2100, 2310, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 11080, 9700, 10000, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 11680, 10000, 10450, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 12580, 11080, 11680, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 13580, cbuffer, 2590, 2740, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 14030, cbuffer, 2740, 2965, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 14705, cbuffer, 2965, 3280, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 15650, 13580, 14030, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 16550, 14030, 14705, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 17900, 15650, 16550, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 19400, cbuffer, 3700, 3910, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 20030, cbuffer, 3910, 4225, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 20975, cbuffer, 4225, 4666, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 22298, 19400, 20030, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 23558, 20030, 20975, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 25448, 22298, 23558, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 27548, cbuffer, 5254, 5534, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 28388, cbuffer, 5534, 5954, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 29648, cbuffer, 5954, 6542, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 31412, 27548, 28388, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 33092, 28388, 29648, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 35612, 31412, 33092, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 38412, cbuffer, 7326, 7426, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 38712, cbuffer, 7426, 7576, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 39162, cbuffer, 7576, 7786, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 39792, 38412, 38712, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 40392, 38712, 39162, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 41292, 39792, 40392, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 42292, cbuffer, 8066, 8216, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 42742, cbuffer, 8216, 8441, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 43417, cbuffer, 8441, 8756, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 44362, 42292, 42742, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 45262, 42742, 43417, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 46612, 44362, 45262, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 48112, cbuffer, 9176, 9386, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 48742, cbuffer, 9386, 9701, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 49687, cbuffer, 9701, 10142, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 51010, 48112, 48742, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 52270, 48742, 49687, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 54160, 51010, 52270, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 56260, cbuffer, 10730, 11010, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 57100, cbuffer, 11010, 11430, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 58360, cbuffer, 11430, 12018, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 60124, 56260, 57100, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 61804, 57100, 58360, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 64324, 60124, 61804, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 67124, cbuffer, 12802, 13162, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 68204, cbuffer, 13162, 13702, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 69824, cbuffer, 13702, 14458, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 72092, 67124, 68204, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 74252, 68204, 69824, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 77492, 72092, 74252, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 81092, cbuffer, 15466, 15916, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 82442, cbuffer, 15916, 16591, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 84467, cbuffer, 16591, 17536, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 87302, 81092, 82442, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 90002, 82442, 84467, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 94052, 87302, 90002, cfactors, 6, 0, 8);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 2880, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 3430, ckbuffer, 8200, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 158368, ckbuffer, 12580, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 158858, ckbuffer, 17900, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 159593, ckbuffer, 25448, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 160622, ckbuffer, 35612, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 168756, ckbuffer, 41292, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 169246, ckbuffer, 46612, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 169981, ckbuffer, 54160, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 171010, ckbuffer, 64324, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 172382, ckbuffer, 77492, 0, 7);

            t4cfunc::ket_transform<3, 3>(skbuffer, 174146, ckbuffer, 94052, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 22981, 0, 3430, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 161994, 158368, 158858, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 163464, 158858, 159593, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 165669, 159593, 160622, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 176351, 168756, 169246, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 177821, 169246, 169981, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 180026, 169981, 171010, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pixx(skbuffer, 183113, 171010, 172382, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pkxx(skbuffer, 187229, 172382, 174146, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 192521, 176351, 177821, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dgxx(skbuffer, 195461, 177821, 180026, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dhxx(skbuffer, 199871, 180026, 183113, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dixx(skbuffer, 206045, 183113, 187229, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 24451, 0, 161994, 163464, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 37681, 3430, 163464, 165669, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 76048, 22981, 24451, 37681, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 490, 158368, 192521, 3, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 4165, 158858, 195461, 4, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 8575, 159593, 199871, 5, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 14749, 160622, 206045, 6, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_pfxx(skbuffer, 28861, 161994, 490, 4165, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_pgxx(skbuffer, 44296, 163464, 4165, 8575, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_phxx(skbuffer, 57526, 165669, 8575, 14749, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dfxx(skbuffer, 84868, 24451, 28861, 44296, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dgxx(skbuffer, 102508, 37681, 44296, 57526, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_ffxx(skbuffer, 128968, 76048, 84868, 102508, r_ab, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 128968, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2401, skbuffer, 133868, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 4802, skbuffer, 138768, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 7203, skbuffer, 143668, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 9604, skbuffer, 148568, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 12005, skbuffer, 153468, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecFFFF_hpp */
