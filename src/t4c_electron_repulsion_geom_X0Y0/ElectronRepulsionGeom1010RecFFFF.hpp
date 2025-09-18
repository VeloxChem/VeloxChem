#ifndef ElectronRepulsionGeom1010RecFFFF_hpp
#define ElectronRepulsionGeom1010RecFFFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXDF.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXFF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPH.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSH.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSI.hpp"
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPHXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSIXX.hpp"
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
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSI.hpp"
#include "ElectronRepulsionPrimRecSHSK.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
#include "ElectronRepulsionPrimRecSISI.hpp"
#include "ElectronRepulsionPrimRecSISK.hpp"
#include "ElectronRepulsionPrimRecSKSF.hpp"
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
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FF|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_ffff(T& distributor,
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

    CSimdArray<double> pbuffer(38233, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(24336, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(237276, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(246519, 1);

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

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 882, 420, 441, 630, 658, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 918, 441, 462, 658, 686, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 954, 462, 483, 686, 714, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 990, 483, 504, 714, 742, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1026, 504, 525, 742, 770, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1062, 525, 546, 770, 798, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1098, 546, 567, 798, 826, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1134, 567, 588, 826, 854, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1170, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1173, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1176, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1179, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1182, 3, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1191, 4, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1200, 5, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1209, 6, 30, 33, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1218, 7, 33, 36, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1227, 21, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1245, 24, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1263, 27, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1281, 30, 81, 87, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1299, 33, 87, 93, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1317, 36, 93, 99, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1335, 63, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1365, 69, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1395, 75, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1425, 81, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1455, 87, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1485, 93, 185, 195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1515, 99, 195, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1545, 145, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1590, 155, 270, 285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1635, 165, 285, 300, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1680, 175, 300, 315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1725, 185, 315, 330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1770, 195, 330, 345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1815, 205, 345, 360, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1860, 270, 420, 441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1923, 285, 441, 462, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1986, 300, 462, 483, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2049, 315, 483, 504, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2112, 330, 504, 525, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2175, 345, 525, 546, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 2238, 360, 546, 567, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2301, 441, 630, 658, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2385, 462, 658, 686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2469, 483, 686, 714, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2553, 504, 714, 742, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2637, 525, 742, 770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2721, 546, 770, 798, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2805, 567, 798, 826, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2889, 658, 882, 918, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2997, 686, 918, 954, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3105, 714, 954, 990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3213, 742, 990, 1026, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3321, 770, 1026, 1062, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3429, 798, 1062, 1098, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3537, 826, 1098, 1134, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 3645, 3, 4, 1170, 1173, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 3651, 4, 5, 1173, 1176, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 3657, 5, 6, 1176, 1179, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3663, 21, 24, 1170, 1182, 1191, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3681, 24, 27, 1173, 1191, 1200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3699, 27, 30, 1176, 1200, 1209, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3717, 30, 33, 1179, 1209, 1218, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3735, 63, 69, 1182, 1227, 1245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3771, 69, 75, 1191, 1245, 1263, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3807, 75, 81, 1200, 1263, 1281, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3843, 81, 87, 1209, 1281, 1299, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3879, 87, 93, 1218, 1299, 1317, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3915, 135, 145, 1227, 1335, 1365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3975, 145, 155, 1245, 1365, 1395, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 4035, 155, 165, 1263, 1395, 1425, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 4095, 165, 175, 1281, 1425, 1455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 4155, 175, 185, 1299, 1455, 1485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 4215, 185, 195, 1317, 1485, 1515, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4275, 255, 270, 1365, 1545, 1590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4365, 270, 285, 1395, 1590, 1635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4455, 285, 300, 1425, 1635, 1680, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4545, 300, 315, 1455, 1680, 1725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4635, 315, 330, 1485, 1725, 1770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4725, 330, 345, 1515, 1770, 1815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4815, 420, 441, 1590, 1860, 1923, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4941, 441, 462, 1635, 1923, 1986, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5067, 462, 483, 1680, 1986, 2049, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5193, 483, 504, 1725, 2049, 2112, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5319, 504, 525, 1770, 2112, 2175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 5445, 525, 546, 1815, 2175, 2238, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5571, 630, 658, 1923, 2301, 2385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5739, 658, 686, 1986, 2385, 2469, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5907, 686, 714, 2049, 2469, 2553, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6075, 714, 742, 2112, 2553, 2637, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6243, 742, 770, 2175, 2637, 2721, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 6411, 770, 798, 2238, 2721, 2805, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6579, 882, 918, 2385, 2889, 2997, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6795, 918, 954, 2469, 2997, 3105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7011, 954, 990, 2553, 3105, 3213, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7227, 990, 1026, 2637, 3213, 3321, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7443, 1026, 1062, 2721, 3321, 3429, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 7659, 1062, 1098, 2805, 3429, 3537, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 7875, 1170, 1173, 3645, 3651, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 7885, 1173, 1176, 3651, 3657, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 7895, 1182, 1191, 3645, 3663, 3681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 7925, 1191, 1200, 3651, 3681, 3699, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 7955, 1200, 1209, 3657, 3699, 3717, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 7985, 1227, 1245, 3663, 3735, 3771, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 8045, 1245, 1263, 3681, 3771, 3807, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 8105, 1263, 1281, 3699, 3807, 3843, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 8165, 1281, 1299, 3717, 3843, 3879, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 8225, 1335, 1365, 3735, 3915, 3975, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 8325, 1365, 1395, 3771, 3975, 4035, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 8425, 1395, 1425, 3807, 4035, 4095, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 8525, 1425, 1455, 3843, 4095, 4155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 8625, 1455, 1485, 3879, 4155, 4215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8725, 1545, 1590, 3975, 4275, 4365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 8875, 1590, 1635, 4035, 4365, 4455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 9025, 1635, 1680, 4095, 4455, 4545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 9175, 1680, 1725, 4155, 4545, 4635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 9325, 1725, 1770, 4215, 4635, 4725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9475, 1860, 1923, 4365, 4815, 4941, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9685, 1923, 1986, 4455, 4941, 5067, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 9895, 1986, 2049, 4545, 5067, 5193, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 10105, 2049, 2112, 4635, 5193, 5319, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 10315, 2112, 2175, 4725, 5319, 5445, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10525, 2301, 2385, 4941, 5571, 5739, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 10805, 2385, 2469, 5067, 5739, 5907, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 11085, 2469, 2553, 5193, 5907, 6075, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 11365, 2553, 2637, 5319, 6075, 6243, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 11645, 2637, 2721, 5445, 6243, 6411, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 11925, 2889, 2997, 5739, 6579, 6795, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 12285, 2997, 3105, 5907, 6795, 7011, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 12645, 3105, 3213, 6075, 7011, 7227, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 13005, 3213, 3321, 6243, 7227, 7443, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 13365, 3321, 3429, 6411, 7443, 7659, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 13725, 3645, 3651, 7875, 7885, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 13740, 3663, 3681, 7875, 7895, 7925, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 13785, 3681, 3699, 7885, 7925, 7955, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 13830, 3735, 3771, 7895, 7985, 8045, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 13920, 3771, 3807, 7925, 8045, 8105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 14010, 3807, 3843, 7955, 8105, 8165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 14100, 3915, 3975, 7985, 8225, 8325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 14250, 3975, 4035, 8045, 8325, 8425, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 14400, 4035, 4095, 8105, 8425, 8525, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 14550, 4095, 4155, 8165, 8525, 8625, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 14700, 4275, 4365, 8325, 8725, 8875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 14925, 4365, 4455, 8425, 8875, 9025, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 15150, 4455, 4545, 8525, 9025, 9175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 15375, 4545, 4635, 8625, 9175, 9325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 15600, 4815, 4941, 8875, 9475, 9685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 15915, 4941, 5067, 9025, 9685, 9895, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 16230, 5067, 5193, 9175, 9895, 10105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 16545, 5193, 5319, 9325, 10105, 10315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 16860, 5571, 5739, 9685, 10525, 10805, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 17280, 5739, 5907, 9895, 10805, 11085, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 17700, 5907, 6075, 10105, 11085, 11365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 18120, 6075, 6243, 10315, 11365, 11645, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 18540, 6579, 6795, 10805, 11925, 12285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 19080, 6795, 7011, 11085, 12285, 12645, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 19620, 7011, 7227, 11365, 12645, 13005, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 20160, 7227, 7443, 11645, 13005, 13365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 20700, 7895, 7925, 13725, 13740, 13785, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 20763, 7985, 8045, 13740, 13830, 13920, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 20889, 8045, 8105, 13785, 13920, 14010, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 21015, 8225, 8325, 13830, 14100, 14250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 21225, 8325, 8425, 13920, 14250, 14400, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 21435, 8425, 8525, 14010, 14400, 14550, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 21645, 8725, 8875, 14250, 14700, 14925, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 21960, 8875, 9025, 14400, 14925, 15150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 22275, 9025, 9175, 14550, 15150, 15375, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 22590, 9475, 9685, 14925, 15600, 15915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 23031, 9685, 9895, 15150, 15915, 16230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 23472, 9895, 10105, 15375, 16230, 16545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 23913, 10525, 10805, 15915, 16860, 17280, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 24501, 10805, 11085, 16230, 17280, 17700, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 25089, 11085, 11365, 16545, 17700, 18120, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsk(pbuffer, 25677, 11925, 12285, 17280, 18540, 19080, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsk(pbuffer, 26433, 12285, 12645, 17700, 19080, 19620, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsk(pbuffer, 27189, 12645, 13005, 18120, 19620, 20160, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 27945, 13830, 13920, 20700, 20763, 20889, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 28113, 14100, 14250, 20763, 21015, 21225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 28393, 14250, 14400, 20889, 21225, 21435, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 28673, 14700, 14925, 21225, 21645, 21960, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 29093, 14925, 15150, 21435, 21960, 22275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 29513, 15600, 15915, 21960, 22590, 23031, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 30101, 15915, 16230, 22275, 23031, 23472, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 30689, 16860, 17280, 23031, 23913, 24501, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 31473, 17280, 17700, 23472, 24501, 25089, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisk(pbuffer, 32257, 18540, 19080, 24501, 25677, 26433, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisk(pbuffer, 33265, 19080, 19620, 25089, 26433, 27189, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 34273, 21015, 21225, 27945, 28113, 28393, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 34633, 21645, 21960, 28393, 28673, 29093, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksh(pbuffer, 35173, 22590, 23031, 29093, 29513, 30101, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksi(pbuffer, 35929, 23913, 24501, 30101, 30689, 31473, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksk(pbuffer, 36937, 25677, 26433, 31473, 32257, 33265, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 8225, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 8725, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 250, pbuffer, 9475, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 460, pbuffer, 14100, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 610, pbuffer, 14700, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 835, pbuffer, 15600, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1150, pbuffer, 21015, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1360, pbuffer, 21645, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1675, pbuffer, 22590, 441, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {8225, 8325});

                pbuffer.scale(2.0 * a_exp, {8725, 8875});

                pbuffer.scale(2.0 * a_exp, {9475, 9685});

                pbuffer.scale(2.0 * a_exp, {14100, 14250});

                pbuffer.scale(2.0 * a_exp, {14700, 14925});

                pbuffer.scale(2.0 * a_exp, {15600, 15915});

                pbuffer.scale(2.0 * a_exp, {21015, 21225});

                pbuffer.scale(2.0 * a_exp, {21645, 21960});

                pbuffer.scale(2.0 * a_exp, {22590, 23031});

                pbuffer.scale(2.0 * a_exp, {28113, 28393});

                pbuffer.scale(2.0 * a_exp, {28673, 29093});

                pbuffer.scale(2.0 * a_exp, {29513, 30101});

                pbuffer.scale(2.0 * a_exp, {34273, 34633});

                pbuffer.scale(2.0 * a_exp, {34633, 35173});

                pbuffer.scale(2.0 * a_exp, {35173, 35929});

                t2cfunc::reduce(cbuffer, 7176, pbuffer, 8225, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7276, pbuffer, 8725, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7426, pbuffer, 9475, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7636, pbuffer, 14100, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7786, pbuffer, 14700, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8011, pbuffer, 15600, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8326, pbuffer, 21015, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8536, pbuffer, 21645, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8851, pbuffer, 22590, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9292, pbuffer, 28113, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9572, pbuffer, 28673, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9992, pbuffer, 29513, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10580, pbuffer, 34273, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10940, pbuffer, 34633, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11480, pbuffer, 35173, 756, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {8225, 8325});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {8725, 8875});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {9475, 9685});

                pbuffer.scale(pfactors, 0, 2.0, {10525, 10805});

                pbuffer.scale(pfactors, 0, 2.0, {11925, 12285});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {14100, 14250});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {14700, 14925});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {15600, 15915});

                pbuffer.scale(pfactors, 0, 2.0, {16860, 17280});

                pbuffer.scale(pfactors, 0, 2.0, {18540, 19080});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {21015, 21225});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {21645, 21960});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {22590, 23031});

                pbuffer.scale(pfactors, 0, 2.0, {23913, 24501});

                pbuffer.scale(pfactors, 0, 2.0, {25677, 26433});

                t2cfunc::reduce(cbuffer, 2116, pbuffer, 8225, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2216, pbuffer, 8725, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2366, pbuffer, 9475, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2576, pbuffer, 10525, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2856, pbuffer, 11925, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3216, pbuffer, 14100, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3366, pbuffer, 14700, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3591, pbuffer, 15600, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3906, pbuffer, 16860, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4326, pbuffer, 18540, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4866, pbuffer, 21015, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5076, pbuffer, 21645, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5391, pbuffer, 22590, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5832, pbuffer, 23913, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6420, pbuffer, 25677, 756, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {8225, 8325});

                pbuffer.scale(2.0 * a_exp, {8725, 8875});

                pbuffer.scale(2.0 * a_exp, {9475, 9685});

                pbuffer.scale(2.0 * a_exp, {10525, 10805});

                pbuffer.scale(2.0 * a_exp, {11925, 12285});

                pbuffer.scale(2.0 * a_exp, {14100, 14250});

                pbuffer.scale(2.0 * a_exp, {14700, 14925});

                pbuffer.scale(2.0 * a_exp, {15600, 15915});

                pbuffer.scale(2.0 * a_exp, {16860, 17280});

                pbuffer.scale(2.0 * a_exp, {18540, 19080});

                pbuffer.scale(2.0 * a_exp, {21015, 21225});

                pbuffer.scale(2.0 * a_exp, {21645, 21960});

                pbuffer.scale(2.0 * a_exp, {22590, 23031});

                pbuffer.scale(2.0 * a_exp, {23913, 24501});

                pbuffer.scale(2.0 * a_exp, {25677, 26433});

                pbuffer.scale(pfactors, 0, 2.0, {28113, 28393});

                pbuffer.scale(pfactors, 0, 2.0, {28673, 29093});

                pbuffer.scale(pfactors, 0, 2.0, {29513, 30101});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {30689, 31473});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {32257, 33265});

                pbuffer.scale(pfactors, 0, 2.0, {34273, 34633});

                pbuffer.scale(pfactors, 0, 2.0, {34633, 35173});

                pbuffer.scale(pfactors, 0, 2.0, {35173, 35929});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {35929, 36937});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {36937, 38233});

                t2cfunc::reduce(cbuffer, 12236, pbuffer, 8225, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12336, pbuffer, 8725, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12486, pbuffer, 9475, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12696, pbuffer, 10525, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12976, pbuffer, 11925, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13336, pbuffer, 14100, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13486, pbuffer, 14700, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13711, pbuffer, 15600, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14026, pbuffer, 16860, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14446, pbuffer, 18540, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14986, pbuffer, 21015, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15196, pbuffer, 21645, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15511, pbuffer, 22590, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15952, pbuffer, 23913, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16540, pbuffer, 25677, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 17296, pbuffer, 28113, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 17576, pbuffer, 28673, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 17996, pbuffer, 29513, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18584, pbuffer, 30689, 784, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 19368, pbuffer, 32257, 1008, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 20376, pbuffer, 34273, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 20736, pbuffer, 34633, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 21276, pbuffer, 35173, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 22032, pbuffer, 35929, 1008, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 23040, pbuffer, 36937, 1296, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2220, cbuffer, 0, 100, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 3420, cbuffer, 100, 250, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 7110, 2220, 3420, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 18540, cbuffer, 460, 610, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 20340, cbuffer, 610, 835, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 25875, 18540, 20340, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 42687, cbuffer, 1150, 1360, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 45207, cbuffer, 1360, 1675, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 52956, 42687, 45207, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 72186, cbuffer, 7176, 7276, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 73386, cbuffer, 7276, 7426, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 77076, 72186, 73386, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 88506, cbuffer, 7636, 7786, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 90306, cbuffer, 7786, 8011, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 95841, 88506, 90306, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 112653, cbuffer, 8326, 8536, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 115173, cbuffer, 8536, 8851, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 122922, 112653, 115173, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 146148, cbuffer, 9292, 9572, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 149508, cbuffer, 9572, 9992, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 159840, 146148, 149508, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 190512, cbuffer, 10580, 10940, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 194832, cbuffer, 10940, 11480, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 208116, 190512, 194832, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 0, cbuffer, 2116, 2216, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 300, cbuffer, 2216, 2366, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 750, cbuffer, 2366, 2576, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 1380, cbuffer, 2576, 2856, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2520, cbuffer, 0, 0, 300, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 3870, cbuffer, 100, 300, 750, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 5220, cbuffer, 250, 750, 1380, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 7710, 2220, 2520, 3870, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 9510, 3420, 3870, 5220, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 12210, 7110, 7710, 9510, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 15210, cbuffer, 3216, 3366, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 15660, cbuffer, 3366, 3591, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 16335, cbuffer, 3591, 3906, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 17280, cbuffer, 3906, 4326, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 18990, cbuffer, 460, 15210, 15660, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 21015, cbuffer, 610, 15660, 16335, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 23040, cbuffer, 835, 16335, 17280, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 26775, 18540, 18990, 21015, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 29475, 20340, 21015, 23040, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 33525, 25875, 26775, 29475, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 38025, cbuffer, 4866, 5076, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 38655, cbuffer, 5076, 5391, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 39600, cbuffer, 5391, 5832, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 40923, cbuffer, 5832, 6420, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 43317, cbuffer, 1150, 38025, 38655, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 46152, cbuffer, 1360, 38655, 39600, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 48987, cbuffer, 1675, 39600, 40923, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 54216, 42687, 43317, 46152, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 57996, 45207, 46152, 48987, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 63666, 52956, 54216, 57996, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 69966, cbuffer, 12236, 12336, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 70266, cbuffer, 12336, 12486, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 70716, cbuffer, 12486, 12696, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 71346, cbuffer, 12696, 12976, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 72486, cbuffer, 7176, 69966, 70266, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 73836, cbuffer, 7276, 70266, 70716, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 75186, cbuffer, 7426, 70716, 71346, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 77676, 72186, 72486, 73836, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 79476, 73386, 73836, 75186, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 82176, 77076, 77676, 79476, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 85176, cbuffer, 13336, 13486, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 85626, cbuffer, 13486, 13711, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 86301, cbuffer, 13711, 14026, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 87246, cbuffer, 14026, 14446, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 88956, cbuffer, 7636, 85176, 85626, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 90981, cbuffer, 7786, 85626, 86301, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 93006, cbuffer, 8011, 86301, 87246, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 96741, 88506, 88956, 90981, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 99441, 90306, 90981, 93006, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 103491, 95841, 96741, 99441, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 107991, cbuffer, 14986, 15196, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 108621, cbuffer, 15196, 15511, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 109566, cbuffer, 15511, 15952, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 110889, cbuffer, 15952, 16540, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 113283, cbuffer, 8326, 107991, 108621, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 116118, cbuffer, 8536, 108621, 109566, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 118953, cbuffer, 8851, 109566, 110889, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 124182, 112653, 113283, 116118, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 127962, 115173, 116118, 118953, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 133632, 122922, 124182, 127962, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 139932, cbuffer, 17296, 17576, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 140772, cbuffer, 17576, 17996, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 142032, cbuffer, 17996, 18584, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 143796, cbuffer, 18584, 19368, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 146988, cbuffer, 9292, 139932, 140772, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 150768, cbuffer, 9572, 140772, 142032, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 154548, cbuffer, 9992, 142032, 143796, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 161520, 146148, 146988, 150768, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 166560, 149508, 150768, 154548, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 174120, 159840, 161520, 166560, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 182520, cbuffer, 20376, 20736, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 183600, cbuffer, 20736, 21276, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 185220, cbuffer, 21276, 22032, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 187488, cbuffer, 22032, 23040, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 191592, cbuffer, 10580, 182520, 183600, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 196452, cbuffer, 10940, 183600, 185220, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 201312, cbuffer, 11480, 185220, 187488, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 210276, 190512, 191592, 196452, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 216756, 194832, 196452, 201312, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 226476, 208116, 210276, 216756, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 12210, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 490, ckbuffer, 13210, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 980, ckbuffer, 14210, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 5880, ckbuffer, 33525, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 6615, ckbuffer, 35025, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 7350, ckbuffer, 36525, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 14700, ckbuffer, 63666, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 15729, ckbuffer, 65766, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 16758, ckbuffer, 67866, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 230349, ckbuffer, 82176, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 230839, ckbuffer, 83176, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 231329, ckbuffer, 84176, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 231819, ckbuffer, 103491, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 232554, ckbuffer, 104991, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 233289, ckbuffer, 106491, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 234024, ckbuffer, 133632, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 235053, ckbuffer, 135732, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 236082, ckbuffer, 137832, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 237111, ckbuffer, 174120, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 238483, ckbuffer, 176920, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 239855, ckbuffer, 179720, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 241227, ckbuffer, 226476, 0, 7);

            t4cfunc::ket_transform<3, 3>(skbuffer, 242991, ckbuffer, 230076, 0, 7);

            t4cfunc::ket_transform<3, 3>(skbuffer, 244755, ckbuffer, 233676, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 39396, 0, 5880, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 40866, 490, 6615, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 42336, 980, 7350, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 57036, 5880, 14700, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 59241, 6615, 15729, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 61446, 7350, 16758, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 111279, 39396, 57036, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 114219, 40866, 59241, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 117159, 42336, 61446, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 1470, 230349, 231819, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 8085, 231819, 234024, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 17787, 234024, 237111, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sixx(skbuffer, 27048, 237111, 241227, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 43806, 0, 1470, 8085, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 63651, 5880, 8085, 17787, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_phxx(skbuffer, 83496, 14700, 17787, 27048, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 120099, 39396, 43806, 63651, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dgxx(skbuffer, 146559, 57036, 63651, 83496, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ffxx(skbuffer, 186249, 111279, 120099, 146559, r_ab, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 186249, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2401, skbuffer, 191149, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 4802, skbuffer, 196049, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 7203, skbuffer, 200949, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 9604, skbuffer, 205849, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 12005, skbuffer, 210749, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 14406, skbuffer, 215649, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 16807, skbuffer, 220549, 3, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 19208, skbuffer, 225449, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFFFF_hpp */
