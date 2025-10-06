#ifndef ElectronRepulsionGeom1010RecFFFD_hpp
#define ElectronRepulsionGeom1010RecFFFD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXDD.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXFD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSH.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FF|1/|r-r'||FD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fffd(T& distributor,
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

    CSimdArray<double> pbuffer(27287, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(17316, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(153972, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(176085, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(15435, 1);

    // setup Boys fuction data

    const CBoysFunc<13> bf_table;

    CSimdArray<double> bf_data(15, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 14, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 14, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 14, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 14);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 29, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 32, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 35, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 38, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 41, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 44, 10, 11, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 47, 11, 12, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 50, 12, 13, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 0, 1, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 1, 2, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 65, 2, 3, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 71, 3, 4, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 77, 4, 5, 26, 29, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 83, 5, 6, 29, 32, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 89, 6, 7, 32, 35, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 95, 7, 8, 35, 38, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 101, 8, 9, 38, 41, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 107, 9, 10, 41, 44, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 113, 10, 11, 44, 47, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 119, 11, 12, 47, 50, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 14, 17, 53, 59, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 17, 20, 59, 65, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 20, 23, 65, 71, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 23, 26, 71, 77, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 26, 29, 77, 83, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 175, 29, 32, 83, 89, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 185, 32, 35, 89, 95, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 195, 35, 38, 95, 101, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 205, 38, 41, 101, 107, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 215, 41, 44, 107, 113, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 225, 44, 47, 113, 119, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 235, 53, 59, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 250, 59, 65, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 265, 65, 71, 145, 155, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 280, 71, 77, 155, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 295, 77, 83, 165, 175, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 310, 83, 89, 175, 185, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 325, 89, 95, 185, 195, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 340, 95, 101, 195, 205, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 355, 101, 107, 205, 215, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 370, 107, 113, 215, 225, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 385, 125, 135, 235, 250, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 406, 135, 145, 250, 265, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 427, 145, 155, 265, 280, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 448, 155, 165, 280, 295, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 469, 165, 175, 295, 310, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 490, 175, 185, 310, 325, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 511, 185, 195, 325, 340, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 532, 195, 205, 340, 355, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 553, 205, 215, 355, 370, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 574, 235, 250, 385, 406, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 602, 250, 265, 406, 427, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 630, 265, 280, 427, 448, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 658, 280, 295, 448, 469, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 686, 295, 310, 469, 490, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 714, 310, 325, 490, 511, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 742, 325, 340, 511, 532, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 770, 340, 355, 532, 553, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 798, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 801, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 804, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 807, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 810, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 813, 2, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 822, 3, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 831, 4, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 840, 5, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 849, 6, 29, 32, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 858, 7, 32, 35, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 867, 17, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 885, 20, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 903, 23, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 921, 26, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 939, 29, 77, 83, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 957, 32, 83, 89, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 975, 35, 89, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 993, 59, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1023, 65, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1053, 71, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1083, 77, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1113, 83, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1143, 89, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1173, 95, 185, 195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1203, 135, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1248, 145, 250, 265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1293, 155, 265, 280, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1338, 165, 280, 295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1383, 175, 295, 310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1428, 185, 310, 325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1473, 195, 325, 340, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1518, 250, 385, 406, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1581, 265, 406, 427, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1644, 280, 427, 448, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1707, 295, 448, 469, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1770, 310, 469, 490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1833, 325, 490, 511, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1896, 340, 511, 532, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1959, 406, 574, 602, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2043, 427, 602, 630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2127, 448, 630, 658, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2211, 469, 658, 686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2295, 490, 686, 714, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2379, 511, 714, 742, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2463, 532, 742, 770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2547, 2, 3, 798, 801, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2553, 3, 4, 801, 804, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2559, 4, 5, 804, 807, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2565, 5, 6, 807, 810, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2571, 17, 20, 798, 813, 822, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2589, 20, 23, 801, 822, 831, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2607, 23, 26, 804, 831, 840, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2625, 26, 29, 807, 840, 849, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2643, 29, 32, 810, 849, 858, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2661, 53, 59, 813, 867, 885, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2697, 59, 65, 822, 885, 903, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2733, 65, 71, 831, 903, 921, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2769, 71, 77, 840, 921, 939, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2805, 77, 83, 849, 939, 957, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2841, 83, 89, 858, 957, 975, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2877, 125, 135, 885, 993, 1023, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2937, 135, 145, 903, 1023, 1053, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2997, 145, 155, 921, 1053, 1083, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3057, 155, 165, 939, 1083, 1113, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3117, 165, 175, 957, 1113, 1143, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3177, 175, 185, 975, 1143, 1173, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3237, 235, 250, 1023, 1203, 1248, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3327, 250, 265, 1053, 1248, 1293, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3417, 265, 280, 1083, 1293, 1338, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3507, 280, 295, 1113, 1338, 1383, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3597, 295, 310, 1143, 1383, 1428, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3687, 310, 325, 1173, 1428, 1473, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3777, 385, 406, 1248, 1518, 1581, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3903, 406, 427, 1293, 1581, 1644, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4029, 427, 448, 1338, 1644, 1707, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4155, 448, 469, 1383, 1707, 1770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4281, 469, 490, 1428, 1770, 1833, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4407, 490, 511, 1473, 1833, 1896, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4533, 574, 602, 1581, 1959, 2043, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4701, 602, 630, 1644, 2043, 2127, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4869, 630, 658, 1707, 2127, 2211, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5037, 658, 686, 1770, 2211, 2295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5205, 686, 714, 1833, 2295, 2379, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5373, 714, 742, 1896, 2379, 2463, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 5541, 798, 801, 2547, 2553, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 5551, 801, 804, 2553, 2559, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 5561, 804, 807, 2559, 2565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 5571, 813, 822, 2547, 2571, 2589, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 5601, 822, 831, 2553, 2589, 2607, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 5631, 831, 840, 2559, 2607, 2625, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 5661, 840, 849, 2565, 2625, 2643, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 5691, 867, 885, 2571, 2661, 2697, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 5751, 885, 903, 2589, 2697, 2733, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 5811, 903, 921, 2607, 2733, 2769, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 5871, 921, 939, 2625, 2769, 2805, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 5931, 939, 957, 2643, 2805, 2841, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 5991, 993, 1023, 2697, 2877, 2937, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 6091, 1023, 1053, 2733, 2937, 2997, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 6191, 1053, 1083, 2769, 2997, 3057, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 6291, 1083, 1113, 2805, 3057, 3117, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 6391, 1113, 1143, 2841, 3117, 3177, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 6491, 1203, 1248, 2937, 3237, 3327, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 6641, 1248, 1293, 2997, 3327, 3417, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 6791, 1293, 1338, 3057, 3417, 3507, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 6941, 1338, 1383, 3117, 3507, 3597, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 7091, 1383, 1428, 3177, 3597, 3687, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 7241, 1518, 1581, 3327, 3777, 3903, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 7451, 1581, 1644, 3417, 3903, 4029, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 7661, 1644, 1707, 3507, 4029, 4155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 7871, 1707, 1770, 3597, 4155, 4281, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 8081, 1770, 1833, 3687, 4281, 4407, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 8291, 1959, 2043, 3903, 4533, 4701, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 8571, 2043, 2127, 4029, 4701, 4869, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 8851, 2127, 2211, 4155, 4869, 5037, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 9131, 2211, 2295, 4281, 5037, 5205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 9411, 2295, 2379, 4407, 5205, 5373, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 9691, 2547, 2553, 5541, 5551, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 9706, 2553, 2559, 5551, 5561, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 9721, 2571, 2589, 5541, 5571, 5601, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 9766, 2589, 2607, 5551, 5601, 5631, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 9811, 2607, 2625, 5561, 5631, 5661, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 9856, 2661, 2697, 5571, 5691, 5751, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 9946, 2697, 2733, 5601, 5751, 5811, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 10036, 2733, 2769, 5631, 5811, 5871, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 10126, 2769, 2805, 5661, 5871, 5931, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 10216, 2877, 2937, 5751, 5991, 6091, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 10366, 2937, 2997, 5811, 6091, 6191, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 10516, 2997, 3057, 5871, 6191, 6291, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 10666, 3057, 3117, 5931, 6291, 6391, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 10816, 3237, 3327, 6091, 6491, 6641, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 11041, 3327, 3417, 6191, 6641, 6791, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 11266, 3417, 3507, 6291, 6791, 6941, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 11491, 3507, 3597, 6391, 6941, 7091, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 11716, 3777, 3903, 6641, 7241, 7451, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 12031, 3903, 4029, 6791, 7451, 7661, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 12346, 4029, 4155, 6941, 7661, 7871, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 12661, 4155, 4281, 7091, 7871, 8081, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 12976, 4533, 4701, 7451, 8291, 8571, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 13396, 4701, 4869, 7661, 8571, 8851, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 13816, 4869, 5037, 7871, 8851, 9131, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 14236, 5037, 5205, 8081, 9131, 9411, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 14656, 5541, 5551, 9691, 9706, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 14677, 5571, 5601, 9691, 9721, 9766, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 14740, 5601, 5631, 9706, 9766, 9811, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 14803, 5691, 5751, 9721, 9856, 9946, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 14929, 5751, 5811, 9766, 9946, 10036, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 15055, 5811, 5871, 9811, 10036, 10126, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 15181, 5991, 6091, 9946, 10216, 10366, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 15391, 6091, 6191, 10036, 10366, 10516, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 15601, 6191, 6291, 10126, 10516, 10666, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 15811, 6491, 6641, 10366, 10816, 11041, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 16126, 6641, 6791, 10516, 11041, 11266, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 16441, 6791, 6941, 10666, 11266, 11491, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 16756, 7241, 7451, 11041, 11716, 12031, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 17197, 7451, 7661, 11266, 12031, 12346, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 17638, 7661, 7871, 11491, 12346, 12661, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 18079, 8291, 8571, 12031, 12976, 13396, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 18667, 8571, 8851, 12346, 13396, 13816, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 19255, 8851, 9131, 12661, 13816, 14236, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 19843, 9721, 9766, 14656, 14677, 14740, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 19927, 9856, 9946, 14677, 14803, 14929, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 20095, 9946, 10036, 14740, 14929, 15055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 20263, 10216, 10366, 14929, 15181, 15391, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 20543, 10366, 10516, 15055, 15391, 15601, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 20823, 10816, 11041, 15391, 15811, 16126, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 21243, 11041, 11266, 15601, 16126, 16441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 21663, 11716, 12031, 16126, 16756, 17197, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 22251, 12031, 12346, 16441, 17197, 17638, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 22839, 12976, 13396, 17197, 18079, 18667, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 23623, 13396, 13816, 17638, 18667, 19255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 24407, 14803, 14929, 19843, 19927, 20095, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 24623, 15181, 15391, 20095, 20263, 20543, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 24983, 15811, 16126, 20543, 20823, 21243, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksh(pbuffer, 25523, 16756, 17197, 21243, 21663, 22251, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksi(pbuffer, 26279, 18079, 18667, 22251, 22839, 23623, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 5691, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 5991, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 6491, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 310, pbuffer, 9856, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 400, pbuffer, 10216, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 550, pbuffer, 10816, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 775, pbuffer, 14803, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 901, pbuffer, 15181, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1111, pbuffer, 15811, 315, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {5691, 5751});

                pbuffer.scale(2.0 * a_exp, {5991, 6091});

                pbuffer.scale(2.0 * a_exp, {6491, 6641});

                pbuffer.scale(2.0 * a_exp, {9856, 9946});

                pbuffer.scale(2.0 * a_exp, {10216, 10366});

                pbuffer.scale(2.0 * a_exp, {10816, 11041});

                pbuffer.scale(2.0 * a_exp, {14803, 14929});

                pbuffer.scale(2.0 * a_exp, {15181, 15391});

                pbuffer.scale(2.0 * a_exp, {15811, 16126});

                pbuffer.scale(2.0 * a_exp, {19927, 20095});

                pbuffer.scale(2.0 * a_exp, {20263, 20543});

                pbuffer.scale(2.0 * a_exp, {20823, 21243});

                pbuffer.scale(2.0 * a_exp, {24407, 24623});

                pbuffer.scale(2.0 * a_exp, {24623, 24983});

                pbuffer.scale(2.0 * a_exp, {24983, 25523});

                t2cfunc::reduce(cbuffer, 5106, pbuffer, 5691, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5166, pbuffer, 5991, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5266, pbuffer, 6491, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5416, pbuffer, 9856, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5506, pbuffer, 10216, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5656, pbuffer, 10816, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5881, pbuffer, 14803, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6007, pbuffer, 15181, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6217, pbuffer, 15811, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6532, pbuffer, 19927, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6700, pbuffer, 20263, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6980, pbuffer, 20823, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7400, pbuffer, 24407, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7616, pbuffer, 24623, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7976, pbuffer, 24983, 540, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {5691, 5751});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {5991, 6091});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {6491, 6641});

                pbuffer.scale(pfactors, 0, 2.0, {7241, 7451});

                pbuffer.scale(pfactors, 0, 2.0, {8291, 8571});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {9856, 9946});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {10216, 10366});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {10816, 11041});

                pbuffer.scale(pfactors, 0, 2.0, {11716, 12031});

                pbuffer.scale(pfactors, 0, 2.0, {12976, 13396});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {14803, 14929});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {15181, 15391});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {15811, 16126});

                pbuffer.scale(pfactors, 0, 2.0, {16756, 17197});

                pbuffer.scale(pfactors, 0, 2.0, {18079, 18667});

                t2cfunc::reduce(cbuffer, 1426, pbuffer, 5691, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1486, pbuffer, 5991, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1586, pbuffer, 6491, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1736, pbuffer, 7241, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1946, pbuffer, 8291, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2226, pbuffer, 9856, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2316, pbuffer, 10216, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2466, pbuffer, 10816, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2691, pbuffer, 11716, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3006, pbuffer, 12976, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3426, pbuffer, 14803, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3552, pbuffer, 15181, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3762, pbuffer, 15811, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4077, pbuffer, 16756, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4518, pbuffer, 18079, 588, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {5691, 5751});

                pbuffer.scale(2.0 * a_exp, {5991, 6091});

                pbuffer.scale(2.0 * a_exp, {6491, 6641});

                pbuffer.scale(2.0 * a_exp, {7241, 7451});

                pbuffer.scale(2.0 * a_exp, {8291, 8571});

                pbuffer.scale(2.0 * a_exp, {9856, 9946});

                pbuffer.scale(2.0 * a_exp, {10216, 10366});

                pbuffer.scale(2.0 * a_exp, {10816, 11041});

                pbuffer.scale(2.0 * a_exp, {11716, 12031});

                pbuffer.scale(2.0 * a_exp, {12976, 13396});

                pbuffer.scale(2.0 * a_exp, {14803, 14929});

                pbuffer.scale(2.0 * a_exp, {15181, 15391});

                pbuffer.scale(2.0 * a_exp, {15811, 16126});

                pbuffer.scale(2.0 * a_exp, {16756, 17197});

                pbuffer.scale(2.0 * a_exp, {18079, 18667});

                pbuffer.scale(pfactors, 0, 2.0, {19927, 20095});

                pbuffer.scale(pfactors, 0, 2.0, {20263, 20543});

                pbuffer.scale(pfactors, 0, 2.0, {20823, 21243});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {21663, 22251});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {22839, 23623});

                pbuffer.scale(pfactors, 0, 2.0, {24407, 24623});

                pbuffer.scale(pfactors, 0, 2.0, {24623, 24983});

                pbuffer.scale(pfactors, 0, 2.0, {24983, 25523});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {25523, 26279});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {26279, 27287});

                t2cfunc::reduce(cbuffer, 8516, pbuffer, 5691, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8576, pbuffer, 5991, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8676, pbuffer, 6491, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8826, pbuffer, 7241, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9036, pbuffer, 8291, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9316, pbuffer, 9856, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9406, pbuffer, 10216, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9556, pbuffer, 10816, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9781, pbuffer, 11716, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10096, pbuffer, 12976, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10516, pbuffer, 14803, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10642, pbuffer, 15181, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10852, pbuffer, 15811, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11167, pbuffer, 16756, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11608, pbuffer, 18079, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12196, pbuffer, 19927, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12364, pbuffer, 20263, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12644, pbuffer, 20823, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13064, pbuffer, 21663, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13652, pbuffer, 22839, 784, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14436, pbuffer, 24407, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14652, pbuffer, 24623, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15012, pbuffer, 24983, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15552, pbuffer, 25523, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16308, pbuffer, 26279, 1008, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1560, cbuffer, 0, 60, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2280, cbuffer, 60, 160, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 4830, 1560, 2280, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 12210, cbuffer, 310, 400, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 13290, cbuffer, 400, 550, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 17115, 12210, 13290, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 27951, cbuffer, 775, 901, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 29463, cbuffer, 901, 1111, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 34818, 27951, 29463, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 46962, cbuffer, 5106, 5166, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 47682, cbuffer, 5166, 5266, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 50232, 46962, 47682, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 57612, cbuffer, 5416, 5506, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 58692, cbuffer, 5506, 5656, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 62517, 57612, 58692, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 73353, cbuffer, 5881, 6007, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 74865, cbuffer, 6007, 6217, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 80220, 73353, 74865, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 95172, cbuffer, 6532, 6700, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 97188, cbuffer, 6700, 6980, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 104328, 95172, 97188, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 124056, cbuffer, 7400, 7616, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 126648, cbuffer, 7616, 7976, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 135828, 124056, 126648, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 0, cbuffer, 1426, 1486, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 180, cbuffer, 1486, 1586, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 480, cbuffer, 1586, 1736, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 930, cbuffer, 1736, 1946, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1740, cbuffer, 0, 0, 180, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2580, cbuffer, 60, 180, 480, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 3480, cbuffer, 160, 480, 930, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 5190, 1560, 1740, 2580, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 6270, 2280, 2580, 3480, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 8070, 4830, 5190, 6270, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9870, cbuffer, 2226, 2316, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 10140, cbuffer, 2316, 2466, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 10590, cbuffer, 2466, 2691, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 11265, cbuffer, 2691, 3006, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 12480, cbuffer, 310, 9870, 10140, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 13740, cbuffer, 400, 10140, 10590, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 15090, cbuffer, 550, 10590, 11265, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 17655, 12210, 12480, 13740, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 19275, 13290, 13740, 15090, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 21975, 17115, 17655, 19275, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 24675, cbuffer, 3426, 3552, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 25053, cbuffer, 3552, 3762, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 25683, cbuffer, 3762, 4077, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 26628, cbuffer, 4077, 4518, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 28329, cbuffer, 775, 24675, 25053, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 30093, cbuffer, 901, 25053, 25683, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 31983, cbuffer, 1111, 25683, 26628, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 35574, 27951, 28329, 30093, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 37842, 29463, 30093, 31983, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 41622, 34818, 35574, 37842, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 45402, cbuffer, 8516, 8576, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 45582, cbuffer, 8576, 8676, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 45882, cbuffer, 8676, 8826, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 46332, cbuffer, 8826, 9036, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 47142, cbuffer, 5106, 45402, 45582, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 47982, cbuffer, 5166, 45582, 45882, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 48882, cbuffer, 5266, 45882, 46332, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 50592, 46962, 47142, 47982, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 51672, 47682, 47982, 48882, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 53472, 50232, 50592, 51672, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 55272, cbuffer, 9316, 9406, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 55542, cbuffer, 9406, 9556, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 55992, cbuffer, 9556, 9781, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 56667, cbuffer, 9781, 10096, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 57882, cbuffer, 5416, 55272, 55542, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 59142, cbuffer, 5506, 55542, 55992, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 60492, cbuffer, 5656, 55992, 56667, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 63057, 57612, 57882, 59142, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 64677, 58692, 59142, 60492, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 67377, 62517, 63057, 64677, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 70077, cbuffer, 10516, 10642, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 70455, cbuffer, 10642, 10852, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 71085, cbuffer, 10852, 11167, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 72030, cbuffer, 11167, 11608, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 73731, cbuffer, 5881, 70077, 70455, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 75495, cbuffer, 6007, 70455, 71085, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 77385, cbuffer, 6217, 71085, 72030, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 80976, 73353, 73731, 75495, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 83244, 74865, 75495, 77385, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 87024, 80220, 80976, 83244, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 90804, cbuffer, 12196, 12364, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 91308, cbuffer, 12364, 12644, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 92148, cbuffer, 12644, 13064, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 93408, cbuffer, 13064, 13652, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 95676, cbuffer, 6532, 90804, 91308, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 98028, cbuffer, 6700, 91308, 92148, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 100548, cbuffer, 6980, 92148, 93408, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 105336, 95172, 95676, 98028, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 108360, 97188, 98028, 100548, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 113400, 104328, 105336, 108360, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 118440, cbuffer, 14436, 14652, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 119088, cbuffer, 14652, 15012, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 120168, cbuffer, 15012, 15552, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 121788, cbuffer, 15552, 16308, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 124704, cbuffer, 7400, 118440, 119088, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 127728, cbuffer, 7616, 119088, 120168, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 130968, cbuffer, 7976, 120168, 121788, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 137124, 124056, 124704, 127728, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 141012, 126648, 127728, 130968, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 147492, 135828, 137124, 141012, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<3, 2>(skbuffer, 0, ckbuffer, 8070, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 350, ckbuffer, 8670, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 700, ckbuffer, 9270, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 4200, ckbuffer, 21975, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 4725, ckbuffer, 22875, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 5250, ckbuffer, 23775, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 10500, ckbuffer, 41622, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 11235, ckbuffer, 42882, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 11970, ckbuffer, 44142, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 164535, ckbuffer, 53472, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 164885, ckbuffer, 54072, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 165235, ckbuffer, 54672, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 165585, ckbuffer, 67377, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 166110, ckbuffer, 68277, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 166635, ckbuffer, 69177, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 167160, ckbuffer, 87024, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 167895, ckbuffer, 88284, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 168630, ckbuffer, 89544, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 169365, ckbuffer, 113400, 0, 6);

            t4cfunc::ket_transform<3, 2>(skbuffer, 170345, ckbuffer, 115080, 0, 6);

            t4cfunc::ket_transform<3, 2>(skbuffer, 171325, ckbuffer, 116760, 0, 6);

            t4cfunc::ket_transform<3, 2>(skbuffer, 172305, ckbuffer, 147492, 0, 7);

            t4cfunc::ket_transform<3, 2>(skbuffer, 173565, ckbuffer, 149652, 0, 7);

            t4cfunc::ket_transform<3, 2>(skbuffer, 174825, ckbuffer, 151812, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 28140, 0, 4200, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 29190, 350, 4725, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 30240, 700, 5250, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 40740, 4200, 10500, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 42315, 4725, 11235, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 43890, 5250, 11970, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 79485, 28140, 40740, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 81585, 29190, 42315, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 83685, 30240, 43890, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 1050, 164535, 165585, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 5775, 165585, 167160, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 12705, 167160, 169365, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sixx(skbuffer, 19320, 169365, 172305, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 31290, 0, 1050, 5775, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 45465, 4200, 5775, 12705, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_phxx(skbuffer, 59640, 10500, 12705, 19320, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 85785, 28140, 31290, 45465, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dgxx(skbuffer, 104685, 40740, 45465, 59640, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ffxx(skbuffer, 133035, 79485, 85785, 104685, r_ab, 3, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 133035, 3, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1715, skbuffer, 136535, 3, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 3430, skbuffer, 140035, 3, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 5145, skbuffer, 143535, 3, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 6860, skbuffer, 147035, 3, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 8575, skbuffer, 150535, 3, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 10290, skbuffer, 154035, 3, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 12005, skbuffer, 157535, 3, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 13720, skbuffer, 161035, 3, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 3, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFFFD_hpp */
