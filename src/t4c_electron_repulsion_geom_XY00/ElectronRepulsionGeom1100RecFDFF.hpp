#ifndef ElectronRepulsionGeom1100RecFDFF_hpp
#define ElectronRepulsionGeom1100RecFDFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecXXDF.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXFF.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
#include "ElectronRepulsionGeom0100ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSIXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecFDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSHXX.hpp"
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
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
#include "ElectronRepulsionPrimRecSISI.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FD|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_fdff(T& distributor,
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

    CSimdArray<double> pbuffer(26269, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(17020, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(122540, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(190463, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 798, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 801, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 804, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 807, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 810, 3, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 819, 4, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 828, 5, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 837, 6, 29, 32, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 846, 7, 32, 35, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 855, 20, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 873, 23, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 891, 26, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 909, 29, 77, 83, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 927, 32, 83, 89, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 945, 35, 89, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 963, 59, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 993, 65, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1023, 71, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1053, 77, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1083, 83, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1113, 89, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1143, 95, 185, 195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1173, 135, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1218, 145, 250, 265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1263, 155, 265, 280, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1308, 165, 280, 295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1353, 175, 295, 310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1398, 185, 310, 325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1443, 195, 325, 340, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1488, 250, 385, 406, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1551, 265, 406, 427, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1614, 280, 427, 448, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1677, 295, 448, 469, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1740, 310, 469, 490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1803, 325, 490, 511, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1866, 340, 511, 532, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1929, 406, 574, 602, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2013, 427, 602, 630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2097, 448, 630, 658, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2181, 469, 658, 686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2265, 490, 686, 714, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2349, 511, 714, 742, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2433, 532, 742, 770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2517, 3, 4, 798, 801, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2523, 4, 5, 801, 804, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2529, 5, 6, 804, 807, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2535, 20, 23, 798, 810, 819, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2553, 23, 26, 801, 819, 828, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2571, 26, 29, 804, 828, 837, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2589, 29, 32, 807, 837, 846, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2607, 59, 65, 810, 855, 873, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2643, 65, 71, 819, 873, 891, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2679, 71, 77, 828, 891, 909, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2715, 77, 83, 837, 909, 927, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2751, 83, 89, 846, 927, 945, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2787, 125, 135, 855, 963, 993, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2847, 135, 145, 873, 993, 1023, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2907, 145, 155, 891, 1023, 1053, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2967, 155, 165, 909, 1053, 1083, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3027, 165, 175, 927, 1083, 1113, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3087, 175, 185, 945, 1113, 1143, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3147, 235, 250, 993, 1173, 1218, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3237, 250, 265, 1023, 1218, 1263, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3327, 265, 280, 1053, 1263, 1308, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3417, 280, 295, 1083, 1308, 1353, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3507, 295, 310, 1113, 1353, 1398, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3597, 310, 325, 1143, 1398, 1443, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3687, 385, 406, 1218, 1488, 1551, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3813, 406, 427, 1263, 1551, 1614, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3939, 427, 448, 1308, 1614, 1677, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4065, 448, 469, 1353, 1677, 1740, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4191, 469, 490, 1398, 1740, 1803, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4317, 490, 511, 1443, 1803, 1866, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4443, 574, 602, 1551, 1929, 2013, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4611, 602, 630, 1614, 2013, 2097, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4779, 630, 658, 1677, 2097, 2181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4947, 658, 686, 1740, 2181, 2265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5115, 686, 714, 1803, 2265, 2349, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5283, 714, 742, 1866, 2349, 2433, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 5451, 798, 801, 2517, 2523, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 5461, 801, 804, 2523, 2529, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 5471, 810, 819, 2517, 2535, 2553, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 5501, 819, 828, 2523, 2553, 2571, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 5531, 828, 837, 2529, 2571, 2589, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 5561, 855, 873, 2535, 2607, 2643, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 5621, 873, 891, 2553, 2643, 2679, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 5681, 891, 909, 2571, 2679, 2715, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 5741, 909, 927, 2589, 2715, 2751, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 5801, 963, 993, 2607, 2787, 2847, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 5901, 993, 1023, 2643, 2847, 2907, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 6001, 1023, 1053, 2679, 2907, 2967, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 6101, 1053, 1083, 2715, 2967, 3027, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 6201, 1083, 1113, 2751, 3027, 3087, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 6301, 1173, 1218, 2847, 3147, 3237, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 6451, 1218, 1263, 2907, 3237, 3327, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 6601, 1263, 1308, 2967, 3327, 3417, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 6751, 1308, 1353, 3027, 3417, 3507, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 6901, 1353, 1398, 3087, 3507, 3597, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 7051, 1488, 1551, 3237, 3687, 3813, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 7261, 1551, 1614, 3327, 3813, 3939, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 7471, 1614, 1677, 3417, 3939, 4065, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 7681, 1677, 1740, 3507, 4065, 4191, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 7891, 1740, 1803, 3597, 4191, 4317, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 8101, 1929, 2013, 3813, 4443, 4611, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 8381, 2013, 2097, 3939, 4611, 4779, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 8661, 2097, 2181, 4065, 4779, 4947, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 8941, 2181, 2265, 4191, 4947, 5115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 9221, 2265, 2349, 4317, 5115, 5283, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 9501, 2517, 2523, 5451, 5461, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 9516, 2535, 2553, 5451, 5471, 5501, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 9561, 2553, 2571, 5461, 5501, 5531, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 9606, 2607, 2643, 5471, 5561, 5621, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 9696, 2643, 2679, 5501, 5621, 5681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 9786, 2679, 2715, 5531, 5681, 5741, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 9876, 2787, 2847, 5561, 5801, 5901, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 10026, 2847, 2907, 5621, 5901, 6001, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 10176, 2907, 2967, 5681, 6001, 6101, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 10326, 2967, 3027, 5741, 6101, 6201, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 10476, 3147, 3237, 5901, 6301, 6451, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 10701, 3237, 3327, 6001, 6451, 6601, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 10926, 3327, 3417, 6101, 6601, 6751, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 11151, 3417, 3507, 6201, 6751, 6901, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 11376, 3687, 3813, 6451, 7051, 7261, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 11691, 3813, 3939, 6601, 7261, 7471, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 12006, 3939, 4065, 6751, 7471, 7681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 12321, 4065, 4191, 6901, 7681, 7891, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 12636, 4443, 4611, 7261, 8101, 8381, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 13056, 4611, 4779, 7471, 8381, 8661, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 13476, 4779, 4947, 7681, 8661, 8941, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 13896, 4947, 5115, 7891, 8941, 9221, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 14316, 5471, 5501, 9501, 9516, 9561, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 14379, 5561, 5621, 9516, 9606, 9696, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 14505, 5621, 5681, 9561, 9696, 9786, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 14631, 5801, 5901, 9606, 9876, 10026, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 14841, 5901, 6001, 9696, 10026, 10176, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 15051, 6001, 6101, 9786, 10176, 10326, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 15261, 6301, 6451, 10026, 10476, 10701, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 15576, 6451, 6601, 10176, 10701, 10926, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 15891, 6601, 6751, 10326, 10926, 11151, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 16206, 7051, 7261, 10701, 11376, 11691, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 16647, 7261, 7471, 10926, 11691, 12006, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 17088, 7471, 7681, 11151, 12006, 12321, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 17529, 8101, 8381, 11691, 12636, 13056, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 18117, 8381, 8661, 12006, 13056, 13476, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 18705, 8661, 8941, 12321, 13476, 13896, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 19293, 9606, 9696, 14316, 14379, 14505, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 19461, 9876, 10026, 14379, 14631, 14841, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 19741, 10026, 10176, 14505, 14841, 15051, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 20021, 10476, 10701, 14841, 15261, 15576, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 20441, 10701, 10926, 15051, 15576, 15891, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 20861, 11376, 11691, 15576, 16206, 16647, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 21449, 11691, 12006, 15891, 16647, 17088, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 22037, 12636, 13056, 16647, 17529, 18117, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 22821, 13056, 13476, 17088, 18117, 18705, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 23605, 14631, 14841, 19293, 19461, 19741, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 23965, 15261, 15576, 19741, 20021, 20441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksh(pbuffer, 24505, 16206, 16647, 20441, 20861, 21449, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksi(pbuffer, 25261, 17529, 18117, 21449, 22037, 22821, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 963, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 1173, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 75, pbuffer, 1488, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 138, pbuffer, 1929, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 222, pbuffer, 2787, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 282, pbuffer, 3147, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 372, pbuffer, 3687, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 498, pbuffer, 4443, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 666, pbuffer, 5801, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 766, pbuffer, 6301, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 916, pbuffer, 7051, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1126, pbuffer, 8101, 280, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {5801, 5901});

                pbuffer.scale(2.0 * b_exp, {6301, 6451});

                pbuffer.scale(2.0 * b_exp, {7051, 7261});

                pbuffer.scale(2.0 * b_exp, {8101, 8381});

                pbuffer.scale(2.0 * b_exp, {9876, 10026});

                pbuffer.scale(2.0 * b_exp, {10476, 10701});

                pbuffer.scale(2.0 * b_exp, {11376, 11691});

                pbuffer.scale(2.0 * b_exp, {12636, 13056});

                pbuffer.scale(2.0 * b_exp, {14631, 14841});

                pbuffer.scale(2.0 * b_exp, {15261, 15576});

                pbuffer.scale(2.0 * b_exp, {16206, 16647});

                pbuffer.scale(2.0 * b_exp, {17529, 18117});

                t2cfunc::reduce(cbuffer, 1406, pbuffer, 5801, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1506, pbuffer, 6301, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1656, pbuffer, 7051, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1866, pbuffer, 8101, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2146, pbuffer, 9876, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2296, pbuffer, 10476, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2521, pbuffer, 11376, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2836, pbuffer, 12636, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3256, pbuffer, 14631, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3466, pbuffer, 15261, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3781, pbuffer, 16206, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4222, pbuffer, 17529, 588, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {963, 993});

                pbuffer.scale(2.0 * a_exp, {1173, 1218});

                pbuffer.scale(2.0 * a_exp, {1488, 1551});

                pbuffer.scale(2.0 * a_exp, {1929, 2013});

                pbuffer.scale(2.0 * a_exp, {2787, 2847});

                pbuffer.scale(2.0 * a_exp, {3147, 3237});

                pbuffer.scale(2.0 * a_exp, {3687, 3813});

                pbuffer.scale(2.0 * a_exp, {4443, 4611});

                pbuffer.scale(a_exp / b_exp, {5801, 5901});

                pbuffer.scale(a_exp / b_exp, {6301, 6451});

                pbuffer.scale(a_exp / b_exp, {7051, 7261});

                pbuffer.scale(a_exp / b_exp, {8101, 8381});

                pbuffer.scale(a_exp / b_exp, {9876, 10026});

                pbuffer.scale(a_exp / b_exp, {10476, 10701});

                pbuffer.scale(a_exp / b_exp, {11376, 11691});

                pbuffer.scale(a_exp / b_exp, {12636, 13056});

                pbuffer.scale(a_exp / b_exp, {14631, 14841});

                pbuffer.scale(a_exp / b_exp, {15261, 15576});

                pbuffer.scale(a_exp / b_exp, {16206, 16647});

                pbuffer.scale(a_exp / b_exp, {17529, 18117});

                t2cfunc::reduce(cbuffer, 4810, pbuffer, 963, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4840, pbuffer, 1173, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4885, pbuffer, 1488, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4948, pbuffer, 1929, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5032, pbuffer, 2787, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5092, pbuffer, 3147, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5182, pbuffer, 3687, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5308, pbuffer, 4443, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5476, pbuffer, 5801, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5576, pbuffer, 6301, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5726, pbuffer, 7051, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5936, pbuffer, 8101, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6216, pbuffer, 9876, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6366, pbuffer, 10476, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6591, pbuffer, 11376, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6906, pbuffer, 12636, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7326, pbuffer, 14631, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7536, pbuffer, 15261, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7851, pbuffer, 16206, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8292, pbuffer, 17529, 588, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {5801, 5901});

                pbuffer.scale(2.0 * b_exp, {6301, 6451});

                pbuffer.scale(2.0 * b_exp, {7051, 7261});

                pbuffer.scale(2.0 * b_exp, {8101, 8381});

                pbuffer.scale(2.0 * b_exp, {9876, 10026});

                pbuffer.scale(2.0 * b_exp, {10476, 10701});

                pbuffer.scale(2.0 * b_exp, {11376, 11691});

                pbuffer.scale(2.0 * b_exp, {12636, 13056});

                pbuffer.scale(2.0 * b_exp, {14631, 14841});

                pbuffer.scale(2.0 * b_exp, {15261, 15576});

                pbuffer.scale(2.0 * b_exp, {16206, 16647});

                pbuffer.scale(2.0 * b_exp, {17529, 18117});

                pbuffer.scale(4.0 * a_exp * b_exp, {19461, 19741});

                pbuffer.scale(4.0 * a_exp * b_exp, {20021, 20441});

                pbuffer.scale(4.0 * a_exp * b_exp, {20861, 21449});

                pbuffer.scale(4.0 * a_exp * b_exp, {22037, 22821});

                pbuffer.scale(4.0 * a_exp * b_exp, {23605, 23965});

                pbuffer.scale(4.0 * a_exp * b_exp, {23965, 24505});

                pbuffer.scale(4.0 * a_exp * b_exp, {24505, 25261});

                pbuffer.scale(4.0 * a_exp * b_exp, {25261, 26269});

                t2cfunc::reduce(cbuffer, 8880, pbuffer, 5801, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8980, pbuffer, 6301, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9130, pbuffer, 7051, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9340, pbuffer, 8101, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9620, pbuffer, 9876, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9770, pbuffer, 10476, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9995, pbuffer, 11376, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10310, pbuffer, 12636, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10730, pbuffer, 14631, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10940, pbuffer, 15261, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11255, pbuffer, 16206, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11696, pbuffer, 17529, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12284, pbuffer, 19461, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12564, pbuffer, 20021, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12984, pbuffer, 20861, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13572, pbuffer, 22037, 784, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14356, pbuffer, 23605, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14716, pbuffer, 23965, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15256, pbuffer, 24505, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16012, pbuffer, 25261, 1008, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 30, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 90, cbuffer, 30, 75, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 225, cbuffer, 75, 138, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 414, 0, 90, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 594, 90, 225, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 864, 414, 594, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1164, cbuffer, 222, 282, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 1344, cbuffer, 282, 372, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 1614, cbuffer, 372, 498, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 1992, 1164, 1344, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 2352, 1344, 1614, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 2892, 1992, 2352, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5292, cbuffer, 666, 766, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 5592, cbuffer, 766, 916, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 6042, cbuffer, 916, 1126, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 6672, 5292, 5592, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 7272, 5592, 6042, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 8172, 6672, 7272, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 16672, cbuffer, 1406, 1506, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 16972, cbuffer, 1506, 1656, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 17422, cbuffer, 1656, 1866, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 18052, 16672, 16972, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 18652, 16972, 17422, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 19552, 18052, 18652, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 20552, cbuffer, 2146, 2296, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 21002, cbuffer, 2296, 2521, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 21677, cbuffer, 2521, 2836, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 22622, 20552, 21002, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 23522, 21002, 21677, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 24872, 22622, 23522, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 26372, cbuffer, 3256, 3466, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 27002, cbuffer, 3466, 3781, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 27947, cbuffer, 3781, 4222, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 29270, 26372, 27002, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 30530, 27002, 27947, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 32420, 29270, 30530, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 34520, cbuffer, 4810, 4840, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 34610, cbuffer, 4840, 4885, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 34745, cbuffer, 4885, 4948, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 34934, 34520, 34610, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 35114, 34610, 34745, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 35384, 34934, 35114, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 35684, cbuffer, 5032, 5092, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 35864, cbuffer, 5092, 5182, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 36134, cbuffer, 5182, 5308, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 36512, 35684, 35864, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 36872, 35864, 36134, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 37412, 36512, 36872, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 39812, cbuffer, 5476, 5576, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 40112, cbuffer, 5576, 5726, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 40562, cbuffer, 5726, 5936, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 41192, 39812, 40112, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 41792, 40112, 40562, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 42692, 41192, 41792, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 46692, cbuffer, 6216, 6366, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 47142, cbuffer, 6366, 6591, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 47817, cbuffer, 6591, 6906, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 48762, 46692, 47142, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 49662, 47142, 47817, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 51012, 48762, 49662, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 57012, cbuffer, 7326, 7536, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 57642, cbuffer, 7536, 7851, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 58587, cbuffer, 7851, 8292, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 59910, 57012, 57642, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 61170, 57642, 58587, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 63060, 59910, 61170, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 79860, cbuffer, 8880, 8980, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 80160, cbuffer, 8980, 9130, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 80610, cbuffer, 9130, 9340, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 81240, 79860, 80160, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 81840, 80160, 80610, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 82740, 81240, 81840, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 83740, cbuffer, 9620, 9770, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 84190, cbuffer, 9770, 9995, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 84865, cbuffer, 9995, 10310, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 85810, 83740, 84190, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 86710, 84190, 84865, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 88060, 85810, 86710, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 89560, cbuffer, 10730, 10940, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 90190, cbuffer, 10940, 11255, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 91135, cbuffer, 11255, 11696, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 92458, 89560, 90190, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 93718, 90190, 91135, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 95608, 92458, 93718, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 97708, cbuffer, 12284, 12564, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 98548, cbuffer, 12564, 12984, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 99808, cbuffer, 12984, 13572, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 101572, 97708, 98548, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 103252, 98548, 99808, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 105772, 101572, 103252, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 108572, cbuffer, 14356, 14716, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 109652, cbuffer, 14716, 15256, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 111272, cbuffer, 15256, 16012, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 113540, 108572, 109652, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 115700, 109652, 111272, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 118940, 113540, 115700, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 864, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 147, ckbuffer, 2892, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 3969, ckbuffer, 8172, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 163807, ckbuffer, 19552, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 164297, ckbuffer, 24872, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 165032, ckbuffer, 32420, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 166061, ckbuffer, 35384, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 166208, ckbuffer, 37412, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 167384, ckbuffer, 42692, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 169344, ckbuffer, 51012, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 172284, ckbuffer, 63060, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 185073, ckbuffer, 82740, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 185563, ckbuffer, 88060, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 186298, ckbuffer, 95608, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 187327, ckbuffer, 105772, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 188699, ckbuffer, 118940, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 28420, 147, 3969, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 180516, 166208, 167384, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 181398, 167384, 169344, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 182868, 169344, 172284, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 31948, 147, 180516, 181398, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 46942, 3969, 181398, 182868, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ddxx(skbuffer, 89719, 28420, 31948, 46942, r_ab, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 441, 0, 163807, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 4459, 147, 164297, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 10339, 3969, 165032, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 29302, 147, 441, 4459, r_ab, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 42532, 3969, 4459, 10339, r_ab, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ddxx(skbuffer, 84427, 28420, 29302, 42532, r_ab, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 166502, 166061, 185073, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 167874, 166208, 185563, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 170079, 167384, 186298, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 173313, 169344, 187327, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 176400, 172284, 188699, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 1323, 166208, 166502, 167874, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 5929, 167384, 167874, 170079, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 12544, 169344, 170079, 173313, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 19159, 172284, 173313, 176400, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 34594, 441, 180516, 1323, 5929, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 51352, 4459, 181398, 5929, 12544, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 64582, 10339, 182868, 12544, 19159, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 95011, 29302, 31948, 34594, 51352, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 110887, 42532, 46942, 51352, 64582, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_fdxx(skbuffer, 137347, 84427, 89719, 95011, 110887, r_ab, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 137347, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1715, skbuffer, 140287, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 3430, skbuffer, 143227, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 5145, skbuffer, 146167, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 6860, skbuffer, 149107, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 8575, skbuffer, 152047, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 10290, skbuffer, 154987, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 12005, skbuffer, 157927, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 13720, skbuffer, 160867, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFDFF_hpp */
