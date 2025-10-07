#ifndef ElectronRepulsionGeom1010RecFDFF_hpp
#define ElectronRepulsionGeom1010RecFDFF_hpp

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
#include "ElectronRepulsionGeom1010ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FD|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fdff(T& distributor,
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

    CSimdArray<double> pbuffer(24326, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(17316, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(168831, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(161406, 1);

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

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 798, 385, 406, 574, 602, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 834, 406, 427, 602, 630, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 870, 427, 448, 630, 658, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 906, 448, 469, 658, 686, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 942, 469, 490, 686, 714, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 978, 490, 511, 714, 742, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 1014, 511, 532, 742, 770, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1050, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1053, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 1056, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1059, 3, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1068, 4, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1077, 5, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 1086, 6, 29, 32, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1095, 20, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1113, 23, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1131, 26, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1149, 29, 77, 83, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 1167, 32, 83, 89, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1185, 59, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1215, 65, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1245, 71, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1275, 77, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1305, 83, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1335, 89, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1365, 135, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1410, 145, 250, 265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1455, 155, 265, 280, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1500, 165, 280, 295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1545, 175, 295, 310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1590, 185, 310, 325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1635, 250, 385, 406, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1698, 265, 406, 427, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1761, 280, 427, 448, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1824, 295, 448, 469, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1887, 310, 469, 490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1950, 325, 490, 511, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2013, 406, 574, 602, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2097, 427, 602, 630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2181, 448, 630, 658, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2265, 469, 658, 686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2349, 490, 686, 714, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2433, 511, 714, 742, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2517, 602, 798, 834, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2625, 630, 834, 870, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2733, 658, 870, 906, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2841, 686, 906, 942, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2949, 714, 942, 978, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 3057, 742, 978, 1014, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 3165, 3, 4, 1050, 1053, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 3171, 4, 5, 1053, 1056, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3177, 20, 23, 1050, 1059, 1068, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3195, 23, 26, 1053, 1068, 1077, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 3213, 26, 29, 1056, 1077, 1086, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3231, 59, 65, 1059, 1095, 1113, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3267, 65, 71, 1068, 1113, 1131, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3303, 71, 77, 1077, 1131, 1149, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 3339, 77, 83, 1086, 1149, 1167, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3375, 125, 135, 1095, 1185, 1215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3435, 135, 145, 1113, 1215, 1245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3495, 145, 155, 1131, 1245, 1275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3555, 155, 165, 1149, 1275, 1305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 3615, 165, 175, 1167, 1305, 1335, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3675, 235, 250, 1215, 1365, 1410, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3765, 250, 265, 1245, 1410, 1455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3855, 265, 280, 1275, 1455, 1500, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3945, 280, 295, 1305, 1500, 1545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 4035, 295, 310, 1335, 1545, 1590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4125, 385, 406, 1410, 1635, 1698, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4251, 406, 427, 1455, 1698, 1761, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4377, 427, 448, 1500, 1761, 1824, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4503, 448, 469, 1545, 1824, 1887, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 4629, 469, 490, 1590, 1887, 1950, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4755, 574, 602, 1698, 2013, 2097, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4923, 602, 630, 1761, 2097, 2181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5091, 630, 658, 1824, 2181, 2265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5259, 658, 686, 1887, 2265, 2349, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 5427, 686, 714, 1950, 2349, 2433, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 5595, 798, 834, 2097, 2517, 2625, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 5811, 834, 870, 2181, 2625, 2733, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6027, 870, 906, 2265, 2733, 2841, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6243, 906, 942, 2349, 2841, 2949, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 6459, 942, 978, 2433, 2949, 3057, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 6675, 1050, 1053, 3165, 3171, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 6685, 1059, 1068, 3165, 3177, 3195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 6715, 1068, 1077, 3171, 3195, 3213, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 6745, 1095, 1113, 3177, 3231, 3267, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 6805, 1113, 1131, 3195, 3267, 3303, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 6865, 1131, 1149, 3213, 3303, 3339, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 6925, 1185, 1215, 3231, 3375, 3435, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 7025, 1215, 1245, 3267, 3435, 3495, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 7125, 1245, 1275, 3303, 3495, 3555, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 7225, 1275, 1305, 3339, 3555, 3615, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 7325, 1365, 1410, 3435, 3675, 3765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 7475, 1410, 1455, 3495, 3765, 3855, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 7625, 1455, 1500, 3555, 3855, 3945, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 7775, 1500, 1545, 3615, 3945, 4035, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 7925, 1635, 1698, 3765, 4125, 4251, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 8135, 1698, 1761, 3855, 4251, 4377, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 8345, 1761, 1824, 3945, 4377, 4503, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 8555, 1824, 1887, 4035, 4503, 4629, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 8765, 2013, 2097, 4251, 4755, 4923, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 9045, 2097, 2181, 4377, 4923, 5091, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 9325, 2181, 2265, 4503, 5091, 5259, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 9605, 2265, 2349, 4629, 5259, 5427, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 9885, 2517, 2625, 4923, 5595, 5811, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 10245, 2625, 2733, 5091, 5811, 6027, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 10605, 2733, 2841, 5259, 6027, 6243, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 10965, 2841, 2949, 5427, 6243, 6459, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 11325, 3177, 3195, 6675, 6685, 6715, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 11370, 3231, 3267, 6685, 6745, 6805, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 11460, 3267, 3303, 6715, 6805, 6865, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 11550, 3375, 3435, 6745, 6925, 7025, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 11700, 3435, 3495, 6805, 7025, 7125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 11850, 3495, 3555, 6865, 7125, 7225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 12000, 3675, 3765, 7025, 7325, 7475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 12225, 3765, 3855, 7125, 7475, 7625, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 12450, 3855, 3945, 7225, 7625, 7775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 12675, 4125, 4251, 7475, 7925, 8135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 12990, 4251, 4377, 7625, 8135, 8345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 13305, 4377, 4503, 7775, 8345, 8555, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 13620, 4755, 4923, 8135, 8765, 9045, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 14040, 4923, 5091, 8345, 9045, 9325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 14460, 5091, 5259, 8555, 9325, 9605, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 14880, 5595, 5811, 9045, 9885, 10245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 15420, 5811, 6027, 9325, 10245, 10605, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 15960, 6027, 6243, 9605, 10605, 10965, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 16500, 6745, 6805, 11325, 11370, 11460, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 16626, 6925, 7025, 11370, 11550, 11700, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 16836, 7025, 7125, 11460, 11700, 11850, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 17046, 7325, 7475, 11700, 12000, 12225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 17361, 7475, 7625, 11850, 12225, 12450, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 17676, 7925, 8135, 12225, 12675, 12990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 18117, 8135, 8345, 12450, 12990, 13305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 18558, 8765, 9045, 12990, 13620, 14040, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 19146, 9045, 9325, 13305, 14040, 14460, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsk(pbuffer, 19734, 9885, 10245, 14040, 14880, 15420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsk(pbuffer, 20490, 10245, 10605, 14460, 15420, 15960, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 21246, 11550, 11700, 16500, 16626, 16836, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 21526, 12000, 12225, 16836, 17046, 17361, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 21946, 12675, 12990, 17361, 17676, 18117, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 22534, 13620, 14040, 18117, 18558, 19146, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisk(pbuffer, 23318, 14880, 15420, 19146, 19734, 20490, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 3375, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 3675, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 150, pbuffer, 4125, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 276, pbuffer, 6925, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 376, pbuffer, 7325, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 526, pbuffer, 7925, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 736, pbuffer, 11550, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 886, pbuffer, 12000, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1111, pbuffer, 12675, 315, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {3375, 3435});

                pbuffer.scale(2.0 * a_exp, {3675, 3765});

                pbuffer.scale(2.0 * a_exp, {4125, 4251});

                pbuffer.scale(2.0 * a_exp, {6925, 7025});

                pbuffer.scale(2.0 * a_exp, {7325, 7475});

                pbuffer.scale(2.0 * a_exp, {7925, 8135});

                pbuffer.scale(2.0 * a_exp, {11550, 11700});

                pbuffer.scale(2.0 * a_exp, {12000, 12225});

                pbuffer.scale(2.0 * a_exp, {12675, 12990});

                pbuffer.scale(2.0 * a_exp, {16626, 16836});

                pbuffer.scale(2.0 * a_exp, {17046, 17361});

                pbuffer.scale(2.0 * a_exp, {17676, 18117});

                pbuffer.scale(2.0 * a_exp, {21246, 21526});

                pbuffer.scale(2.0 * a_exp, {21526, 21946});

                pbuffer.scale(2.0 * a_exp, {21946, 22534});

                t2cfunc::reduce(cbuffer, 4836, pbuffer, 3375, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4896, pbuffer, 3675, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4986, pbuffer, 4125, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5112, pbuffer, 6925, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5212, pbuffer, 7325, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5362, pbuffer, 7925, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5572, pbuffer, 11550, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5722, pbuffer, 12000, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5947, pbuffer, 12675, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6262, pbuffer, 16626, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6472, pbuffer, 17046, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6787, pbuffer, 17676, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7228, pbuffer, 21246, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7508, pbuffer, 21526, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7928, pbuffer, 21946, 588, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3375, 3435});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3675, 3765});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {4125, 4251});

                pbuffer.scale(pfactors, 0, 2.0, {4755, 4923});

                pbuffer.scale(pfactors, 0, 2.0, {5595, 5811});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {6925, 7025});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {7325, 7475});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {7925, 8135});

                pbuffer.scale(pfactors, 0, 2.0, {8765, 9045});

                pbuffer.scale(pfactors, 0, 2.0, {9885, 10245});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {11550, 11700});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {12000, 12225});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {12675, 12990});

                pbuffer.scale(pfactors, 0, 2.0, {13620, 14040});

                pbuffer.scale(pfactors, 0, 2.0, {14880, 15420});

                t2cfunc::reduce(cbuffer, 1426, pbuffer, 3375, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1486, pbuffer, 3675, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1576, pbuffer, 4125, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1702, pbuffer, 4755, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1870, pbuffer, 5595, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2086, pbuffer, 6925, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2186, pbuffer, 7325, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2336, pbuffer, 7925, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2546, pbuffer, 8765, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2826, pbuffer, 9885, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3186, pbuffer, 11550, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3336, pbuffer, 12000, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3561, pbuffer, 12675, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3876, pbuffer, 13620, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4296, pbuffer, 14880, 540, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {3375, 3435});

                pbuffer.scale(2.0 * a_exp, {3675, 3765});

                pbuffer.scale(2.0 * a_exp, {4125, 4251});

                pbuffer.scale(2.0 * a_exp, {4755, 4923});

                pbuffer.scale(2.0 * a_exp, {5595, 5811});

                pbuffer.scale(2.0 * a_exp, {6925, 7025});

                pbuffer.scale(2.0 * a_exp, {7325, 7475});

                pbuffer.scale(2.0 * a_exp, {7925, 8135});

                pbuffer.scale(2.0 * a_exp, {8765, 9045});

                pbuffer.scale(2.0 * a_exp, {9885, 10245});

                pbuffer.scale(2.0 * a_exp, {11550, 11700});

                pbuffer.scale(2.0 * a_exp, {12000, 12225});

                pbuffer.scale(2.0 * a_exp, {12675, 12990});

                pbuffer.scale(2.0 * a_exp, {13620, 14040});

                pbuffer.scale(2.0 * a_exp, {14880, 15420});

                pbuffer.scale(pfactors, 0, 2.0, {16626, 16836});

                pbuffer.scale(pfactors, 0, 2.0, {17046, 17361});

                pbuffer.scale(pfactors, 0, 2.0, {17676, 18117});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {18558, 19146});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {19734, 20490});

                pbuffer.scale(pfactors, 0, 2.0, {21246, 21526});

                pbuffer.scale(pfactors, 0, 2.0, {21526, 21946});

                pbuffer.scale(pfactors, 0, 2.0, {21946, 22534});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {22534, 23318});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {23318, 24326});

                t2cfunc::reduce(cbuffer, 8516, pbuffer, 3375, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8576, pbuffer, 3675, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8666, pbuffer, 4125, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8792, pbuffer, 4755, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8960, pbuffer, 5595, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9176, pbuffer, 6925, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9276, pbuffer, 7325, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9426, pbuffer, 7925, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9636, pbuffer, 8765, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9916, pbuffer, 9885, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10276, pbuffer, 11550, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10426, pbuffer, 12000, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10651, pbuffer, 12675, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10966, pbuffer, 13620, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11386, pbuffer, 14880, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11926, pbuffer, 16626, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12136, pbuffer, 17046, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12451, pbuffer, 17676, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12892, pbuffer, 18558, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13480, pbuffer, 19734, 756, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14236, pbuffer, 21246, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14516, pbuffer, 21526, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14936, pbuffer, 21946, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 15524, pbuffer, 22534, 784, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16308, pbuffer, 23318, 1008, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1332, cbuffer, 0, 60, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 2052, cbuffer, 60, 150, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 4266, 1332, 2052, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 11346, cbuffer, 276, 376, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 12546, cbuffer, 376, 526, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 16236, 11346, 12546, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 27666, cbuffer, 736, 886, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 29466, cbuffer, 886, 1111, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 35001, 27666, 29466, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 48483, cbuffer, 4836, 4896, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 49203, cbuffer, 4896, 4986, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 51417, 48483, 49203, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 58497, cbuffer, 5112, 5212, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 59697, cbuffer, 5212, 5362, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 63387, 58497, 59697, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 74817, cbuffer, 5572, 5722, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 76617, cbuffer, 5722, 5947, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 82152, 74817, 76617, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 98964, cbuffer, 6262, 6472, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 101484, cbuffer, 6472, 6787, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 109233, 98964, 101484, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 132459, cbuffer, 7228, 7508, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 135819, cbuffer, 7508, 7928, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 146151, 132459, 135819, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 0, cbuffer, 1426, 1486, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 180, cbuffer, 1486, 1576, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 450, cbuffer, 1576, 1702, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 828, cbuffer, 1702, 1870, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 1512, cbuffer, 0, 0, 180, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 2322, cbuffer, 60, 180, 450, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 3132, cbuffer, 150, 450, 828, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 4626, 1332, 1512, 2322, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 5706, 2052, 2322, 3132, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 7326, 4266, 4626, 5706, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 9126, cbuffer, 2086, 2186, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 9426, cbuffer, 2186, 2336, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 9876, cbuffer, 2336, 2546, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 10506, cbuffer, 2546, 2826, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 11646, cbuffer, 276, 9126, 9426, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 12996, cbuffer, 376, 9426, 9876, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 14346, cbuffer, 526, 9876, 10506, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 16836, 11346, 11646, 12996, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 18636, 12546, 12996, 14346, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 21336, 16236, 16836, 18636, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 24336, cbuffer, 3186, 3336, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 24786, cbuffer, 3336, 3561, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 25461, cbuffer, 3561, 3876, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 26406, cbuffer, 3876, 4296, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 28116, cbuffer, 736, 24336, 24786, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 30141, cbuffer, 886, 24786, 25461, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 32166, cbuffer, 1111, 25461, 26406, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 35901, 27666, 28116, 30141, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 38601, 29466, 30141, 32166, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 42651, 35001, 35901, 38601, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 47151, cbuffer, 8516, 8576, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 47331, cbuffer, 8576, 8666, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 47601, cbuffer, 8666, 8792, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 47979, cbuffer, 8792, 8960, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 48663, cbuffer, 4836, 47151, 47331, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 49473, cbuffer, 4896, 47331, 47601, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 50283, cbuffer, 4986, 47601, 47979, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 51777, 48483, 48663, 49473, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 52857, 49203, 49473, 50283, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 54477, 51417, 51777, 52857, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 56277, cbuffer, 9176, 9276, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 56577, cbuffer, 9276, 9426, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 57027, cbuffer, 9426, 9636, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 57657, cbuffer, 9636, 9916, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 58797, cbuffer, 5112, 56277, 56577, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 60147, cbuffer, 5212, 56577, 57027, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 61497, cbuffer, 5362, 57027, 57657, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 63987, 58497, 58797, 60147, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 65787, 59697, 60147, 61497, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 68487, 63387, 63987, 65787, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 71487, cbuffer, 10276, 10426, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 71937, cbuffer, 10426, 10651, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 72612, cbuffer, 10651, 10966, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 73557, cbuffer, 10966, 11386, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 75267, cbuffer, 5572, 71487, 71937, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 77292, cbuffer, 5722, 71937, 72612, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 79317, cbuffer, 5947, 72612, 73557, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 83052, 74817, 75267, 77292, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 85752, 76617, 77292, 79317, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 89802, 82152, 83052, 85752, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 94302, cbuffer, 11926, 12136, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 94932, cbuffer, 12136, 12451, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 95877, cbuffer, 12451, 12892, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 97200, cbuffer, 12892, 13480, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 99594, cbuffer, 6262, 94302, 94932, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 102429, cbuffer, 6472, 94932, 95877, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 105264, cbuffer, 6787, 95877, 97200, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 110493, 98964, 99594, 102429, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 114273, 101484, 102429, 105264, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 119943, 109233, 110493, 114273, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 126243, cbuffer, 14236, 14516, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 127083, cbuffer, 14516, 14936, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 128343, cbuffer, 14936, 15524, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 130107, cbuffer, 15524, 16308, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 133299, cbuffer, 7228, 126243, 127083, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 137079, cbuffer, 7508, 127083, 128343, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 140859, cbuffer, 7928, 128343, 130107, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 147831, 132459, 133299, 137079, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 152871, 135819, 137079, 140859, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 160431, 146151, 147831, 152871, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 7326, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 294, ckbuffer, 7926, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 588, ckbuffer, 8526, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 3528, ckbuffer, 21336, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 4018, ckbuffer, 22336, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 4508, ckbuffer, 23336, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 9408, ckbuffer, 42651, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 10143, ckbuffer, 44151, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 10878, ckbuffer, 45651, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 149646, ckbuffer, 54477, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 149940, ckbuffer, 55077, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 150234, ckbuffer, 55677, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 150528, ckbuffer, 68487, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 151018, ckbuffer, 69487, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 151508, ckbuffer, 70487, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 151998, ckbuffer, 89802, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 152733, ckbuffer, 91302, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 153468, ckbuffer, 92802, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 154203, ckbuffer, 119943, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 155232, ckbuffer, 122043, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 156261, ckbuffer, 124143, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 157290, ckbuffer, 160431, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 158662, ckbuffer, 163231, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 160034, ckbuffer, 166031, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 27489, 0, 3528, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 28371, 294, 4018, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 29253, 588, 4508, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 38073, 3528, 9408, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 39543, 4018, 10143, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 41013, 4508, 10878, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 75558, 27489, 38073, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 77322, 28371, 39543, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 79086, 29253, 41013, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 882, 149646, 150528, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 4998, 150528, 151998, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 11613, 151998, 154203, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 18228, 154203, 157290, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 30135, 0, 882, 4998, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 42483, 3528, 4998, 11613, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 55713, 9408, 11613, 18228, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 80850, 27489, 30135, 42483, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 96726, 38073, 42483, 55713, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fdxx(skbuffer, 123186, 75558, 80850, 96726, r_ab, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 123186, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1715, skbuffer, 126126, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 3430, skbuffer, 129066, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 5145, skbuffer, 132006, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 6860, skbuffer, 134946, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 8575, skbuffer, 137886, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 10290, skbuffer, 140826, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 12005, skbuffer, 143766, 3, 3);

            t4cfunc::bra_transform<3, 2>(sbuffer, 13720, skbuffer, 146706, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFDFF_hpp */
