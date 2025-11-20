#ifndef ElectronRepulsionGeom1010RecDFFD_hpp
#define ElectronRepulsionGeom1010RecDFFD_hpp

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
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DF|1/|r-r'||FD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dffd(T& distributor,
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

    CSimdArray<double> pbuffer(17379, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(10989, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(97713, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(70560, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(11025, 1);

    // setup Boys fuction data

    const CBoysFunc<12> bf_table;

    CSimdArray<double> bf_data(14, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 13, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 13, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 13, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 13);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 25, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 28, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 31, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 34, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 37, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 40, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 43, 10, 11, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 46, 11, 12, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 49, 0, 1, 13, 16, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 55, 1, 2, 16, 19, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 61, 2, 3, 19, 22, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 67, 3, 4, 22, 25, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 73, 4, 5, 25, 28, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 79, 5, 6, 28, 31, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 85, 6, 7, 31, 34, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 91, 7, 8, 34, 37, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 97, 8, 9, 37, 40, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 103, 9, 10, 40, 43, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 109, 10, 11, 43, 46, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 13, 16, 49, 55, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 16, 19, 55, 61, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 19, 22, 61, 67, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 22, 25, 67, 73, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 25, 28, 73, 79, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 28, 31, 79, 85, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 175, 31, 34, 85, 91, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 185, 34, 37, 91, 97, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 195, 37, 40, 97, 103, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 205, 40, 43, 103, 109, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 215, 49, 55, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 230, 55, 61, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 245, 61, 67, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 260, 67, 73, 145, 155, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 275, 73, 79, 155, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 290, 79, 85, 165, 175, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 305, 85, 91, 175, 185, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 320, 91, 97, 185, 195, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 335, 97, 103, 195, 205, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 350, 115, 125, 215, 230, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 371, 125, 135, 230, 245, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 392, 135, 145, 245, 260, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 413, 145, 155, 260, 275, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 434, 155, 165, 275, 290, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 455, 165, 175, 290, 305, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 476, 175, 185, 305, 320, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 497, 185, 195, 320, 335, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 518, 215, 230, 350, 371, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 546, 230, 245, 371, 392, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 574, 245, 260, 392, 413, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 602, 260, 275, 413, 434, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 630, 275, 290, 434, 455, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 658, 290, 305, 455, 476, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 686, 305, 320, 476, 497, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 714, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 717, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 720, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 723, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 726, 2, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 735, 3, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 744, 4, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 753, 5, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 762, 6, 28, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 771, 16, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 789, 19, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 807, 22, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 825, 25, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 843, 28, 73, 79, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 861, 31, 79, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 879, 55, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 909, 61, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 939, 67, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 969, 73, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 999, 79, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 1029, 85, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1059, 125, 215, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1104, 135, 230, 245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1149, 145, 245, 260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1194, 155, 260, 275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1239, 165, 275, 290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1284, 175, 290, 305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1329, 230, 350, 371, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1392, 245, 371, 392, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1455, 260, 392, 413, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1518, 275, 413, 434, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1581, 290, 434, 455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1644, 305, 455, 476, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1707, 371, 518, 546, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1791, 392, 546, 574, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1875, 413, 574, 602, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1959, 434, 602, 630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2043, 455, 630, 658, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2127, 476, 658, 686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2211, 2, 3, 714, 717, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2217, 3, 4, 717, 720, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2223, 4, 5, 720, 723, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2229, 16, 19, 714, 726, 735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2247, 19, 22, 717, 735, 744, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2265, 22, 25, 720, 744, 753, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2283, 25, 28, 723, 753, 762, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2301, 49, 55, 726, 771, 789, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2337, 55, 61, 735, 789, 807, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2373, 61, 67, 744, 807, 825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2409, 67, 73, 753, 825, 843, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2445, 73, 79, 762, 843, 861, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2481, 115, 125, 789, 879, 909, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2541, 125, 135, 807, 909, 939, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2601, 135, 145, 825, 939, 969, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2661, 145, 155, 843, 969, 999, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2721, 155, 165, 861, 999, 1029, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2781, 215, 230, 909, 1059, 1104, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2871, 230, 245, 939, 1104, 1149, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2961, 245, 260, 969, 1149, 1194, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3051, 260, 275, 999, 1194, 1239, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3141, 275, 290, 1029, 1239, 1284, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3231, 350, 371, 1104, 1329, 1392, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3357, 371, 392, 1149, 1392, 1455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3483, 392, 413, 1194, 1455, 1518, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3609, 413, 434, 1239, 1518, 1581, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3735, 434, 455, 1284, 1581, 1644, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3861, 518, 546, 1392, 1707, 1791, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4029, 546, 574, 1455, 1791, 1875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4197, 574, 602, 1518, 1875, 1959, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4365, 602, 630, 1581, 1959, 2043, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4533, 630, 658, 1644, 2043, 2127, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 4701, 714, 717, 2211, 2217, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 4711, 717, 720, 2217, 2223, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 4721, 726, 735, 2211, 2229, 2247, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 4751, 735, 744, 2217, 2247, 2265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 4781, 744, 753, 2223, 2265, 2283, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4811, 771, 789, 2229, 2301, 2337, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4871, 789, 807, 2247, 2337, 2373, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4931, 807, 825, 2265, 2373, 2409, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4991, 825, 843, 2283, 2409, 2445, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 5051, 879, 909, 2337, 2481, 2541, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 5151, 909, 939, 2373, 2541, 2601, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 5251, 939, 969, 2409, 2601, 2661, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 5351, 969, 999, 2445, 2661, 2721, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5451, 1059, 1104, 2541, 2781, 2871, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5601, 1104, 1149, 2601, 2871, 2961, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5751, 1149, 1194, 2661, 2961, 3051, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5901, 1194, 1239, 2721, 3051, 3141, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 6051, 1329, 1392, 2871, 3231, 3357, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 6261, 1392, 1455, 2961, 3357, 3483, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 6471, 1455, 1518, 3051, 3483, 3609, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 6681, 1518, 1581, 3141, 3609, 3735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 6891, 1707, 1791, 3357, 3861, 4029, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 7171, 1791, 1875, 3483, 4029, 4197, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 7451, 1875, 1959, 3609, 4197, 4365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 7731, 1959, 2043, 3735, 4365, 4533, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 8011, 2211, 2217, 4701, 4711, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 8026, 2229, 2247, 4701, 4721, 4751, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 8071, 2247, 2265, 4711, 4751, 4781, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 8116, 2301, 2337, 4721, 4811, 4871, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 8206, 2337, 2373, 4751, 4871, 4931, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 8296, 2373, 2409, 4781, 4931, 4991, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 8386, 2481, 2541, 4871, 5051, 5151, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 8536, 2541, 2601, 4931, 5151, 5251, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 8686, 2601, 2661, 4991, 5251, 5351, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 8836, 2781, 2871, 5151, 5451, 5601, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 9061, 2871, 2961, 5251, 5601, 5751, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 9286, 2961, 3051, 5351, 5751, 5901, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 9511, 3231, 3357, 5601, 6051, 6261, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 9826, 3357, 3483, 5751, 6261, 6471, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 10141, 3483, 3609, 5901, 6471, 6681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 10456, 3861, 4029, 6261, 6891, 7171, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 10876, 4029, 4197, 6471, 7171, 7451, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 11296, 4197, 4365, 6681, 7451, 7731, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 11716, 4721, 4751, 8011, 8026, 8071, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 11779, 4811, 4871, 8026, 8116, 8206, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 11905, 4871, 4931, 8071, 8206, 8296, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 12031, 5051, 5151, 8206, 8386, 8536, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 12241, 5151, 5251, 8296, 8536, 8686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 12451, 5451, 5601, 8536, 8836, 9061, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 12766, 5601, 5751, 8686, 9061, 9286, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 13081, 6051, 6261, 9061, 9511, 9826, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 13522, 6261, 6471, 9286, 9826, 10141, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 13963, 6891, 7171, 9826, 10456, 10876, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 14551, 7171, 7451, 10141, 10876, 11296, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 15139, 8116, 8206, 11716, 11779, 11905, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 15307, 8386, 8536, 11905, 12031, 12241, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 15587, 8836, 9061, 12241, 12451, 12766, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 16007, 9511, 9826, 12766, 13081, 13522, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 16595, 10456, 10876, 13522, 13963, 14551, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 4811, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 5051, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 5451, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 310, pbuffer, 8116, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 400, pbuffer, 8386, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 550, pbuffer, 8836, 225, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {4811, 4871});

                pbuffer.scale(2.0 * a_exp, {5051, 5151});

                pbuffer.scale(2.0 * a_exp, {5451, 5601});

                pbuffer.scale(2.0 * a_exp, {8116, 8206});

                pbuffer.scale(2.0 * a_exp, {8386, 8536});

                pbuffer.scale(2.0 * a_exp, {8836, 9061});

                pbuffer.scale(2.0 * a_exp, {11779, 11905});

                pbuffer.scale(2.0 * a_exp, {12031, 12241});

                pbuffer.scale(2.0 * a_exp, {12451, 12766});

                pbuffer.scale(2.0 * a_exp, {15139, 15307});

                pbuffer.scale(2.0 * a_exp, {15307, 15587});

                pbuffer.scale(2.0 * a_exp, {15587, 16007});

                t2cfunc::reduce(cbuffer, 2775, pbuffer, 4811, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2835, pbuffer, 5051, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2935, pbuffer, 5451, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3085, pbuffer, 8116, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3175, pbuffer, 8386, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3325, pbuffer, 8836, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3550, pbuffer, 11779, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3676, pbuffer, 12031, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3886, pbuffer, 12451, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4201, pbuffer, 15139, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4369, pbuffer, 15307, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4649, pbuffer, 15587, 420, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {4811, 4871});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {5051, 5151});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {5451, 5601});

                pbuffer.scale(pfactors, 0, 2.0, {6051, 6261});

                pbuffer.scale(pfactors, 0, 2.0, {6891, 7171});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {8116, 8206});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {8386, 8536});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {8836, 9061});

                pbuffer.scale(pfactors, 0, 2.0, {9511, 9826});

                pbuffer.scale(pfactors, 0, 2.0, {10456, 10876});

                t2cfunc::reduce(cbuffer, 775, pbuffer, 4811, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 835, pbuffer, 5051, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 935, pbuffer, 5451, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1085, pbuffer, 6051, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1295, pbuffer, 6891, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1575, pbuffer, 8116, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1665, pbuffer, 8386, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1815, pbuffer, 8836, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2040, pbuffer, 9511, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2355, pbuffer, 10456, 420, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {4811, 4871});

                pbuffer.scale(2.0 * a_exp, {5051, 5151});

                pbuffer.scale(2.0 * a_exp, {5451, 5601});

                pbuffer.scale(2.0 * a_exp, {6051, 6261});

                pbuffer.scale(2.0 * a_exp, {6891, 7171});

                pbuffer.scale(2.0 * a_exp, {8116, 8206});

                pbuffer.scale(2.0 * a_exp, {8386, 8536});

                pbuffer.scale(2.0 * a_exp, {8836, 9061});

                pbuffer.scale(2.0 * a_exp, {9511, 9826});

                pbuffer.scale(2.0 * a_exp, {10456, 10876});

                pbuffer.scale(pfactors, 0, 2.0, {11779, 11905});

                pbuffer.scale(pfactors, 0, 2.0, {12031, 12241});

                pbuffer.scale(pfactors, 0, 2.0, {12451, 12766});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {13081, 13522});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {13963, 14551});

                pbuffer.scale(pfactors, 0, 2.0, {15139, 15307});

                pbuffer.scale(pfactors, 0, 2.0, {15307, 15587});

                pbuffer.scale(pfactors, 0, 2.0, {15587, 16007});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {16007, 16595});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {16595, 17379});

                t2cfunc::reduce(cbuffer, 5069, pbuffer, 4811, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5129, pbuffer, 5051, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5229, pbuffer, 5451, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5379, pbuffer, 6051, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5589, pbuffer, 6891, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5869, pbuffer, 8116, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5959, pbuffer, 8386, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6109, pbuffer, 8836, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6334, pbuffer, 9511, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6649, pbuffer, 10456, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7069, pbuffer, 11779, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7195, pbuffer, 12031, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7405, pbuffer, 12451, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7720, pbuffer, 13081, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8161, pbuffer, 13963, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8749, pbuffer, 15139, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8917, pbuffer, 15307, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9197, pbuffer, 15587, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9617, pbuffer, 16007, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10205, pbuffer, 16595, 784, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1560, cbuffer, 0, 60, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2280, cbuffer, 60, 160, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 4830, 1560, 2280, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 12210, cbuffer, 310, 400, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 13290, cbuffer, 400, 550, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 17115, 12210, 13290, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 26235, cbuffer, 2775, 2835, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 26955, cbuffer, 2835, 2935, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 29505, 26235, 26955, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 36885, cbuffer, 3085, 3175, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 37965, cbuffer, 3175, 3325, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 41790, 36885, 37965, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 52626, cbuffer, 3550, 3676, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 54138, cbuffer, 3676, 3886, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 59493, 52626, 54138, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 74445, cbuffer, 4201, 4369, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 76461, cbuffer, 4369, 4649, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 83601, 74445, 76461, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 0, cbuffer, 775, 835, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 180, cbuffer, 835, 935, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 480, cbuffer, 935, 1085, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 930, cbuffer, 1085, 1295, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1740, cbuffer, 0, 0, 180, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2580, cbuffer, 60, 180, 480, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 3480, cbuffer, 160, 480, 930, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 5190, 1560, 1740, 2580, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 6270, 2280, 2580, 3480, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 8070, 4830, 5190, 6270, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9870, cbuffer, 1575, 1665, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 10140, cbuffer, 1665, 1815, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 10590, cbuffer, 1815, 2040, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 11265, cbuffer, 2040, 2355, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 12480, cbuffer, 310, 9870, 10140, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 13740, cbuffer, 400, 10140, 10590, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 15090, cbuffer, 550, 10590, 11265, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 17655, 12210, 12480, 13740, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 19275, 13290, 13740, 15090, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 21975, 17115, 17655, 19275, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 24675, cbuffer, 5069, 5129, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 24855, cbuffer, 5129, 5229, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 25155, cbuffer, 5229, 5379, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 25605, cbuffer, 5379, 5589, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 26415, cbuffer, 2775, 24675, 24855, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 27255, cbuffer, 2835, 24855, 25155, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 28155, cbuffer, 2935, 25155, 25605, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 29865, 26235, 26415, 27255, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 30945, 26955, 27255, 28155, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 32745, 29505, 29865, 30945, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 34545, cbuffer, 5869, 5959, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 34815, cbuffer, 5959, 6109, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 35265, cbuffer, 6109, 6334, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 35940, cbuffer, 6334, 6649, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 37155, cbuffer, 3085, 34545, 34815, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 38415, cbuffer, 3175, 34815, 35265, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 39765, cbuffer, 3325, 35265, 35940, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 42330, 36885, 37155, 38415, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 43950, 37965, 38415, 39765, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 46650, 41790, 42330, 43950, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 49350, cbuffer, 7069, 7195, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 49728, cbuffer, 7195, 7405, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 50358, cbuffer, 7405, 7720, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 51303, cbuffer, 7720, 8161, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 53004, cbuffer, 3550, 49350, 49728, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 54768, cbuffer, 3676, 49728, 50358, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 56658, cbuffer, 3886, 50358, 51303, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 60249, 52626, 53004, 54768, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 62517, 54138, 54768, 56658, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 66297, 59493, 60249, 62517, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 70077, cbuffer, 8749, 8917, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 70581, cbuffer, 8917, 9197, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 71421, cbuffer, 9197, 9617, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 72681, cbuffer, 9617, 10205, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 74949, cbuffer, 4201, 70077, 70581, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 77301, cbuffer, 4369, 70581, 71421, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 79821, cbuffer, 4649, 71421, 72681, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 84609, 74445, 74949, 77301, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 87633, 76461, 77301, 79821, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 92673, 83601, 84609, 87633, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<3, 2>(skbuffer, 0, ckbuffer, 8070, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 350, ckbuffer, 8670, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 700, ckbuffer, 9270, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 4200, ckbuffer, 21975, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 4725, ckbuffer, 22875, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 5250, ckbuffer, 23775, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 62790, ckbuffer, 32745, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 63140, ckbuffer, 33345, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 63490, ckbuffer, 33945, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 63840, ckbuffer, 46650, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 64365, ckbuffer, 47550, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 64890, ckbuffer, 48450, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 65415, ckbuffer, 66297, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 66150, ckbuffer, 67557, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 66885, ckbuffer, 68817, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 67620, ckbuffer, 92673, 0, 6);

            t4cfunc::ket_transform<3, 2>(skbuffer, 68600, ckbuffer, 94353, 0, 6);

            t4cfunc::ket_transform<3, 2>(skbuffer, 69580, ckbuffer, 96033, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 17115, 0, 4200, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 18165, 350, 4725, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 19215, 700, 5250, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 1050, 62790, 63840, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 5775, 63840, 65415, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 10500, 65415, 67620, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 20265, 0, 1050, 5775, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 29715, 4200, 5775, 10500, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 43890, 17115, 20265, 29715, r_ab, 3, 2);

            t4cfunc::bra_transform<2, 3>(sbuffer, 0, skbuffer, 43890, 3, 2);

            t4cfunc::bra_transform<2, 3>(sbuffer, 1225, skbuffer, 45990, 3, 2);

            t4cfunc::bra_transform<2, 3>(sbuffer, 2450, skbuffer, 48090, 3, 2);

            t4cfunc::bra_transform<2, 3>(sbuffer, 3675, skbuffer, 50190, 3, 2);

            t4cfunc::bra_transform<2, 3>(sbuffer, 4900, skbuffer, 52290, 3, 2);

            t4cfunc::bra_transform<2, 3>(sbuffer, 6125, skbuffer, 54390, 3, 2);

            t4cfunc::bra_transform<2, 3>(sbuffer, 7350, skbuffer, 56490, 3, 2);

            t4cfunc::bra_transform<2, 3>(sbuffer, 8575, skbuffer, 58590, 3, 2);

            t4cfunc::bra_transform<2, 3>(sbuffer, 9800, skbuffer, 60690, 3, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 3, 3, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDFFD_hpp */
