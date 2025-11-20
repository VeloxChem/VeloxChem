#ifndef ElectronRepulsionGeom1010RecFDFD_hpp
#define ElectronRepulsionGeom1010RecFDFD_hpp

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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FD|1/|r-r'||FD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fdfd(T& distributor,
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

    CSimdArray<double> cbuffer(12321, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(109557, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(115290, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 2301, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 2481, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 96, pbuffer, 2781, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 186, pbuffer, 4811, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 246, pbuffer, 5051, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 346, pbuffer, 5451, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 496, pbuffer, 8116, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 586, pbuffer, 8386, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 736, pbuffer, 8836, 225, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {2301, 2337});

                pbuffer.scale(2.0 * a_exp, {2481, 2541});

                pbuffer.scale(2.0 * a_exp, {2781, 2871});

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

                t2cfunc::reduce(cbuffer, 3441, pbuffer, 2301, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3477, pbuffer, 2481, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3537, pbuffer, 2781, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3627, pbuffer, 4811, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3687, pbuffer, 5051, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3787, pbuffer, 5451, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3937, pbuffer, 8116, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4027, pbuffer, 8386, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4177, pbuffer, 8836, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4402, pbuffer, 11779, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4528, pbuffer, 12031, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4738, pbuffer, 12451, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5053, pbuffer, 15139, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5221, pbuffer, 15307, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5501, pbuffer, 15587, 420, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2301, 2337});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2481, 2541});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2781, 2871});

                pbuffer.scale(pfactors, 0, 2.0, {3231, 3357});

                pbuffer.scale(pfactors, 0, 2.0, {3861, 4029});

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

                t2cfunc::reduce(cbuffer, 961, pbuffer, 2301, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 997, pbuffer, 2481, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1057, pbuffer, 2781, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1147, pbuffer, 3231, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1273, pbuffer, 3861, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1441, pbuffer, 4811, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1501, pbuffer, 5051, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1601, pbuffer, 5451, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1751, pbuffer, 6051, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1961, pbuffer, 6891, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2241, pbuffer, 8116, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2331, pbuffer, 8386, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2481, pbuffer, 8836, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2706, pbuffer, 9511, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3021, pbuffer, 10456, 420, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {2301, 2337});

                pbuffer.scale(2.0 * a_exp, {2481, 2541});

                pbuffer.scale(2.0 * a_exp, {2781, 2871});

                pbuffer.scale(2.0 * a_exp, {3231, 3357});

                pbuffer.scale(2.0 * a_exp, {3861, 4029});

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

                t2cfunc::reduce(cbuffer, 5921, pbuffer, 2301, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5957, pbuffer, 2481, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6017, pbuffer, 2781, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6107, pbuffer, 3231, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6233, pbuffer, 3861, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6401, pbuffer, 4811, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6461, pbuffer, 5051, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6561, pbuffer, 5451, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6711, pbuffer, 6051, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6921, pbuffer, 6891, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7201, pbuffer, 8116, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7291, pbuffer, 8386, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7441, pbuffer, 8836, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7666, pbuffer, 9511, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7981, pbuffer, 10456, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8401, pbuffer, 11779, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8527, pbuffer, 12031, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8737, pbuffer, 12451, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9052, pbuffer, 13081, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9493, pbuffer, 13963, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10081, pbuffer, 15139, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10249, pbuffer, 15307, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10529, pbuffer, 15587, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10949, pbuffer, 16007, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11537, pbuffer, 16595, 784, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 936, cbuffer, 0, 36, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1368, cbuffer, 36, 96, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 2898, 936, 1368, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 7482, cbuffer, 186, 246, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 8202, cbuffer, 246, 346, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 10752, 7482, 8202, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 18132, cbuffer, 496, 586, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 19212, cbuffer, 586, 736, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 23037, 18132, 19212, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 31533, cbuffer, 3441, 3477, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 31965, cbuffer, 3477, 3537, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 33495, 31533, 31965, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 38079, cbuffer, 3627, 3687, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 38799, cbuffer, 3687, 3787, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 41349, 38079, 38799, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 48729, cbuffer, 3937, 4027, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 49809, cbuffer, 4027, 4177, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 53634, 48729, 49809, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 64470, cbuffer, 4402, 4528, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 65982, cbuffer, 4528, 4738, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 71337, 64470, 65982, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 86289, cbuffer, 5053, 5221, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 88305, cbuffer, 5221, 5501, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 95445, 86289, 88305, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 0, cbuffer, 961, 997, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 108, cbuffer, 997, 1057, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 288, cbuffer, 1057, 1147, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 558, cbuffer, 1147, 1273, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1044, cbuffer, 0, 0, 108, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 1548, cbuffer, 36, 108, 288, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 2088, cbuffer, 96, 288, 558, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 3114, 936, 1044, 1548, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 3762, 1368, 1548, 2088, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 4842, 2898, 3114, 3762, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5922, cbuffer, 1441, 1501, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6102, cbuffer, 1501, 1601, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 6402, cbuffer, 1601, 1751, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 6852, cbuffer, 1751, 1961, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 7662, cbuffer, 186, 5922, 6102, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 8502, cbuffer, 246, 6102, 6402, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 9402, cbuffer, 346, 6402, 6852, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 11112, 7482, 7662, 8502, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 12192, 8202, 8502, 9402, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 13992, 10752, 11112, 12192, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 15792, cbuffer, 2241, 2331, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 16062, cbuffer, 2331, 2481, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 16512, cbuffer, 2481, 2706, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 17187, cbuffer, 2706, 3021, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 18402, cbuffer, 496, 15792, 16062, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 19662, cbuffer, 586, 16062, 16512, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 21012, cbuffer, 736, 16512, 17187, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 23577, 18132, 18402, 19662, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 25197, 19212, 19662, 21012, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 27897, 23037, 23577, 25197, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 30597, cbuffer, 5921, 5957, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 30705, cbuffer, 5957, 6017, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 30885, cbuffer, 6017, 6107, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 31155, cbuffer, 6107, 6233, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 31641, cbuffer, 3441, 30597, 30705, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 32145, cbuffer, 3477, 30705, 30885, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 32685, cbuffer, 3537, 30885, 31155, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 33711, 31533, 31641, 32145, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 34359, 31965, 32145, 32685, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 35439, 33495, 33711, 34359, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 36519, cbuffer, 6401, 6461, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 36699, cbuffer, 6461, 6561, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 36999, cbuffer, 6561, 6711, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 37449, cbuffer, 6711, 6921, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 38259, cbuffer, 3627, 36519, 36699, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 39099, cbuffer, 3687, 36699, 36999, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 39999, cbuffer, 3787, 36999, 37449, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 41709, 38079, 38259, 39099, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 42789, 38799, 39099, 39999, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 44589, 41349, 41709, 42789, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 46389, cbuffer, 7201, 7291, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 46659, cbuffer, 7291, 7441, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 47109, cbuffer, 7441, 7666, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 47784, cbuffer, 7666, 7981, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 48999, cbuffer, 3937, 46389, 46659, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 50259, cbuffer, 4027, 46659, 47109, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 51609, cbuffer, 4177, 47109, 47784, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 54174, 48729, 48999, 50259, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 55794, 49809, 50259, 51609, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 58494, 53634, 54174, 55794, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 61194, cbuffer, 8401, 8527, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 61572, cbuffer, 8527, 8737, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 62202, cbuffer, 8737, 9052, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 63147, cbuffer, 9052, 9493, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 64848, cbuffer, 4402, 61194, 61572, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 66612, cbuffer, 4528, 61572, 62202, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 68502, cbuffer, 4738, 62202, 63147, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 72093, 64470, 64848, 66612, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 74361, 65982, 66612, 68502, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 78141, 71337, 72093, 74361, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 81921, cbuffer, 10081, 10249, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 82425, cbuffer, 10249, 10529, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 83265, cbuffer, 10529, 10949, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 84525, cbuffer, 10949, 11537, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 86793, cbuffer, 5053, 81921, 82425, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 89145, cbuffer, 5221, 82425, 83265, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 91665, cbuffer, 5501, 83265, 84525, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 96453, 86289, 86793, 89145, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 99477, 88305, 89145, 91665, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 104517, 95445, 96453, 99477, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<3, 2>(skbuffer, 0, ckbuffer, 4842, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 210, ckbuffer, 5202, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 420, ckbuffer, 5562, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 2520, ckbuffer, 13992, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 2870, ckbuffer, 14592, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 3220, ckbuffer, 15192, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 6720, ckbuffer, 27897, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 7245, ckbuffer, 28797, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 7770, ckbuffer, 29697, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 106890, ckbuffer, 35439, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 107100, ckbuffer, 35799, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 107310, ckbuffer, 36159, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 107520, ckbuffer, 44589, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 107870, ckbuffer, 45189, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 108220, ckbuffer, 45789, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 108570, ckbuffer, 58494, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 109095, ckbuffer, 59394, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 109620, ckbuffer, 60294, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 110145, ckbuffer, 78141, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 110880, ckbuffer, 79401, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 111615, ckbuffer, 80661, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 112350, ckbuffer, 104517, 0, 6);

            t4cfunc::ket_transform<3, 2>(skbuffer, 113330, ckbuffer, 106197, 0, 6);

            t4cfunc::ket_transform<3, 2>(skbuffer, 114310, ckbuffer, 107877, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 19635, 0, 2520, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 20265, 210, 2870, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 20895, 420, 3220, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 27195, 2520, 6720, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 28245, 2870, 7245, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 29295, 3220, 7770, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 53970, 19635, 27195, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 55230, 20265, 28245, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 56490, 20895, 29295, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 630, 106890, 107520, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 3570, 107520, 108570, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 8295, 108570, 110145, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 13020, 110145, 112350, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 21525, 0, 630, 3570, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 30345, 2520, 3570, 8295, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 39795, 6720, 8295, 13020, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 57750, 19635, 21525, 30345, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 69090, 27195, 30345, 39795, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fdxx(skbuffer, 87990, 53970, 57750, 69090, r_ab, 3, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 87990, 3, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1225, skbuffer, 90090, 3, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2450, skbuffer, 92190, 3, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 3675, skbuffer, 94290, 3, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 4900, skbuffer, 96390, 3, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 6125, skbuffer, 98490, 3, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 7350, skbuffer, 100590, 3, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 8575, skbuffer, 102690, 3, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 9800, skbuffer, 104790, 3, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 3, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFDFD_hpp */
