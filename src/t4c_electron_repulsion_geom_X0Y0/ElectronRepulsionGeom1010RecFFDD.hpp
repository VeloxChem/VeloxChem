#ifndef ElectronRepulsionGeom1010RecFFDD_hpp
#define ElectronRepulsionGeom1010RecFFDD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"
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
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSF.hpp"
#include "ElectronRepulsionPrimRecSKSG.hpp"
#include "ElectronRepulsionPrimRecSKSH.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSH.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FF|1/|r-r'||DD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_ffdd(T& distributor,
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

    CSimdArray<double> pbuffer(17991, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(10608, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(56628, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(125775, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 518, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 521, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 524, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 527, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 530, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 533, 2, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 542, 3, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 551, 4, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 560, 5, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 569, 6, 28, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 578, 7, 31, 34, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 587, 16, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 605, 19, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 623, 22, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 641, 25, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 659, 28, 73, 79, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 677, 31, 79, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 695, 34, 85, 91, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 713, 55, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 743, 61, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 773, 67, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 803, 73, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 833, 79, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 863, 85, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 893, 91, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 923, 125, 215, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 968, 135, 230, 245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1013, 145, 245, 260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1058, 155, 260, 275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1103, 165, 275, 290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1148, 175, 290, 305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1193, 185, 305, 320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1238, 230, 350, 371, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1301, 245, 371, 392, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1364, 260, 392, 413, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1427, 275, 413, 434, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1490, 290, 434, 455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1553, 305, 455, 476, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1616, 320, 476, 497, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1679, 2, 3, 518, 521, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1685, 3, 4, 521, 524, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1691, 4, 5, 524, 527, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1697, 5, 6, 527, 530, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1703, 16, 19, 518, 533, 542, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1721, 19, 22, 521, 542, 551, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1739, 22, 25, 524, 551, 560, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1757, 25, 28, 527, 560, 569, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1775, 28, 31, 530, 569, 578, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1793, 49, 55, 533, 587, 605, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1829, 55, 61, 542, 605, 623, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1865, 61, 67, 551, 623, 641, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1901, 67, 73, 560, 641, 659, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1937, 73, 79, 569, 659, 677, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1973, 79, 85, 578, 677, 695, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2009, 115, 125, 605, 713, 743, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2069, 125, 135, 623, 743, 773, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2129, 135, 145, 641, 773, 803, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2189, 145, 155, 659, 803, 833, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2249, 155, 165, 677, 833, 863, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2309, 165, 175, 695, 863, 893, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2369, 215, 230, 743, 923, 968, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2459, 230, 245, 773, 968, 1013, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2549, 245, 260, 803, 1013, 1058, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2639, 260, 275, 833, 1058, 1103, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2729, 275, 290, 863, 1103, 1148, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2819, 290, 305, 893, 1148, 1193, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2909, 350, 371, 968, 1238, 1301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3035, 371, 392, 1013, 1301, 1364, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3161, 392, 413, 1058, 1364, 1427, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3287, 413, 434, 1103, 1427, 1490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3413, 434, 455, 1148, 1490, 1553, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3539, 455, 476, 1193, 1553, 1616, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3665, 518, 521, 1679, 1685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3675, 521, 524, 1685, 1691, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3685, 524, 527, 1691, 1697, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3695, 533, 542, 1679, 1703, 1721, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3725, 542, 551, 1685, 1721, 1739, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3755, 551, 560, 1691, 1739, 1757, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3785, 560, 569, 1697, 1757, 1775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3815, 587, 605, 1703, 1793, 1829, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3875, 605, 623, 1721, 1829, 1865, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3935, 623, 641, 1739, 1865, 1901, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3995, 641, 659, 1757, 1901, 1937, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4055, 659, 677, 1775, 1937, 1973, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4115, 713, 743, 1829, 2009, 2069, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4215, 743, 773, 1865, 2069, 2129, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4315, 773, 803, 1901, 2129, 2189, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4415, 803, 833, 1937, 2189, 2249, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4515, 833, 863, 1973, 2249, 2309, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4615, 923, 968, 2069, 2369, 2459, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4765, 968, 1013, 2129, 2459, 2549, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4915, 1013, 1058, 2189, 2549, 2639, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5065, 1058, 1103, 2249, 2639, 2729, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5215, 1103, 1148, 2309, 2729, 2819, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5365, 1238, 1301, 2459, 2909, 3035, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5575, 1301, 1364, 2549, 3035, 3161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5785, 1364, 1427, 2639, 3161, 3287, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5995, 1427, 1490, 2729, 3287, 3413, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 6205, 1490, 1553, 2819, 3413, 3539, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 6415, 1679, 1685, 3665, 3675, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 6430, 1685, 1691, 3675, 3685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 6445, 1703, 1721, 3665, 3695, 3725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 6490, 1721, 1739, 3675, 3725, 3755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 6535, 1739, 1757, 3685, 3755, 3785, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 6580, 1793, 1829, 3695, 3815, 3875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 6670, 1829, 1865, 3725, 3875, 3935, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 6760, 1865, 1901, 3755, 3935, 3995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 6850, 1901, 1937, 3785, 3995, 4055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 6940, 2009, 2069, 3875, 4115, 4215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 7090, 2069, 2129, 3935, 4215, 4315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 7240, 2129, 2189, 3995, 4315, 4415, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 7390, 2189, 2249, 4055, 4415, 4515, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 7540, 2369, 2459, 4215, 4615, 4765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 7765, 2459, 2549, 4315, 4765, 4915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 7990, 2549, 2639, 4415, 4915, 5065, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 8215, 2639, 2729, 4515, 5065, 5215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 8440, 2909, 3035, 4765, 5365, 5575, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 8755, 3035, 3161, 4915, 5575, 5785, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 9070, 3161, 3287, 5065, 5785, 5995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 9385, 3287, 3413, 5215, 5995, 6205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 9700, 3665, 3675, 6415, 6430, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 9721, 3695, 3725, 6415, 6445, 6490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 9784, 3725, 3755, 6430, 6490, 6535, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 9847, 3815, 3875, 6445, 6580, 6670, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 9973, 3875, 3935, 6490, 6670, 6760, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 10099, 3935, 3995, 6535, 6760, 6850, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 10225, 4115, 4215, 6670, 6940, 7090, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 10435, 4215, 4315, 6760, 7090, 7240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 10645, 4315, 4415, 6850, 7240, 7390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 10855, 4615, 4765, 7090, 7540, 7765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 11170, 4765, 4915, 7240, 7765, 7990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 11485, 4915, 5065, 7390, 7990, 8215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 11800, 5365, 5575, 7765, 8440, 8755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 12241, 5575, 5785, 7990, 8755, 9070, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 12682, 5785, 5995, 8215, 9070, 9385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 13123, 6445, 6490, 9700, 9721, 9784, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 13207, 6580, 6670, 9721, 9847, 9973, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 13375, 6670, 6760, 9784, 9973, 10099, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 13543, 6940, 7090, 9973, 10225, 10435, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 13823, 7090, 7240, 10099, 10435, 10645, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 14103, 7540, 7765, 10435, 10855, 11170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 14523, 7765, 7990, 10645, 11170, 11485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 14943, 8440, 8755, 11170, 11800, 12241, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 15531, 8755, 9070, 11485, 12241, 12682, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 16119, 9847, 9973, 13123, 13207, 13375, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 16335, 10225, 10435, 13375, 13543, 13823, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 16695, 10855, 11170, 13823, 14103, 14523, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksh(pbuffer, 17235, 11800, 12241, 14523, 14943, 15531, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 3815, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 4115, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 6580, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 250, pbuffer, 6940, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 400, pbuffer, 9847, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 526, pbuffer, 10225, 210, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {3815, 3875});

                pbuffer.scale(2.0 * a_exp, {4115, 4215});

                pbuffer.scale(2.0 * a_exp, {6580, 6670});

                pbuffer.scale(2.0 * a_exp, {6940, 7090});

                pbuffer.scale(2.0 * a_exp, {9847, 9973});

                pbuffer.scale(2.0 * a_exp, {10225, 10435});

                pbuffer.scale(2.0 * a_exp, {13207, 13375});

                pbuffer.scale(2.0 * a_exp, {13543, 13823});

                pbuffer.scale(2.0 * a_exp, {16119, 16335});

                pbuffer.scale(2.0 * a_exp, {16335, 16695});

                t2cfunc::reduce(cbuffer, 3128, pbuffer, 3815, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3188, pbuffer, 4115, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3288, pbuffer, 6580, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3378, pbuffer, 6940, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3528, pbuffer, 9847, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3654, pbuffer, 10225, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3864, pbuffer, 13207, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4032, pbuffer, 13543, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4312, pbuffer, 16119, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4528, pbuffer, 16335, 360, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3815, 3875});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {4115, 4215});

                pbuffer.scale(pfactors, 0, 2.0, {4615, 4765});

                pbuffer.scale(pfactors, 0, 2.0, {5365, 5575});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {6580, 6670});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {6940, 7090});

                pbuffer.scale(pfactors, 0, 2.0, {7540, 7765});

                pbuffer.scale(pfactors, 0, 2.0, {8440, 8755});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {9847, 9973});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {10225, 10435});

                pbuffer.scale(pfactors, 0, 2.0, {10855, 11170});

                pbuffer.scale(pfactors, 0, 2.0, {11800, 12241});

                t2cfunc::reduce(cbuffer, 736, pbuffer, 3815, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 796, pbuffer, 4115, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 896, pbuffer, 4615, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1046, pbuffer, 5365, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1256, pbuffer, 6580, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1346, pbuffer, 6940, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1496, pbuffer, 7540, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1721, pbuffer, 8440, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2036, pbuffer, 9847, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2162, pbuffer, 10225, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2372, pbuffer, 10855, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2687, pbuffer, 11800, 441, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {3815, 3875});

                pbuffer.scale(2.0 * a_exp, {4115, 4215});

                pbuffer.scale(2.0 * a_exp, {4615, 4765});

                pbuffer.scale(2.0 * a_exp, {5365, 5575});

                pbuffer.scale(2.0 * a_exp, {6580, 6670});

                pbuffer.scale(2.0 * a_exp, {6940, 7090});

                pbuffer.scale(2.0 * a_exp, {7540, 7765});

                pbuffer.scale(2.0 * a_exp, {8440, 8755});

                pbuffer.scale(2.0 * a_exp, {9847, 9973});

                pbuffer.scale(2.0 * a_exp, {10225, 10435});

                pbuffer.scale(2.0 * a_exp, {10855, 11170});

                pbuffer.scale(2.0 * a_exp, {11800, 12241});

                pbuffer.scale(pfactors, 0, 2.0, {13207, 13375});

                pbuffer.scale(pfactors, 0, 2.0, {13543, 13823});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {14103, 14523});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {14943, 15531});

                pbuffer.scale(pfactors, 0, 2.0, {16119, 16335});

                pbuffer.scale(pfactors, 0, 2.0, {16335, 16695});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {16695, 17235});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {17235, 17991});

                t2cfunc::reduce(cbuffer, 4888, pbuffer, 3815, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4948, pbuffer, 4115, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5048, pbuffer, 4615, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5198, pbuffer, 5365, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5408, pbuffer, 6580, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5498, pbuffer, 6940, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5648, pbuffer, 7540, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5873, pbuffer, 8440, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6188, pbuffer, 9847, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6314, pbuffer, 10225, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6524, pbuffer, 10855, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6839, pbuffer, 11800, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7280, pbuffer, 13207, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7448, pbuffer, 13543, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7728, pbuffer, 14103, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8148, pbuffer, 14943, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8736, pbuffer, 16119, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8952, pbuffer, 16335, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9312, pbuffer, 16695, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9852, pbuffer, 17235, 756, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 930, cbuffer, 0, 60, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5025, cbuffer, 160, 250, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 11028, cbuffer, 400, 526, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 17628, cbuffer, 3128, 3188, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 21723, cbuffer, 3288, 3378, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 27726, cbuffer, 3528, 3654, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 36000, cbuffer, 3864, 4032, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 46908, cbuffer, 4312, 4528, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 0, cbuffer, 736, 796, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 180, cbuffer, 796, 896, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 480, cbuffer, 896, 1046, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1110, cbuffer, 0, 0, 180, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 1650, cbuffer, 60, 180, 480, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 2550, 930, 1110, 1650, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3630, cbuffer, 1256, 1346, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3900, cbuffer, 1346, 1496, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 4350, cbuffer, 1496, 1721, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 5295, cbuffer, 160, 3630, 3900, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 6105, cbuffer, 250, 3900, 4350, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 7455, 5025, 5295, 6105, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9075, cbuffer, 2036, 2162, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 9453, cbuffer, 2162, 2372, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 10083, cbuffer, 2372, 2687, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 11406, cbuffer, 400, 9075, 9453, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 12540, cbuffer, 526, 9453, 10083, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 14430, 11028, 11406, 12540, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 16698, cbuffer, 4888, 4948, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 16878, cbuffer, 4948, 5048, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 17178, cbuffer, 5048, 5198, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 17808, cbuffer, 3128, 16698, 16878, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 18348, cbuffer, 3188, 16878, 17178, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 19248, 17628, 17808, 18348, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 20328, cbuffer, 5408, 5498, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 20598, cbuffer, 5498, 5648, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 21048, cbuffer, 5648, 5873, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 21993, cbuffer, 3288, 20328, 20598, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 22803, cbuffer, 3378, 20598, 21048, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 24153, 21723, 21993, 22803, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 25773, cbuffer, 6188, 6314, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 26151, cbuffer, 6314, 6524, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 26781, cbuffer, 6524, 6839, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 28104, cbuffer, 3528, 25773, 26151, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 29238, cbuffer, 3654, 26151, 26781, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 31128, 27726, 28104, 29238, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 33396, cbuffer, 7280, 7448, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 33900, cbuffer, 7448, 7728, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 34740, cbuffer, 7728, 8148, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 36504, cbuffer, 3864, 33396, 33900, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 38016, cbuffer, 4032, 33900, 34740, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 40536, 36000, 36504, 38016, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 43560, cbuffer, 8736, 8952, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 44208, cbuffer, 8952, 9312, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 45288, cbuffer, 9312, 9852, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 47556, cbuffer, 4312, 43560, 44208, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 49500, cbuffer, 4528, 44208, 45288, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 52740, 46908, 47556, 49500, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<2, 2>(skbuffer, 0, ckbuffer, 2550, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 250, ckbuffer, 2910, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 500, ckbuffer, 3270, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 3000, ckbuffer, 7455, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 3375, ckbuffer, 7995, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 3750, ckbuffer, 8535, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 7500, ckbuffer, 14430, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 8025, ckbuffer, 15186, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 8550, ckbuffer, 15942, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 117525, ckbuffer, 19248, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 117775, ckbuffer, 19608, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 118025, ckbuffer, 19968, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 118275, ckbuffer, 24153, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 118650, ckbuffer, 24693, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 119025, ckbuffer, 25233, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 119400, ckbuffer, 31128, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 119925, ckbuffer, 31884, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 120450, ckbuffer, 32640, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 120975, ckbuffer, 40536, 0, 6);

            t4cfunc::ket_transform<2, 2>(skbuffer, 121675, ckbuffer, 41544, 0, 6);

            t4cfunc::ket_transform<2, 2>(skbuffer, 122375, ckbuffer, 42552, 0, 6);

            t4cfunc::ket_transform<2, 2>(skbuffer, 123075, ckbuffer, 52740, 0, 7);

            t4cfunc::ket_transform<2, 2>(skbuffer, 123975, ckbuffer, 54036, 0, 7);

            t4cfunc::ket_transform<2, 2>(skbuffer, 124875, ckbuffer, 55332, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 20100, 0, 3000, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 20850, 250, 3375, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 21600, 500, 3750, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 29100, 3000, 7500, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 30225, 3375, 8025, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 31350, 3750, 8550, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 56775, 20100, 29100, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 58275, 20850, 30225, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 59775, 21600, 31350, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 750, 117525, 118275, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 4125, 118275, 119400, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 9075, 119400, 120975, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sixx(skbuffer, 13800, 120975, 123075, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 22350, 0, 750, 4125, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 32475, 3000, 4125, 9075, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_phxx(skbuffer, 42600, 7500, 9075, 13800, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 61275, 20100, 22350, 32475, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dgxx(skbuffer, 74775, 29100, 32475, 42600, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ffxx(skbuffer, 95025, 56775, 61275, 74775, r_ab, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 95025, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1225, skbuffer, 97525, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2450, skbuffer, 100025, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 3675, skbuffer, 102525, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 4900, skbuffer, 105025, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 6125, skbuffer, 107525, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 7350, skbuffer, 110025, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 8575, skbuffer, 112525, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 9800, skbuffer, 115025, 2, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 2, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFFDD_hpp */
