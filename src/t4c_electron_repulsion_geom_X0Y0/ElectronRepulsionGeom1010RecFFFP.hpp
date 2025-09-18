#ifndef ElectronRepulsionGeom1010RecFFFP_hpp
#define ElectronRepulsionGeom1010RecFFFP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXDP.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXFP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
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
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSF.hpp"
#include "ElectronRepulsionPrimRecSKSG.hpp"
#include "ElectronRepulsionPrimRecSKSH.hpp"
#include "ElectronRepulsionPrimRecSKSP.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FF|1/|r-r'||FP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fffp(T& distributor,
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

    CSimdArray<double> pbuffer(18431, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(11544, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(88920, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(105651, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(9261, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 518, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 521, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 524, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 527, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 530, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 533, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 536, 1, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 545, 2, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 554, 3, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 563, 4, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 572, 5, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 581, 6, 28, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 590, 7, 31, 34, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 599, 16, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 617, 19, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 635, 22, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 653, 25, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 671, 28, 73, 79, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 689, 31, 79, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 707, 34, 85, 91, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 725, 55, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 755, 61, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 785, 67, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 815, 73, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 845, 79, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 875, 85, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 905, 91, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 935, 125, 215, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 980, 135, 230, 245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1025, 145, 245, 260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1070, 155, 260, 275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1115, 165, 275, 290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1160, 175, 290, 305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1205, 185, 305, 320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1250, 230, 350, 371, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1313, 245, 371, 392, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1376, 260, 392, 413, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1439, 275, 413, 434, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1502, 290, 434, 455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1565, 305, 455, 476, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1628, 320, 476, 497, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1691, 1, 2, 518, 521, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1697, 2, 3, 521, 524, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1703, 3, 4, 524, 527, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1709, 4, 5, 527, 530, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1715, 5, 6, 530, 533, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1721, 13, 16, 518, 536, 545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1739, 16, 19, 521, 545, 554, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1757, 19, 22, 524, 554, 563, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1775, 22, 25, 527, 563, 572, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1793, 25, 28, 530, 572, 581, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1811, 28, 31, 533, 581, 590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1829, 49, 55, 545, 599, 617, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1865, 55, 61, 554, 617, 635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1901, 61, 67, 563, 635, 653, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1937, 67, 73, 572, 653, 671, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1973, 73, 79, 581, 671, 689, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2009, 79, 85, 590, 689, 707, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2045, 115, 125, 617, 725, 755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2105, 125, 135, 635, 755, 785, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2165, 135, 145, 653, 785, 815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2225, 145, 155, 671, 815, 845, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2285, 155, 165, 689, 845, 875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2345, 165, 175, 707, 875, 905, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2405, 215, 230, 755, 935, 980, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2495, 230, 245, 785, 980, 1025, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2585, 245, 260, 815, 1025, 1070, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2675, 260, 275, 845, 1070, 1115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2765, 275, 290, 875, 1115, 1160, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2855, 290, 305, 905, 1160, 1205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2945, 350, 371, 980, 1250, 1313, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3071, 371, 392, 1025, 1313, 1376, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3197, 392, 413, 1070, 1376, 1439, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3323, 413, 434, 1115, 1439, 1502, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3449, 434, 455, 1160, 1502, 1565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3575, 455, 476, 1205, 1565, 1628, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3701, 518, 521, 1691, 1697, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3711, 521, 524, 1697, 1703, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3721, 524, 527, 1703, 1709, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3731, 527, 530, 1709, 1715, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3741, 536, 545, 1691, 1721, 1739, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3771, 545, 554, 1697, 1739, 1757, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3801, 554, 563, 1703, 1757, 1775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3831, 563, 572, 1709, 1775, 1793, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3861, 572, 581, 1715, 1793, 1811, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3891, 599, 617, 1739, 1829, 1865, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3951, 617, 635, 1757, 1865, 1901, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4011, 635, 653, 1775, 1901, 1937, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4071, 653, 671, 1793, 1937, 1973, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4131, 671, 689, 1811, 1973, 2009, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4191, 725, 755, 1865, 2045, 2105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4291, 755, 785, 1901, 2105, 2165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4391, 785, 815, 1937, 2165, 2225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4491, 815, 845, 1973, 2225, 2285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4591, 845, 875, 2009, 2285, 2345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4691, 935, 980, 2105, 2405, 2495, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4841, 980, 1025, 2165, 2495, 2585, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4991, 1025, 1070, 2225, 2585, 2675, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5141, 1070, 1115, 2285, 2675, 2765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5291, 1115, 1160, 2345, 2765, 2855, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5441, 1250, 1313, 2495, 2945, 3071, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5651, 1313, 1376, 2585, 3071, 3197, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5861, 1376, 1439, 2675, 3197, 3323, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 6071, 1439, 1502, 2765, 3323, 3449, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 6281, 1502, 1565, 2855, 3449, 3575, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 6491, 1691, 1697, 3701, 3711, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 6506, 1697, 1703, 3711, 3721, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 6521, 1703, 1709, 3721, 3731, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 6536, 1721, 1739, 3701, 3741, 3771, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 6581, 1739, 1757, 3711, 3771, 3801, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 6626, 1757, 1775, 3721, 3801, 3831, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 6671, 1775, 1793, 3731, 3831, 3861, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 6716, 1829, 1865, 3771, 3891, 3951, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 6806, 1865, 1901, 3801, 3951, 4011, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 6896, 1901, 1937, 3831, 4011, 4071, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 6986, 1937, 1973, 3861, 4071, 4131, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 7076, 2045, 2105, 3951, 4191, 4291, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 7226, 2105, 2165, 4011, 4291, 4391, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 7376, 2165, 2225, 4071, 4391, 4491, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 7526, 2225, 2285, 4131, 4491, 4591, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 7676, 2405, 2495, 4291, 4691, 4841, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 7901, 2495, 2585, 4391, 4841, 4991, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 8126, 2585, 2675, 4491, 4991, 5141, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 8351, 2675, 2765, 4591, 5141, 5291, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 8576, 2945, 3071, 4841, 5441, 5651, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 8891, 3071, 3197, 4991, 5651, 5861, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 9206, 3197, 3323, 5141, 5861, 6071, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 9521, 3323, 3449, 5291, 6071, 6281, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 9836, 3701, 3711, 6491, 6506, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 9857, 3711, 3721, 6506, 6521, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 9878, 3741, 3771, 6491, 6536, 6581, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 9941, 3771, 3801, 6506, 6581, 6626, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 10004, 3801, 3831, 6521, 6626, 6671, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 10067, 3891, 3951, 6581, 6716, 6806, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 10193, 3951, 4011, 6626, 6806, 6896, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 10319, 4011, 4071, 6671, 6896, 6986, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 10445, 4191, 4291, 6806, 7076, 7226, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 10655, 4291, 4391, 6896, 7226, 7376, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 10865, 4391, 4491, 6986, 7376, 7526, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 11075, 4691, 4841, 7226, 7676, 7901, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 11390, 4841, 4991, 7376, 7901, 8126, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 11705, 4991, 5141, 7526, 8126, 8351, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 12020, 5441, 5651, 7901, 8576, 8891, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 12461, 5651, 5861, 8126, 8891, 9206, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 12902, 5861, 6071, 8351, 9206, 9521, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 13343, 6491, 6506, 9836, 9857, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 13371, 6536, 6581, 9836, 9878, 9941, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 13455, 6581, 6626, 9857, 9941, 10004, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 13539, 6716, 6806, 9941, 10067, 10193, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 13707, 6806, 6896, 10004, 10193, 10319, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 13875, 7076, 7226, 10193, 10445, 10655, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 14155, 7226, 7376, 10319, 10655, 10865, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 14435, 7676, 7901, 10655, 11075, 11390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 14855, 7901, 8126, 10865, 11390, 11705, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 15275, 8576, 8891, 11390, 12020, 12461, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 15863, 8891, 9206, 11705, 12461, 12902, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 16451, 9878, 9941, 13343, 13371, 13455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 16559, 10067, 10193, 13455, 13539, 13707, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 16775, 10445, 10655, 13707, 13875, 14155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 17135, 11075, 11390, 14155, 14435, 14855, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksh(pbuffer, 17675, 12020, 12461, 14855, 15275, 15863, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 3741, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 3891, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 4191, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 190, pbuffer, 6536, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 235, pbuffer, 6716, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 325, pbuffer, 7076, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 475, pbuffer, 9878, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 538, pbuffer, 10067, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 664, pbuffer, 10445, 210, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {3741, 3771});

                pbuffer.scale(2.0 * a_exp, {3891, 3951});

                pbuffer.scale(2.0 * a_exp, {4191, 4291});

                pbuffer.scale(2.0 * a_exp, {6536, 6581});

                pbuffer.scale(2.0 * a_exp, {6716, 6806});

                pbuffer.scale(2.0 * a_exp, {7076, 7226});

                pbuffer.scale(2.0 * a_exp, {9878, 9941});

                pbuffer.scale(2.0 * a_exp, {10067, 10193});

                pbuffer.scale(2.0 * a_exp, {10445, 10655});

                pbuffer.scale(2.0 * a_exp, {13371, 13455});

                pbuffer.scale(2.0 * a_exp, {13539, 13707});

                pbuffer.scale(2.0 * a_exp, {13875, 14155});

                pbuffer.scale(2.0 * a_exp, {16451, 16559});

                pbuffer.scale(2.0 * a_exp, {16559, 16775});

                pbuffer.scale(2.0 * a_exp, {16775, 17135});

                t2cfunc::reduce(cbuffer, 3404, pbuffer, 3741, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3434, pbuffer, 3891, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3494, pbuffer, 4191, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3594, pbuffer, 6536, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3639, pbuffer, 6716, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3729, pbuffer, 7076, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3879, pbuffer, 9878, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3942, pbuffer, 10067, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4068, pbuffer, 10445, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4278, pbuffer, 13371, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4362, pbuffer, 13539, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4530, pbuffer, 13875, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4810, pbuffer, 16451, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4918, pbuffer, 16559, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5134, pbuffer, 16775, 360, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3741, 3771});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3891, 3951});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {4191, 4291});

                pbuffer.scale(pfactors, 0, 2.0, {4691, 4841});

                pbuffer.scale(pfactors, 0, 2.0, {5441, 5651});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {6536, 6581});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {6716, 6806});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {7076, 7226});

                pbuffer.scale(pfactors, 0, 2.0, {7676, 7901});

                pbuffer.scale(pfactors, 0, 2.0, {8576, 8891});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {9878, 9941});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {10067, 10193});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {10445, 10655});

                pbuffer.scale(pfactors, 0, 2.0, {11075, 11390});

                pbuffer.scale(pfactors, 0, 2.0, {12020, 12461});

                t2cfunc::reduce(cbuffer, 874, pbuffer, 3741, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 904, pbuffer, 3891, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 964, pbuffer, 4191, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1064, pbuffer, 4691, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1214, pbuffer, 5441, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1424, pbuffer, 6536, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1469, pbuffer, 6716, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1559, pbuffer, 7076, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1709, pbuffer, 7676, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1934, pbuffer, 8576, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2249, pbuffer, 9878, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2312, pbuffer, 10067, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2438, pbuffer, 10445, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2648, pbuffer, 11075, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2963, pbuffer, 12020, 441, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {3741, 3771});

                pbuffer.scale(2.0 * a_exp, {3891, 3951});

                pbuffer.scale(2.0 * a_exp, {4191, 4291});

                pbuffer.scale(2.0 * a_exp, {4691, 4841});

                pbuffer.scale(2.0 * a_exp, {5441, 5651});

                pbuffer.scale(2.0 * a_exp, {6536, 6581});

                pbuffer.scale(2.0 * a_exp, {6716, 6806});

                pbuffer.scale(2.0 * a_exp, {7076, 7226});

                pbuffer.scale(2.0 * a_exp, {7676, 7901});

                pbuffer.scale(2.0 * a_exp, {8576, 8891});

                pbuffer.scale(2.0 * a_exp, {9878, 9941});

                pbuffer.scale(2.0 * a_exp, {10067, 10193});

                pbuffer.scale(2.0 * a_exp, {10445, 10655});

                pbuffer.scale(2.0 * a_exp, {11075, 11390});

                pbuffer.scale(2.0 * a_exp, {12020, 12461});

                pbuffer.scale(pfactors, 0, 2.0, {13371, 13455});

                pbuffer.scale(pfactors, 0, 2.0, {13539, 13707});

                pbuffer.scale(pfactors, 0, 2.0, {13875, 14155});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {14435, 14855});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {15275, 15863});

                pbuffer.scale(pfactors, 0, 2.0, {16451, 16559});

                pbuffer.scale(pfactors, 0, 2.0, {16559, 16775});

                pbuffer.scale(pfactors, 0, 2.0, {16775, 17135});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {17135, 17675});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {17675, 18431});

                t2cfunc::reduce(cbuffer, 5494, pbuffer, 3741, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5524, pbuffer, 3891, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5584, pbuffer, 4191, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5684, pbuffer, 4691, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5834, pbuffer, 5441, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6044, pbuffer, 6536, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6089, pbuffer, 6716, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6179, pbuffer, 7076, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6329, pbuffer, 7676, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6554, pbuffer, 8576, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6869, pbuffer, 9878, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6932, pbuffer, 10067, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7058, pbuffer, 10445, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7268, pbuffer, 11075, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7583, pbuffer, 12020, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8024, pbuffer, 13371, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8108, pbuffer, 13539, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8276, pbuffer, 13875, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8556, pbuffer, 14435, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8976, pbuffer, 15275, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9564, pbuffer, 16451, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9672, pbuffer, 16559, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9888, pbuffer, 16775, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10248, pbuffer, 17135, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10788, pbuffer, 17675, 756, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1020, cbuffer, 0, 30, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1380, cbuffer, 30, 90, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 3000, 1020, 1380, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 7230, cbuffer, 190, 235, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 7770, cbuffer, 235, 325, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 10200, 7230, 7770, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 16392, cbuffer, 475, 538, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 17148, cbuffer, 538, 664, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 20550, 16392, 17148, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 27240, cbuffer, 3404, 3434, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 27600, cbuffer, 3434, 3494, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 29220, 27240, 27600, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 33450, cbuffer, 3594, 3639, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 33990, cbuffer, 3639, 3729, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 36420, 33450, 33990, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 42612, cbuffer, 3879, 3942, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 43368, cbuffer, 3942, 4068, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 46770, 42612, 43368, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 55296, cbuffer, 4278, 4362, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 56304, cbuffer, 4362, 4530, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 60840, 55296, 56304, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 72072, cbuffer, 4810, 4918, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 73368, cbuffer, 4918, 5134, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 79200, 72072, 73368, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 874, 904, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 90, cbuffer, 904, 964, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 270, cbuffer, 964, 1064, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 570, cbuffer, 1064, 1214, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1110, cbuffer, 0, 0, 90, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1560, cbuffer, 30, 90, 270, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2100, cbuffer, 90, 270, 570, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 3180, 1020, 1110, 1560, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 3720, 1380, 1560, 2100, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 4800, 3000, 3180, 3720, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5700, cbuffer, 1424, 1469, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5835, cbuffer, 1469, 1559, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6105, cbuffer, 1559, 1709, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 6555, cbuffer, 1709, 1934, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 7365, cbuffer, 190, 5700, 5835, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 8040, cbuffer, 235, 5835, 6105, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 8850, cbuffer, 325, 6105, 6555, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 10470, 7230, 7365, 8040, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 11280, 7770, 8040, 8850, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 12900, 10200, 10470, 11280, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 14250, cbuffer, 2249, 2312, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 14439, cbuffer, 2312, 2438, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 14817, cbuffer, 2438, 2648, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 15447, cbuffer, 2648, 2963, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 16581, cbuffer, 475, 14250, 14439, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 17526, cbuffer, 538, 14439, 14817, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 18660, cbuffer, 664, 14817, 15447, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 20928, 16392, 16581, 17526, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 22062, 17148, 17526, 18660, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 24330, 20550, 20928, 22062, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 26220, cbuffer, 5494, 5524, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 26310, cbuffer, 5524, 5584, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 26490, cbuffer, 5584, 5684, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 26790, cbuffer, 5684, 5834, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 27330, cbuffer, 3404, 26220, 26310, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 27780, cbuffer, 3434, 26310, 26490, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 28320, cbuffer, 3494, 26490, 26790, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 29400, 27240, 27330, 27780, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 29940, 27600, 27780, 28320, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 31020, 29220, 29400, 29940, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 31920, cbuffer, 6044, 6089, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 32055, cbuffer, 6089, 6179, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 32325, cbuffer, 6179, 6329, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 32775, cbuffer, 6329, 6554, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 33585, cbuffer, 3594, 31920, 32055, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 34260, cbuffer, 3639, 32055, 32325, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 35070, cbuffer, 3729, 32325, 32775, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 36690, 33450, 33585, 34260, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 37500, 33990, 34260, 35070, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 39120, 36420, 36690, 37500, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 40470, cbuffer, 6869, 6932, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 40659, cbuffer, 6932, 7058, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 41037, cbuffer, 7058, 7268, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 41667, cbuffer, 7268, 7583, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 42801, cbuffer, 3879, 40470, 40659, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 43746, cbuffer, 3942, 40659, 41037, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 44880, cbuffer, 4068, 41037, 41667, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 47148, 42612, 42801, 43746, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 48282, 43368, 43746, 44880, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 50550, 46770, 47148, 48282, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 52440, cbuffer, 8024, 8108, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 52692, cbuffer, 8108, 8276, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 53196, cbuffer, 8276, 8556, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 54036, cbuffer, 8556, 8976, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 55548, cbuffer, 4278, 52440, 52692, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 56808, cbuffer, 4362, 52692, 53196, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 58320, cbuffer, 4530, 53196, 54036, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 61344, 55296, 55548, 56808, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 62856, 56304, 56808, 58320, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 65880, 60840, 61344, 62856, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 68400, cbuffer, 9564, 9672, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 68724, cbuffer, 9672, 9888, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 69372, cbuffer, 9888, 10248, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 70452, cbuffer, 10248, 10788, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 72396, cbuffer, 4810, 68400, 68724, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 74016, cbuffer, 4918, 68724, 69372, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 75960, cbuffer, 5134, 69372, 70452, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 79848, 72072, 72396, 74016, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 81792, 73368, 74016, 75960, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 85680, 79200, 79848, 81792, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<3, 1>(skbuffer, 0, ckbuffer, 4800, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 210, ckbuffer, 5100, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 420, ckbuffer, 5400, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 2520, ckbuffer, 12900, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 2835, ckbuffer, 13350, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 3150, ckbuffer, 13800, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 6300, ckbuffer, 24330, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 6741, ckbuffer, 24960, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 7182, ckbuffer, 25590, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 98721, ckbuffer, 31020, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 98931, ckbuffer, 31320, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 99141, ckbuffer, 31620, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 99351, ckbuffer, 39120, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 99666, ckbuffer, 39570, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 99981, ckbuffer, 40020, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 100296, ckbuffer, 50550, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 100737, ckbuffer, 51180, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 101178, ckbuffer, 51810, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 101619, ckbuffer, 65880, 0, 6);

            t4cfunc::ket_transform<3, 1>(skbuffer, 102207, ckbuffer, 66720, 0, 6);

            t4cfunc::ket_transform<3, 1>(skbuffer, 102795, ckbuffer, 67560, 0, 6);

            t4cfunc::ket_transform<3, 1>(skbuffer, 103383, ckbuffer, 85680, 0, 7);

            t4cfunc::ket_transform<3, 1>(skbuffer, 104139, ckbuffer, 86760, 0, 7);

            t4cfunc::ket_transform<3, 1>(skbuffer, 104895, ckbuffer, 87840, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 16884, 0, 2520, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 17514, 210, 2835, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 18144, 420, 3150, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 24444, 2520, 6300, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 25389, 2835, 6741, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 26334, 3150, 7182, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 47691, 16884, 24444, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 48951, 17514, 25389, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 50211, 18144, 26334, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 630, 98721, 99351, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 3465, 99351, 100296, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 7623, 100296, 101619, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sixx(skbuffer, 11592, 101619, 103383, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 18774, 0, 630, 3465, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 27279, 2520, 3465, 7623, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_phxx(skbuffer, 35784, 6300, 7623, 11592, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 51471, 16884, 18774, 27279, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dgxx(skbuffer, 62811, 24444, 27279, 35784, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ffxx(skbuffer, 79821, 47691, 51471, 62811, r_ab, 3, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 79821, 3, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1029, skbuffer, 81921, 3, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2058, skbuffer, 84021, 3, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 3087, skbuffer, 86121, 3, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 4116, skbuffer, 88221, 3, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 5145, skbuffer, 90321, 3, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 6174, skbuffer, 92421, 3, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 7203, skbuffer, 94521, 3, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 8232, skbuffer, 96621, 3, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 3, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFFFP_hpp */
