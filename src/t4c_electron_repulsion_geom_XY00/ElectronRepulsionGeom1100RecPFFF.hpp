#ifndef ElectronRepulsionGeom1100RecPFFF_hpp
#define ElectronRepulsionGeom1100RecPFFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecXXDF.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXFF.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
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
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSI.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(PF|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_pfff(T& distributor,
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

    CSimdArray<double> pbuffer(16682, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(8584, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(61808, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(39641, 1);

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

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 518, 215, 230, 350, 371, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 546, 230, 245, 371, 392, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 574, 245, 260, 392, 413, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 602, 260, 275, 413, 434, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 630, 275, 290, 434, 455, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 658, 290, 305, 455, 476, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 686, 305, 320, 476, 497, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 714, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 717, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 720, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 723, 3, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 732, 4, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 741, 5, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 750, 6, 28, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 759, 19, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 777, 22, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 795, 25, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 813, 28, 73, 79, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 831, 31, 79, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 849, 55, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 879, 61, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 909, 67, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 939, 73, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 969, 79, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 999, 85, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1029, 125, 215, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1074, 135, 230, 245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1119, 145, 245, 260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1164, 155, 260, 275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1209, 165, 275, 290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1254, 175, 290, 305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1299, 230, 350, 371, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1362, 245, 371, 392, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1425, 260, 392, 413, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1488, 275, 413, 434, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1551, 290, 434, 455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1614, 305, 455, 476, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1677, 371, 518, 546, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1761, 392, 546, 574, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1845, 413, 574, 602, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1929, 434, 602, 630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2013, 455, 630, 658, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 2097, 476, 658, 686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2181, 3, 4, 714, 717, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 2187, 4, 5, 717, 720, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2193, 19, 22, 714, 723, 732, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2211, 22, 25, 717, 732, 741, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2229, 25, 28, 720, 741, 750, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2247, 55, 61, 723, 759, 777, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2283, 61, 67, 732, 777, 795, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2319, 67, 73, 741, 795, 813, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2355, 73, 79, 750, 813, 831, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2391, 115, 125, 759, 849, 879, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2451, 125, 135, 777, 879, 909, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2511, 135, 145, 795, 909, 939, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2571, 145, 155, 813, 939, 969, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2631, 155, 165, 831, 969, 999, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2691, 215, 230, 879, 1029, 1074, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2781, 230, 245, 909, 1074, 1119, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2871, 245, 260, 939, 1119, 1164, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2961, 260, 275, 969, 1164, 1209, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 3051, 275, 290, 999, 1209, 1254, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3141, 350, 371, 1074, 1299, 1362, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3267, 371, 392, 1119, 1362, 1425, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3393, 392, 413, 1164, 1425, 1488, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3519, 413, 434, 1209, 1488, 1551, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3645, 434, 455, 1254, 1551, 1614, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3771, 518, 546, 1362, 1677, 1761, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3939, 546, 574, 1425, 1761, 1845, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4107, 574, 602, 1488, 1845, 1929, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4275, 602, 630, 1551, 1929, 2013, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 4443, 630, 658, 1614, 2013, 2097, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 4611, 714, 717, 2181, 2187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 4621, 723, 732, 2181, 2193, 2211, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 4651, 732, 741, 2187, 2211, 2229, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4681, 759, 777, 2193, 2247, 2283, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4741, 777, 795, 2211, 2283, 2319, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4801, 795, 813, 2229, 2319, 2355, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4861, 849, 879, 2247, 2391, 2451, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4961, 879, 909, 2283, 2451, 2511, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 5061, 909, 939, 2319, 2511, 2571, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 5161, 939, 969, 2355, 2571, 2631, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5261, 1029, 1074, 2451, 2691, 2781, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5411, 1074, 1119, 2511, 2781, 2871, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5561, 1119, 1164, 2571, 2871, 2961, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 5711, 1164, 1209, 2631, 2961, 3051, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5861, 1299, 1362, 2781, 3141, 3267, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 6071, 1362, 1425, 2871, 3267, 3393, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 6281, 1425, 1488, 2961, 3393, 3519, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 6491, 1488, 1551, 3051, 3519, 3645, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 6701, 1677, 1761, 3267, 3771, 3939, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 6981, 1761, 1845, 3393, 3939, 4107, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 7261, 1845, 1929, 3519, 4107, 4275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 7541, 1929, 2013, 3645, 4275, 4443, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 7821, 2193, 2211, 4611, 4621, 4651, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 7866, 2247, 2283, 4621, 4681, 4741, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 7956, 2283, 2319, 4651, 4741, 4801, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 8046, 2391, 2451, 4681, 4861, 4961, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 8196, 2451, 2511, 4741, 4961, 5061, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 8346, 2511, 2571, 4801, 5061, 5161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 8496, 2691, 2781, 4961, 5261, 5411, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 8721, 2781, 2871, 5061, 5411, 5561, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 8946, 2871, 2961, 5161, 5561, 5711, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 9171, 3141, 3267, 5411, 5861, 6071, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 9486, 3267, 3393, 5561, 6071, 6281, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 9801, 3393, 3519, 5711, 6281, 6491, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 10116, 3771, 3939, 6071, 6701, 6981, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 10536, 3939, 4107, 6281, 6981, 7261, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 10956, 4107, 4275, 6491, 7261, 7541, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 11376, 4681, 4741, 7821, 7866, 7956, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 11502, 4861, 4961, 7866, 8046, 8196, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 11712, 4961, 5061, 7956, 8196, 8346, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 11922, 5261, 5411, 8196, 8496, 8721, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 12237, 5411, 5561, 8346, 8721, 8946, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 12552, 5861, 6071, 8721, 9171, 9486, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 12993, 6071, 6281, 8946, 9486, 9801, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 13434, 6701, 6981, 9486, 10116, 10536, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 14022, 6981, 7261, 9801, 10536, 10956, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 14610, 8046, 8196, 11376, 11502, 11712, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 14890, 8496, 8721, 11712, 11922, 12237, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 15310, 9171, 9486, 12237, 12552, 12993, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisi(pbuffer, 15898, 10116, 10536, 12993, 13434, 14022, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 2391, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 2691, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 150, pbuffer, 3141, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 276, pbuffer, 3771, 168, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {8046, 8196});

                pbuffer.scale(2.0 * b_exp, {8496, 8721});

                pbuffer.scale(2.0 * b_exp, {9171, 9486});

                pbuffer.scale(2.0 * b_exp, {10116, 10536});

                t2cfunc::reduce(cbuffer, 444, pbuffer, 8046, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 594, pbuffer, 8496, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 819, pbuffer, 9171, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1134, pbuffer, 10116, 420, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {2391, 2451});

                pbuffer.scale(2.0 * a_exp, {2691, 2781});

                pbuffer.scale(2.0 * a_exp, {3141, 3267});

                pbuffer.scale(2.0 * a_exp, {3771, 3939});

                pbuffer.scale(2.0 * a_exp, {4861, 4961});

                pbuffer.scale(2.0 * a_exp, {5261, 5411});

                pbuffer.scale(2.0 * a_exp, {5861, 6071});

                pbuffer.scale(2.0 * a_exp, {6701, 6981});

                pbuffer.scale(a_exp / b_exp, {8046, 8196});

                pbuffer.scale(a_exp / b_exp, {8496, 8721});

                pbuffer.scale(a_exp / b_exp, {9171, 9486});

                pbuffer.scale(a_exp / b_exp, {10116, 10536});

                t2cfunc::reduce(cbuffer, 1554, pbuffer, 2391, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1614, pbuffer, 2691, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1704, pbuffer, 3141, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1830, pbuffer, 3771, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1998, pbuffer, 4861, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2098, pbuffer, 5261, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2248, pbuffer, 5861, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2458, pbuffer, 6701, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2738, pbuffer, 8046, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2888, pbuffer, 8496, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3113, pbuffer, 9171, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3428, pbuffer, 10116, 420, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {8046, 8196});

                pbuffer.scale(2.0 * b_exp, {8496, 8721});

                pbuffer.scale(2.0 * b_exp, {9171, 9486});

                pbuffer.scale(2.0 * b_exp, {10116, 10536});

                pbuffer.scale(4.0 * a_exp * b_exp, {11502, 11712});

                pbuffer.scale(4.0 * a_exp * b_exp, {11922, 12237});

                pbuffer.scale(4.0 * a_exp * b_exp, {12552, 12993});

                pbuffer.scale(4.0 * a_exp * b_exp, {13434, 14022});

                pbuffer.scale(4.0 * a_exp * b_exp, {14610, 14890});

                pbuffer.scale(4.0 * a_exp * b_exp, {14890, 15310});

                pbuffer.scale(4.0 * a_exp * b_exp, {15310, 15898});

                pbuffer.scale(4.0 * a_exp * b_exp, {15898, 16682});

                t2cfunc::reduce(cbuffer, 3848, pbuffer, 8046, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3998, pbuffer, 8496, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4223, pbuffer, 9171, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4538, pbuffer, 10116, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4958, pbuffer, 11502, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5168, pbuffer, 11922, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5483, pbuffer, 12552, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5924, pbuffer, 13434, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6512, pbuffer, 14610, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6792, pbuffer, 14890, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7212, pbuffer, 15310, 588, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7800, pbuffer, 15898, 784, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 60, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 180, cbuffer, 60, 150, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 450, cbuffer, 150, 276, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 828, 0, 180, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 1188, 180, 450, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 1728, 828, 1188, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5328, cbuffer, 444, 594, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 5778, cbuffer, 594, 819, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 6453, cbuffer, 819, 1134, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 7398, 5328, 5778, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 8298, 5778, 6453, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 9648, 7398, 8298, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 11148, cbuffer, 1554, 1614, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 11328, cbuffer, 1614, 1704, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 11598, cbuffer, 1704, 1830, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 11976, 11148, 11328, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 12336, 11328, 11598, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 12876, 11976, 12336, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 13476, cbuffer, 1998, 2098, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 13776, cbuffer, 2098, 2248, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 14226, cbuffer, 2248, 2458, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 14856, 13476, 13776, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 15456, 13776, 14226, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 16356, 14856, 15456, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 20356, cbuffer, 2738, 2888, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 20806, cbuffer, 2888, 3113, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 21481, cbuffer, 3113, 3428, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 22426, 20356, 20806, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 23326, 20806, 21481, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 24676, 22426, 23326, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 36976, cbuffer, 3848, 3998, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 37426, cbuffer, 3998, 4223, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 38101, cbuffer, 4223, 4538, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 39046, 36976, 37426, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 39946, 37426, 38101, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 41296, 39046, 39946, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 42796, cbuffer, 4958, 5168, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 43426, cbuffer, 5168, 5483, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 44371, cbuffer, 5483, 5924, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 45694, 42796, 43426, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 46954, 43426, 44371, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 48844, 45694, 46954, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 50944, cbuffer, 6512, 6792, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 51784, cbuffer, 6792, 7212, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 53044, cbuffer, 7212, 7800, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 54808, 50944, 51784, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 56488, 51784, 53044, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 59008, 54808, 56488, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 1728, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 26019, ckbuffer, 9648, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 26754, ckbuffer, 12876, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 27048, ckbuffer, 16356, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 29008, ckbuffer, 24676, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 36505, ckbuffer, 41296, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 37240, ckbuffer, 48844, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 38269, ckbuffer, 59008, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 35035, 27048, 29008, r_ab, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 294, 0, 26019, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 27538, 26754, 36505, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 29743, 27048, 37240, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 31948, 29008, 38269, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 1764, 27048, 27538, 29743, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 6174, 29008, 29743, 31948, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 12789, 294, 35035, 1764, 6174, r_ab, 3, 3);

            t4cfunc::bra_transform<1, 3>(sbuffer, 0, skbuffer, 12789, 3, 3);

            t4cfunc::bra_transform<1, 3>(sbuffer, 1029, skbuffer, 14259, 3, 3);

            t4cfunc::bra_transform<1, 3>(sbuffer, 2058, skbuffer, 15729, 3, 3);

            t4cfunc::bra_transform<1, 3>(sbuffer, 3087, skbuffer, 17199, 3, 3);

            t4cfunc::bra_transform<1, 3>(sbuffer, 4116, skbuffer, 18669, 3, 3);

            t4cfunc::bra_transform<1, 3>(sbuffer, 5145, skbuffer, 20139, 3, 3);

            t4cfunc::bra_transform<1, 3>(sbuffer, 6174, skbuffer, 21609, 3, 3);

            t4cfunc::bra_transform<1, 3>(sbuffer, 7203, skbuffer, 23079, 3, 3);

            t4cfunc::bra_transform<1, 3>(sbuffer, 8232, skbuffer, 24549, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 1, 3, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecPFFF_hpp */
