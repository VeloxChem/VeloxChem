#ifndef ElectronRepulsionGeom1100RecDPPF_hpp
#define ElectronRepulsionGeom1100RecDPPF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DP|1/|r-r'||PF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dppf(T& distributor,
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

    CSimdArray<double> pbuffer(3740, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2300, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(6630, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(18438, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2835, 1);

    // setup Boys fuction data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_data(11, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 10, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 10);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 10, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 25, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 28, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 31, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 34, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 37, 0, 1, 10, 13, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 43, 1, 2, 13, 16, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 49, 2, 3, 16, 19, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 55, 3, 4, 19, 22, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 61, 4, 5, 22, 25, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 67, 5, 6, 25, 28, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 73, 6, 7, 28, 31, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 79, 7, 8, 31, 34, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 10, 13, 37, 43, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 13, 16, 43, 49, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 16, 19, 49, 55, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 19, 22, 55, 61, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 22, 25, 61, 67, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 25, 28, 67, 73, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 28, 31, 73, 79, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 155, 37, 43, 85, 95, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 170, 43, 49, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 185, 49, 55, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 200, 55, 61, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 215, 61, 67, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 230, 67, 73, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 245, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 248, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 251, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 260, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 269, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 278, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 296, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 314, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 332, 25, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 350, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 380, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 410, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 440, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 470, 67, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 500, 95, 155, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 545, 105, 170, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 590, 115, 185, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 635, 125, 200, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 680, 135, 215, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 725, 3, 4, 245, 248, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 731, 16, 19, 245, 251, 260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 749, 19, 22, 248, 260, 269, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 767, 43, 49, 251, 278, 296, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 803, 49, 55, 260, 296, 314, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 839, 55, 61, 269, 314, 332, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 875, 85, 95, 278, 350, 380, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 935, 95, 105, 296, 380, 410, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 995, 105, 115, 314, 410, 440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1055, 115, 125, 332, 440, 470, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1115, 155, 170, 380, 500, 545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1205, 170, 185, 410, 545, 590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1295, 185, 200, 440, 590, 635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1385, 200, 215, 470, 635, 680, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1475, 251, 260, 725, 731, 749, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1505, 278, 296, 731, 767, 803, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1565, 296, 314, 749, 803, 839, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1625, 350, 380, 767, 875, 935, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1725, 380, 410, 803, 935, 995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1825, 410, 440, 839, 995, 1055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1925, 500, 545, 935, 1115, 1205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2075, 545, 590, 995, 1205, 1295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2225, 590, 635, 1055, 1295, 1385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2375, 767, 803, 1475, 1505, 1565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2465, 875, 935, 1505, 1625, 1725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2615, 935, 995, 1565, 1725, 1825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 2765, 1115, 1205, 1725, 1925, 2075, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 2990, 1205, 1295, 1825, 2075, 2225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3215, 1625, 1725, 2375, 2465, 2615, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 3425, 1925, 2075, 2615, 2765, 2990, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 25, pbuffer, 350, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 55, pbuffer, 500, 45, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {875, 935});

                pbuffer.scale(2.0 * b_exp, {1115, 1205});

                pbuffer.scale(2.0 * b_exp, {1625, 1725});

                pbuffer.scale(2.0 * b_exp, {1925, 2075});

                t2cfunc::reduce(cbuffer, 100, pbuffer, 875, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 1115, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 250, pbuffer, 1625, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 350, pbuffer, 1925, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {85, 95});

                pbuffer.scale(2.0 * a_exp, {155, 170});

                pbuffer.scale(2.0 * a_exp, {350, 380});

                pbuffer.scale(2.0 * a_exp, {500, 545});

                pbuffer.scale(a_exp / b_exp, {875, 935});

                pbuffer.scale(a_exp / b_exp, {1115, 1205});

                pbuffer.scale(a_exp / b_exp, {1625, 1725});

                pbuffer.scale(a_exp / b_exp, {1925, 2075});

                t2cfunc::reduce(cbuffer, 500, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 510, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 525, pbuffer, 350, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 555, pbuffer, 500, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 600, pbuffer, 875, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 660, pbuffer, 1115, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 750, pbuffer, 1625, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 850, pbuffer, 1925, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {875, 935});

                pbuffer.scale(2.0 * b_exp, {1115, 1205});

                pbuffer.scale(2.0 * b_exp, {1625, 1725});

                pbuffer.scale(2.0 * b_exp, {1925, 2075});

                pbuffer.scale(4.0 * a_exp * b_exp, {2465, 2615});

                pbuffer.scale(4.0 * a_exp * b_exp, {2765, 2990});

                pbuffer.scale(4.0 * a_exp * b_exp, {3215, 3425});

                pbuffer.scale(4.0 * a_exp * b_exp, {3425, 3740});

                t2cfunc::reduce(cbuffer, 1000, pbuffer, 875, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1060, pbuffer, 1115, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1150, pbuffer, 1625, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1250, pbuffer, 1925, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1400, pbuffer, 2465, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1550, pbuffer, 2765, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1775, pbuffer, 3215, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1985, pbuffer, 3425, 315, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 10, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 30, cbuffer, 25, 55, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 930, cbuffer, 100, 160, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1110, cbuffer, 250, 350, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1410, cbuffer, 500, 510, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1440, cbuffer, 525, 555, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1800, cbuffer, 600, 660, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2520, cbuffer, 750, 850, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5070, cbuffer, 1000, 1060, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5250, cbuffer, 1150, 1250, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5550, cbuffer, 1400, 1550, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 6000, cbuffer, 1775, 1985, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<1, 3>(skbuffer, 0, ckbuffer, 0, 0, 0);

            t4cfunc::ket_transform<1, 3>(skbuffer, 21, ckbuffer, 30, 0, 1);

            t4cfunc::ket_transform<1, 3>(skbuffer, 13881, ckbuffer, 930, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 14007, ckbuffer, 1110, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 14217, ckbuffer, 1410, 0, 0);

            t4cfunc::ket_transform<1, 3>(skbuffer, 14238, ckbuffer, 1440, 0, 1);

            t4cfunc::ket_transform<1, 3>(skbuffer, 14490, ckbuffer, 1800, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 14994, ckbuffer, 2520, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 17346, ckbuffer, 5070, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 17472, ckbuffer, 5250, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 17682, ckbuffer, 5550, 0, 4);

            t4cfunc::ket_transform<1, 3>(skbuffer, 17997, ckbuffer, 6000, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 16779, 14238, 14490, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 16968, 14490, 14994, r_ab, 1, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 4809, 21, 16779, 16968, r_ab, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 84, 0, 13881, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 840, 21, 14007, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 4242, 21, 84, 840, r_ab, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 14301, 14217, 17346, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 14616, 14238, 17472, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 15204, 14490, 17682, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 15834, 14994, 17997, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 273, 14238, 14301, 14616, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 1218, 14490, 14616, 15204, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 2352, 14994, 15204, 15834, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 5376, 84, 16779, 273, 1218, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 7077, 840, 16968, 1218, 2352, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 10479, 4242, 4809, 5376, 7077, r_ab, 1, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 10479, 1, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 315, skbuffer, 10857, 1, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 630, skbuffer, 11235, 1, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 945, skbuffer, 11613, 1, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1260, skbuffer, 11991, 1, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1575, skbuffer, 12369, 1, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1890, skbuffer, 12747, 1, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2205, skbuffer, 13125, 1, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2520, skbuffer, 13503, 1, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 1, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDPPF_hpp */
