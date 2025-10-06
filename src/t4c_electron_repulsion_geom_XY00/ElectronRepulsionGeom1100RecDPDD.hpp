#ifndef ElectronRepulsionGeom1100RecDPDD_hpp
#define ElectronRepulsionGeom1100RecDPDD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecXXDD.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
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
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DP|1/|r-r'||DD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dpdd(T& distributor,
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

    CSimdArray<double> pbuffer(4191, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2852, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(12372, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(21950, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3375, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 245, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 248, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 251, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 254, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 263, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 272, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 281, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 290, 13, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 308, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 326, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 344, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 362, 25, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 380, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 410, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 440, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 470, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 500, 67, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 530, 95, 155, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 575, 105, 170, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 620, 115, 185, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 665, 125, 200, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 710, 135, 215, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 755, 2, 3, 245, 248, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 761, 3, 4, 248, 251, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 767, 13, 16, 245, 254, 263, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 785, 16, 19, 248, 263, 272, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 803, 19, 22, 251, 272, 281, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 821, 37, 43, 254, 290, 308, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 857, 43, 49, 263, 308, 326, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 893, 49, 55, 272, 326, 344, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 929, 55, 61, 281, 344, 362, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 965, 85, 95, 308, 380, 410, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1025, 95, 105, 326, 410, 440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1085, 105, 115, 344, 440, 470, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1145, 115, 125, 362, 470, 500, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1205, 155, 170, 410, 530, 575, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1295, 170, 185, 440, 575, 620, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1385, 185, 200, 470, 620, 665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1475, 200, 215, 500, 665, 710, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1565, 245, 248, 755, 761, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1575, 254, 263, 755, 767, 785, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1605, 263, 272, 761, 785, 803, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1635, 290, 308, 767, 821, 857, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1695, 308, 326, 785, 857, 893, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1755, 326, 344, 803, 893, 929, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1815, 380, 410, 857, 965, 1025, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1915, 410, 440, 893, 1025, 1085, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2015, 440, 470, 929, 1085, 1145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2115, 530, 575, 1025, 1205, 1295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2265, 575, 620, 1085, 1295, 1385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2415, 620, 665, 1145, 1385, 1475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2565, 767, 785, 1565, 1575, 1605, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2610, 821, 857, 1575, 1635, 1695, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2700, 857, 893, 1605, 1695, 1755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2790, 965, 1025, 1695, 1815, 1915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2940, 1025, 1085, 1755, 1915, 2015, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3090, 1205, 1295, 1915, 2115, 2265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3315, 1295, 1385, 2015, 2265, 2415, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3540, 1635, 1695, 2565, 2610, 2700, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3666, 1815, 1915, 2700, 2790, 2940, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 3876, 2115, 2265, 2940, 3090, 3315, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 31, pbuffer, 290, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 49, pbuffer, 380, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 79, pbuffer, 530, 45, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {821, 857});

                pbuffer.scale(2.0 * b_exp, {965, 1025});

                pbuffer.scale(2.0 * b_exp, {1205, 1295});

                pbuffer.scale(2.0 * b_exp, {1635, 1695});

                pbuffer.scale(2.0 * b_exp, {1815, 1915});

                pbuffer.scale(2.0 * b_exp, {2115, 2265});

                t2cfunc::reduce(cbuffer, 124, pbuffer, 821, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 965, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 220, pbuffer, 1205, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 310, pbuffer, 1635, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 370, pbuffer, 1815, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 470, pbuffer, 2115, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {37, 43});

                pbuffer.scale(2.0 * a_exp, {85, 95});

                pbuffer.scale(2.0 * a_exp, {155, 170});

                pbuffer.scale(2.0 * a_exp, {290, 308});

                pbuffer.scale(2.0 * a_exp, {380, 410});

                pbuffer.scale(2.0 * a_exp, {530, 575});

                pbuffer.scale(a_exp / b_exp, {821, 857});

                pbuffer.scale(a_exp / b_exp, {965, 1025});

                pbuffer.scale(a_exp / b_exp, {1205, 1295});

                pbuffer.scale(a_exp / b_exp, {1635, 1695});

                pbuffer.scale(a_exp / b_exp, {1815, 1915});

                pbuffer.scale(a_exp / b_exp, {2115, 2265});

                t2cfunc::reduce(cbuffer, 620, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 626, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 636, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 651, pbuffer, 290, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 669, pbuffer, 380, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 699, pbuffer, 530, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 744, pbuffer, 821, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 780, pbuffer, 965, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 840, pbuffer, 1205, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 930, pbuffer, 1635, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 990, pbuffer, 1815, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1090, pbuffer, 2115, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {821, 857});

                pbuffer.scale(2.0 * b_exp, {965, 1025});

                pbuffer.scale(2.0 * b_exp, {1205, 1295});

                pbuffer.scale(2.0 * b_exp, {1635, 1695});

                pbuffer.scale(2.0 * b_exp, {1815, 1915});

                pbuffer.scale(2.0 * b_exp, {2115, 2265});

                pbuffer.scale(4.0 * a_exp * b_exp, {2610, 2700});

                pbuffer.scale(4.0 * a_exp * b_exp, {2790, 2940});

                pbuffer.scale(4.0 * a_exp * b_exp, {3090, 3315});

                pbuffer.scale(4.0 * a_exp * b_exp, {3540, 3666});

                pbuffer.scale(4.0 * a_exp * b_exp, {3666, 3876});

                pbuffer.scale(4.0 * a_exp * b_exp, {3876, 4191});

                t2cfunc::reduce(cbuffer, 1240, pbuffer, 821, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1276, pbuffer, 965, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1336, pbuffer, 1205, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1426, pbuffer, 1635, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1486, pbuffer, 1815, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1586, pbuffer, 2115, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1736, pbuffer, 2610, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1826, pbuffer, 2790, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1976, pbuffer, 3090, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2201, pbuffer, 3540, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2327, pbuffer, 3666, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2537, pbuffer, 3876, 315, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 6, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 18, cbuffer, 6, 16, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 48, 0, 18, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 84, cbuffer, 31, 49, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 138, cbuffer, 49, 79, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 228, 84, 138, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1308, cbuffer, 124, 160, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1416, cbuffer, 160, 220, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 1596, 1308, 1416, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1812, cbuffer, 310, 370, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1992, cbuffer, 370, 470, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 2292, 1812, 1992, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2652, cbuffer, 620, 626, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2670, cbuffer, 626, 636, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 2700, 2652, 2670, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2736, cbuffer, 651, 669, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2790, cbuffer, 669, 699, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 2880, 2736, 2790, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3312, cbuffer, 744, 780, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 3420, cbuffer, 780, 840, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 3600, 3312, 3420, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4464, cbuffer, 930, 990, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4644, cbuffer, 990, 1090, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 4944, 4464, 4644, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8004, cbuffer, 1240, 1276, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 8112, cbuffer, 1276, 1336, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 8292, 8004, 8112, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8508, cbuffer, 1426, 1486, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 8688, cbuffer, 1486, 1586, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 8988, 8508, 8688, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 9348, cbuffer, 1736, 1826, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9618, cbuffer, 1826, 1976, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 10068, 9348, 9618, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 10608, cbuffer, 2201, 2327, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 10986, cbuffer, 2327, 2537, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 11616, 10608, 10986, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 0, ckbuffer, 48, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 25, ckbuffer, 228, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 16525, ckbuffer, 1596, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 16675, ckbuffer, 2292, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 16925, ckbuffer, 2700, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 16950, ckbuffer, 2880, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 17250, ckbuffer, 3600, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 17850, ckbuffer, 4944, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 20650, ckbuffer, 8292, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 20800, ckbuffer, 8988, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21050, ckbuffer, 10068, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21425, ckbuffer, 11616, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 19975, 16950, 17250, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 20200, 17250, 17850, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 5725, 25, 19975, 20200, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 100, 0, 16525, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 1000, 25, 16675, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 5050, 25, 100, 1000, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 17025, 16925, 20650, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 17400, 16950, 20800, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 18100, 17250, 21050, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 18850, 17850, 21425, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 325, 16950, 17025, 17400, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 1450, 17250, 17400, 18100, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 2800, 17850, 18100, 18850, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 6400, 100, 19975, 325, 1450, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 8425, 1000, 20200, 1450, 2800, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 12475, 5050, 5725, 6400, 8425, r_ab, 2, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 12475, 2, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 375, skbuffer, 12925, 2, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 750, skbuffer, 13375, 2, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1125, skbuffer, 13825, 2, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1500, skbuffer, 14275, 2, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1875, skbuffer, 14725, 2, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2250, skbuffer, 15175, 2, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2625, skbuffer, 15625, 2, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 3000, skbuffer, 16075, 2, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 2, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDPDD_hpp */
