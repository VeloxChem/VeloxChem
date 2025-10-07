#ifndef ElectronRepulsionGeom1010RecFPPF_hpp
#define ElectronRepulsionGeom1010RecFPPF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"
#include "ElectronRepulsionGeom1010ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
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
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FP|1/|r-r'||PF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fppf(T& distributor,
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

    CSimdArray<double> pbuffer(6421, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(4144, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(12210, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(40572, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3969, 1);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 11, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 11);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 29, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 32, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 35, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 38, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 41, 0, 1, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 47, 1, 2, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 2, 3, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 3, 4, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 65, 4, 5, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 71, 5, 6, 26, 29, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 77, 6, 7, 29, 32, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 83, 7, 8, 32, 35, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 89, 8, 9, 35, 38, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 11, 14, 41, 47, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 14, 17, 47, 53, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 17, 20, 53, 59, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 20, 23, 59, 65, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 23, 26, 65, 71, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 26, 29, 71, 77, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 29, 32, 77, 83, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 32, 35, 83, 89, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 175, 41, 47, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 190, 47, 53, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 205, 53, 59, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 220, 59, 65, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 235, 65, 71, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 250, 71, 77, 145, 155, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 265, 77, 83, 155, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 280, 95, 105, 175, 190, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 301, 105, 115, 190, 205, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 322, 115, 125, 205, 220, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 343, 125, 135, 220, 235, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 364, 135, 145, 235, 250, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 385, 145, 155, 250, 265, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 406, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 409, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 412, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 421, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 430, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 439, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 457, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 475, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 493, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 511, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 541, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 571, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 601, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 631, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 661, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 706, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 751, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 796, 135, 220, 235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 841, 145, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 886, 190, 280, 301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 949, 205, 301, 322, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1012, 220, 322, 343, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1075, 235, 343, 364, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1138, 250, 364, 385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1201, 3, 4, 406, 409, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1207, 17, 20, 406, 412, 421, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1225, 20, 23, 409, 421, 430, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1243, 47, 53, 412, 439, 457, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1279, 53, 59, 421, 457, 475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1315, 59, 65, 430, 475, 493, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1351, 95, 105, 439, 511, 541, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1411, 105, 115, 457, 541, 571, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1471, 115, 125, 475, 571, 601, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1531, 125, 135, 493, 601, 631, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1591, 175, 190, 541, 661, 706, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1681, 190, 205, 571, 706, 751, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1771, 205, 220, 601, 751, 796, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1861, 220, 235, 631, 796, 841, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1951, 280, 301, 706, 886, 949, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2077, 301, 322, 751, 949, 1012, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2203, 322, 343, 796, 1012, 1075, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2329, 343, 364, 841, 1075, 1138, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2455, 412, 421, 1201, 1207, 1225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2485, 439, 457, 1207, 1243, 1279, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2545, 457, 475, 1225, 1279, 1315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2605, 511, 541, 1243, 1351, 1411, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2705, 541, 571, 1279, 1411, 1471, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2805, 571, 601, 1315, 1471, 1531, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2905, 661, 706, 1411, 1591, 1681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3055, 706, 751, 1471, 1681, 1771, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3205, 751, 796, 1531, 1771, 1861, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3355, 886, 949, 1681, 1951, 2077, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3565, 949, 1012, 1771, 2077, 2203, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3775, 1012, 1075, 1861, 2203, 2329, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3985, 1243, 1279, 2455, 2485, 2545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4075, 1351, 1411, 2485, 2605, 2705, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4225, 1411, 1471, 2545, 2705, 2805, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4375, 1591, 1681, 2705, 2905, 3055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4600, 1681, 1771, 2805, 3055, 3205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 4825, 1951, 2077, 3055, 3355, 3565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 5140, 2077, 2203, 3205, 3565, 3775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5455, 2605, 2705, 3985, 4075, 4225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 5665, 2905, 3055, 4225, 4375, 4600, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 5980, 3355, 3565, 4600, 4825, 5140, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 511, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 1351, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 2605, 100, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {511, 541});

                pbuffer.scale(2.0 * a_exp, {1351, 1411});

                pbuffer.scale(2.0 * a_exp, {2605, 2705});

                pbuffer.scale(2.0 * a_exp, {4075, 4225});

                pbuffer.scale(2.0 * a_exp, {5455, 5665});

                t2cfunc::reduce(cbuffer, 1064, pbuffer, 511, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1094, pbuffer, 1351, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1154, pbuffer, 2605, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1254, pbuffer, 4075, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1404, pbuffer, 5455, 210, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {511, 541});

                pbuffer.scale(pfactors, 0, 2.0, {661, 706});

                pbuffer.scale(pfactors, 0, 2.0, {886, 949});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1351, 1411});

                pbuffer.scale(pfactors, 0, 2.0, {1591, 1681});

                pbuffer.scale(pfactors, 0, 2.0, {1951, 2077});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2605, 2705});

                pbuffer.scale(pfactors, 0, 2.0, {2905, 3055});

                pbuffer.scale(pfactors, 0, 2.0, {3355, 3565});

                t2cfunc::reduce(cbuffer, 190, pbuffer, 511, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 220, pbuffer, 661, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 265, pbuffer, 886, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 328, pbuffer, 1351, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 388, pbuffer, 1591, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 478, pbuffer, 1951, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 604, pbuffer, 2605, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 704, pbuffer, 2905, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 854, pbuffer, 3355, 210, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {511, 541});

                pbuffer.scale(2.0 * a_exp, {661, 706});

                pbuffer.scale(2.0 * a_exp, {886, 949});

                pbuffer.scale(2.0 * a_exp, {1351, 1411});

                pbuffer.scale(2.0 * a_exp, {1591, 1681});

                pbuffer.scale(2.0 * a_exp, {1951, 2077});

                pbuffer.scale(2.0 * a_exp, {2605, 2705});

                pbuffer.scale(2.0 * a_exp, {2905, 3055});

                pbuffer.scale(2.0 * a_exp, {3355, 3565});

                pbuffer.scale(pfactors, 0, 2.0, {4075, 4225});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {4375, 4600});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {4825, 5140});

                pbuffer.scale(pfactors, 0, 2.0, {5455, 5665});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5665, 5980});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5980, 6421});

                t2cfunc::reduce(cbuffer, 1614, pbuffer, 511, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1644, pbuffer, 661, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1689, pbuffer, 886, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1752, pbuffer, 1351, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1812, pbuffer, 1591, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1902, pbuffer, 1951, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2028, pbuffer, 2605, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2128, pbuffer, 2905, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2278, pbuffer, 3355, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2488, pbuffer, 4075, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2638, pbuffer, 4375, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2863, pbuffer, 4825, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3178, pbuffer, 5455, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3388, pbuffer, 5665, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3703, pbuffer, 5980, 441, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 0, cbuffer, 190, 220, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 90, cbuffer, 220, 265, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 225, cbuffer, 0, 0, 90, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 495, cbuffer, 328, 388, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 675, cbuffer, 388, 478, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 945, cbuffer, 30, 495, 675, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1485, cbuffer, 604, 704, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 1785, cbuffer, 704, 854, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2235, cbuffer, 90, 1485, 1785, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3135, cbuffer, 1614, 1644, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 3225, cbuffer, 1644, 1689, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 3360, cbuffer, 1064, 3135, 3225, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3630, cbuffer, 1752, 1812, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 3810, cbuffer, 1812, 1902, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 4080, cbuffer, 1094, 3630, 3810, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 4620, cbuffer, 2028, 2128, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 4920, cbuffer, 2128, 2278, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 5370, cbuffer, 1154, 4620, 4920, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6270, cbuffer, 2488, 2638, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 6720, cbuffer, 2638, 2863, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 7395, cbuffer, 1254, 6270, 6720, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 8745, cbuffer, 3178, 3388, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 9375, cbuffer, 3388, 3703, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 10320, cbuffer, 1404, 8745, 9375, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<1, 3>(skbuffer, 0, ckbuffer, 225, 0, 1);

            t4cfunc::ket_transform<1, 3>(skbuffer, 63, ckbuffer, 315, 0, 1);

            t4cfunc::ket_transform<1, 3>(skbuffer, 126, ckbuffer, 405, 0, 1);

            t4cfunc::ket_transform<1, 3>(skbuffer, 756, ckbuffer, 945, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 882, ckbuffer, 1125, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 1008, ckbuffer, 1305, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 2268, ckbuffer, 2235, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 2478, ckbuffer, 2535, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 2688, ckbuffer, 2835, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 37107, ckbuffer, 3360, 0, 1);

            t4cfunc::ket_transform<1, 3>(skbuffer, 37170, ckbuffer, 3450, 0, 1);

            t4cfunc::ket_transform<1, 3>(skbuffer, 37233, ckbuffer, 3540, 0, 1);

            t4cfunc::ket_transform<1, 3>(skbuffer, 37296, ckbuffer, 4080, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 37422, ckbuffer, 4260, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 37548, ckbuffer, 4440, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 37674, ckbuffer, 5370, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 37884, ckbuffer, 5670, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 38094, ckbuffer, 5970, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 38304, ckbuffer, 7395, 0, 4);

            t4cfunc::ket_transform<1, 3>(skbuffer, 38619, ckbuffer, 7845, 0, 4);

            t4cfunc::ket_transform<1, 3>(skbuffer, 38934, ckbuffer, 8295, 0, 4);

            t4cfunc::ket_transform<1, 3>(skbuffer, 39249, ckbuffer, 10320, 0, 5);

            t4cfunc::ket_transform<1, 3>(skbuffer, 39690, ckbuffer, 10950, 0, 5);

            t4cfunc::ket_transform<1, 3>(skbuffer, 40131, ckbuffer, 11580, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 7623, 0, 756, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 7812, 63, 882, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 8001, 126, 1008, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 9891, 756, 2268, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 10269, 882, 2478, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 10647, 1008, 2688, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 20097, 7623, 9891, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 20475, 7812, 10269, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 20853, 8001, 10647, r_ab, 1, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 189, 37107, 37296, r_ab, 1, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 1134, 37296, 37674, r_ab, 1, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 2898, 37674, 38304, r_ab, 1, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 4788, 38304, 39249, r_ab, 1, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 8190, 0, 189, 1134, r_ab, 1, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 11025, 756, 1134, 2898, r_ab, 1, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 14427, 2268, 2898, 4788, r_ab, 1, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 21231, 7623, 8190, 11025, r_ab, 1, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 24633, 9891, 11025, 14427, r_ab, 1, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fpxx(skbuffer, 31437, 20097, 21231, 24633, r_ab, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 31437, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 441, skbuffer, 32067, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 882, skbuffer, 32697, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1323, skbuffer, 33327, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1764, skbuffer, 33957, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2205, skbuffer, 34587, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2646, skbuffer, 35217, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 3087, skbuffer, 35847, 1, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 3528, skbuffer, 36477, 1, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 1, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFPPF_hpp */
