#ifndef ElectronRepulsionGeom1010RecFPFP_hpp
#define ElectronRepulsionGeom1010RecFPFP_hpp

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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FP|1/|r-r'||FP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fpfp(T& distributor,
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

    CSimdArray<double> pbuffer(7071, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(5476, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(42180, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 406, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 409, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 412, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 415, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 418, 1, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 427, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 436, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 445, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 454, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 463, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 481, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 499, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 517, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 535, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 553, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 583, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 613, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 643, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 673, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 703, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 748, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 793, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 838, 135, 220, 235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 883, 145, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 928, 190, 280, 301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 991, 205, 301, 322, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1054, 220, 322, 343, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1117, 235, 343, 364, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1180, 250, 364, 385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1243, 1, 2, 406, 409, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1249, 2, 3, 409, 412, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1255, 3, 4, 412, 415, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1261, 11, 14, 406, 418, 427, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1279, 14, 17, 409, 427, 436, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1297, 17, 20, 412, 436, 445, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1315, 20, 23, 415, 445, 454, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1333, 41, 47, 427, 463, 481, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1369, 47, 53, 436, 481, 499, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1405, 53, 59, 445, 499, 517, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1441, 59, 65, 454, 517, 535, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1477, 95, 105, 481, 553, 583, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1537, 105, 115, 499, 583, 613, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1597, 115, 125, 517, 613, 643, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1657, 125, 135, 535, 643, 673, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1717, 175, 190, 583, 703, 748, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1807, 190, 205, 613, 748, 793, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1897, 205, 220, 643, 793, 838, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1987, 220, 235, 673, 838, 883, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2077, 280, 301, 748, 928, 991, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2203, 301, 322, 793, 991, 1054, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2329, 322, 343, 838, 1054, 1117, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2455, 343, 364, 883, 1117, 1180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2581, 406, 409, 1243, 1249, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2591, 409, 412, 1249, 1255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2601, 418, 427, 1243, 1261, 1279, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2631, 427, 436, 1249, 1279, 1297, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2661, 436, 445, 1255, 1297, 1315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2691, 463, 481, 1279, 1333, 1369, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2751, 481, 499, 1297, 1369, 1405, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2811, 499, 517, 1315, 1405, 1441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2871, 553, 583, 1369, 1477, 1537, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2971, 583, 613, 1405, 1537, 1597, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3071, 613, 643, 1441, 1597, 1657, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3171, 703, 748, 1537, 1717, 1807, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3321, 748, 793, 1597, 1807, 1897, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3471, 793, 838, 1657, 1897, 1987, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3621, 928, 991, 1807, 2077, 2203, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3831, 991, 1054, 1897, 2203, 2329, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4041, 1054, 1117, 1987, 2329, 2455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4251, 1243, 1249, 2581, 2591, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4266, 1261, 1279, 2581, 2601, 2631, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4311, 1279, 1297, 2591, 2631, 2661, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4356, 1333, 1369, 2631, 2691, 2751, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4446, 1369, 1405, 2661, 2751, 2811, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4536, 1477, 1537, 2751, 2871, 2971, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4686, 1537, 1597, 2811, 2971, 3071, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4836, 1717, 1807, 2971, 3171, 3321, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5061, 1807, 1897, 3071, 3321, 3471, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 5286, 2077, 2203, 3321, 3621, 3831, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 5601, 2203, 2329, 3471, 3831, 4041, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 5916, 2601, 2631, 4251, 4266, 4311, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 5979, 2691, 2751, 4311, 4356, 4446, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 6105, 2871, 2971, 4446, 4536, 4686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 6315, 3171, 3321, 4686, 4836, 5061, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 6630, 3621, 3831, 5061, 5286, 5601, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 418, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 463, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 27, pbuffer, 553, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 57, pbuffer, 1261, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 75, pbuffer, 1333, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 111, pbuffer, 1477, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 171, pbuffer, 2601, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 201, pbuffer, 2691, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 261, pbuffer, 2871, 100, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {418, 427});

                pbuffer.scale(2.0 * a_exp, {463, 481});

                pbuffer.scale(2.0 * a_exp, {553, 583});

                pbuffer.scale(2.0 * a_exp, {1261, 1279});

                pbuffer.scale(2.0 * a_exp, {1333, 1369});

                pbuffer.scale(2.0 * a_exp, {1477, 1537});

                pbuffer.scale(2.0 * a_exp, {2601, 2631});

                pbuffer.scale(2.0 * a_exp, {2691, 2751});

                pbuffer.scale(2.0 * a_exp, {2871, 2971});

                pbuffer.scale(2.0 * a_exp, {4266, 4311});

                pbuffer.scale(2.0 * a_exp, {4356, 4446});

                pbuffer.scale(2.0 * a_exp, {4536, 4686});

                pbuffer.scale(2.0 * a_exp, {5916, 5979});

                pbuffer.scale(2.0 * a_exp, {5979, 6105});

                pbuffer.scale(2.0 * a_exp, {6105, 6315});

                t2cfunc::reduce(cbuffer, 1406, pbuffer, 418, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1415, pbuffer, 463, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1433, pbuffer, 553, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1463, pbuffer, 1261, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1481, pbuffer, 1333, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1517, pbuffer, 1477, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1577, pbuffer, 2601, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1607, pbuffer, 2691, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1667, pbuffer, 2871, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1767, pbuffer, 4266, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1812, pbuffer, 4356, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1902, pbuffer, 4536, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2052, pbuffer, 5916, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2115, pbuffer, 5979, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2241, pbuffer, 6105, 210, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {418, 427});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {463, 481});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {553, 583});

                pbuffer.scale(pfactors, 0, 2.0, {703, 748});

                pbuffer.scale(pfactors, 0, 2.0, {928, 991});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1261, 1279});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1333, 1369});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1477, 1537});

                pbuffer.scale(pfactors, 0, 2.0, {1717, 1807});

                pbuffer.scale(pfactors, 0, 2.0, {2077, 2203});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2601, 2631});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2691, 2751});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2871, 2971});

                pbuffer.scale(pfactors, 0, 2.0, {3171, 3321});

                pbuffer.scale(pfactors, 0, 2.0, {3621, 3831});

                t2cfunc::reduce(cbuffer, 361, pbuffer, 418, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 370, pbuffer, 463, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 388, pbuffer, 553, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 418, pbuffer, 703, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 463, pbuffer, 928, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 526, pbuffer, 1261, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 544, pbuffer, 1333, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 580, pbuffer, 1477, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 640, pbuffer, 1717, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 730, pbuffer, 2077, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 856, pbuffer, 2601, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 886, pbuffer, 2691, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 946, pbuffer, 2871, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1046, pbuffer, 3171, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1196, pbuffer, 3621, 210, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {418, 427});

                pbuffer.scale(2.0 * a_exp, {463, 481});

                pbuffer.scale(2.0 * a_exp, {553, 583});

                pbuffer.scale(2.0 * a_exp, {703, 748});

                pbuffer.scale(2.0 * a_exp, {928, 991});

                pbuffer.scale(2.0 * a_exp, {1261, 1279});

                pbuffer.scale(2.0 * a_exp, {1333, 1369});

                pbuffer.scale(2.0 * a_exp, {1477, 1537});

                pbuffer.scale(2.0 * a_exp, {1717, 1807});

                pbuffer.scale(2.0 * a_exp, {2077, 2203});

                pbuffer.scale(2.0 * a_exp, {2601, 2631});

                pbuffer.scale(2.0 * a_exp, {2691, 2751});

                pbuffer.scale(2.0 * a_exp, {2871, 2971});

                pbuffer.scale(2.0 * a_exp, {3171, 3321});

                pbuffer.scale(2.0 * a_exp, {3621, 3831});

                pbuffer.scale(pfactors, 0, 2.0, {4266, 4311});

                pbuffer.scale(pfactors, 0, 2.0, {4356, 4446});

                pbuffer.scale(pfactors, 0, 2.0, {4536, 4686});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {4836, 5061});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5286, 5601});

                pbuffer.scale(pfactors, 0, 2.0, {5916, 5979});

                pbuffer.scale(pfactors, 0, 2.0, {5979, 6105});

                pbuffer.scale(pfactors, 0, 2.0, {6105, 6315});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {6315, 6630});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {6630, 7071});

                t2cfunc::reduce(cbuffer, 2451, pbuffer, 418, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2460, pbuffer, 463, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2478, pbuffer, 553, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2508, pbuffer, 703, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2553, pbuffer, 928, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2616, pbuffer, 1261, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2634, pbuffer, 1333, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2670, pbuffer, 1477, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2730, pbuffer, 1717, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2820, pbuffer, 2077, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2946, pbuffer, 2601, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2976, pbuffer, 2691, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3036, pbuffer, 2871, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3136, pbuffer, 3171, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3286, pbuffer, 3621, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3496, pbuffer, 4266, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3541, pbuffer, 4356, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3631, pbuffer, 4536, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3781, pbuffer, 4836, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4006, pbuffer, 5286, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4321, pbuffer, 5916, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4384, pbuffer, 5979, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4510, pbuffer, 6105, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4720, pbuffer, 6315, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5035, pbuffer, 6630, 441, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 306, cbuffer, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 414, cbuffer, 9, 27, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 900, 306, 414, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2322, cbuffer, 57, 75, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2538, cbuffer, 75, 111, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 3510, 2322, 2538, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6150, cbuffer, 171, 201, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 6510, cbuffer, 201, 261, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 8130, 6150, 6510, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 11136, cbuffer, 1406, 1415, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 11244, cbuffer, 1415, 1433, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 11730, 11136, 11244, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 13152, cbuffer, 1463, 1481, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 13368, cbuffer, 1481, 1517, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 14340, 13152, 13368, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 16980, cbuffer, 1577, 1607, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 17340, cbuffer, 1607, 1667, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 18960, 16980, 17340, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 23190, cbuffer, 1767, 1812, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 23730, cbuffer, 1812, 1902, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 26160, 23190, 23730, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 32352, cbuffer, 2052, 2115, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 33108, cbuffer, 2115, 2241, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 36510, 32352, 33108, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 361, 370, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 27, cbuffer, 370, 388, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 81, cbuffer, 388, 418, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 171, cbuffer, 418, 463, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 333, cbuffer, 0, 0, 27, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 468, cbuffer, 9, 27, 81, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 630, cbuffer, 27, 81, 171, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 954, 306, 333, 468, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 1116, 414, 468, 630, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 1440, 900, 954, 1116, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1710, cbuffer, 526, 544, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1764, cbuffer, 544, 580, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1872, cbuffer, 580, 640, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 2052, cbuffer, 640, 730, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2376, cbuffer, 57, 1710, 1764, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2646, cbuffer, 75, 1764, 1872, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2970, cbuffer, 111, 1872, 2052, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 3618, 2322, 2376, 2646, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 3942, 2538, 2646, 2970, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 4590, 3510, 3618, 3942, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5130, cbuffer, 856, 886, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5220, cbuffer, 886, 946, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 5400, cbuffer, 946, 1046, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 5700, cbuffer, 1046, 1196, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 6240, cbuffer, 171, 5130, 5220, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 6690, cbuffer, 201, 5220, 5400, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 7230, cbuffer, 261, 5400, 5700, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 8310, 6150, 6240, 6690, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 8850, 6510, 6690, 7230, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 9930, 8130, 8310, 8850, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 10830, cbuffer, 2451, 2460, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 10857, cbuffer, 2460, 2478, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 10911, cbuffer, 2478, 2508, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 11001, cbuffer, 2508, 2553, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 11163, cbuffer, 1406, 10830, 10857, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 11298, cbuffer, 1415, 10857, 10911, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 11460, cbuffer, 1433, 10911, 11001, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 11784, 11136, 11163, 11298, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 11946, 11244, 11298, 11460, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 12270, 11730, 11784, 11946, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 12540, cbuffer, 2616, 2634, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 12594, cbuffer, 2634, 2670, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 12702, cbuffer, 2670, 2730, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 12882, cbuffer, 2730, 2820, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 13206, cbuffer, 1463, 12540, 12594, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 13476, cbuffer, 1481, 12594, 12702, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 13800, cbuffer, 1517, 12702, 12882, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 14448, 13152, 13206, 13476, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 14772, 13368, 13476, 13800, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 15420, 14340, 14448, 14772, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 15960, cbuffer, 2946, 2976, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 16050, cbuffer, 2976, 3036, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 16230, cbuffer, 3036, 3136, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 16530, cbuffer, 3136, 3286, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 17070, cbuffer, 1577, 15960, 16050, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 17520, cbuffer, 1607, 16050, 16230, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 18060, cbuffer, 1667, 16230, 16530, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 19140, 16980, 17070, 17520, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 19680, 17340, 17520, 18060, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 20760, 18960, 19140, 19680, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 21660, cbuffer, 3496, 3541, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 21795, cbuffer, 3541, 3631, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 22065, cbuffer, 3631, 3781, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 22515, cbuffer, 3781, 4006, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 23325, cbuffer, 1767, 21660, 21795, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 24000, cbuffer, 1812, 21795, 22065, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 24810, cbuffer, 1902, 22065, 22515, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 26430, 23190, 23325, 24000, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 27240, 23730, 24000, 24810, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 28860, 26160, 26430, 27240, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 30210, cbuffer, 4321, 4384, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 30399, cbuffer, 4384, 4510, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 30777, cbuffer, 4510, 4720, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 31407, cbuffer, 4720, 5035, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 32541, cbuffer, 2052, 30210, 30399, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 33486, cbuffer, 2115, 30399, 30777, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 34620, cbuffer, 2241, 30777, 31407, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 36888, 32352, 32541, 33486, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 38022, 33108, 33486, 34620, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 40290, 36510, 36888, 38022, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 0, ckbuffer, 1440, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 63, ckbuffer, 1530, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 126, ckbuffer, 1620, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 756, ckbuffer, 4590, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 882, ckbuffer, 4770, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 1008, ckbuffer, 4950, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 2268, ckbuffer, 9930, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 2478, ckbuffer, 10230, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 2688, ckbuffer, 10530, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 37107, ckbuffer, 12270, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 37170, ckbuffer, 12360, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 37233, ckbuffer, 12450, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 37296, ckbuffer, 15420, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 37422, ckbuffer, 15600, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 37548, ckbuffer, 15780, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 37674, ckbuffer, 20760, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 37884, ckbuffer, 21060, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 38094, ckbuffer, 21360, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 38304, ckbuffer, 28860, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 38619, ckbuffer, 29310, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 38934, ckbuffer, 29760, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 39249, ckbuffer, 40290, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 39690, ckbuffer, 40920, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 40131, ckbuffer, 41550, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 7623, 0, 756, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 7812, 63, 882, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 8001, 126, 1008, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 9891, 756, 2268, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 10269, 882, 2478, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 10647, 1008, 2688, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 20097, 7623, 9891, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 20475, 7812, 10269, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 20853, 8001, 10647, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 189, 37107, 37296, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 1134, 37296, 37674, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 2898, 37674, 38304, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 4788, 38304, 39249, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 8190, 0, 189, 1134, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 11025, 756, 1134, 2898, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 14427, 2268, 2898, 4788, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 21231, 7623, 8190, 11025, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 24633, 9891, 11025, 14427, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fpxx(skbuffer, 31437, 20097, 21231, 24633, r_ab, 3, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 31437, 3, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 441, skbuffer, 32067, 3, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 882, skbuffer, 32697, 3, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1323, skbuffer, 33327, 3, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1764, skbuffer, 33957, 3, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2205, skbuffer, 34587, 3, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2646, skbuffer, 35217, 3, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 3087, skbuffer, 35847, 3, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 3528, skbuffer, 36477, 3, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 3, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFPFP_hpp */
