#ifndef ElectronRepulsionGeom1010RecDSFF_hpp
#define ElectronRepulsionGeom1010RecDSFF_hpp

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
#include "ElectronRepulsionGeom1010ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSSXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSI.hpp"
#include "ElectronRepulsionPrimRecSDSK.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSI.hpp"
#include "ElectronRepulsionPrimRecSFSK.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
#include "ElectronRepulsionPrimRecSPSI.hpp"
#include "ElectronRepulsionPrimRecSPSK.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DS|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dsff(T& distributor,
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

    CSimdArray<double> pbuffer(4181, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(3744, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(36504, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(16317, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2205, 1);

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

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 406, 175, 190, 280, 301, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 434, 190, 205, 301, 322, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 462, 205, 220, 322, 343, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 490, 220, 235, 343, 364, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 518, 235, 250, 364, 385, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 546, 280, 301, 406, 434, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 582, 301, 322, 434, 462, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 618, 322, 343, 462, 490, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 654, 343, 364, 490, 518, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 690, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 699, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 717, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 735, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 765, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 795, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 825, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 870, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 915, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 960, 190, 280, 301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1023, 205, 301, 322, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1086, 220, 322, 343, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1149, 301, 406, 434, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1233, 322, 434, 462, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1317, 343, 462, 490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 1401, 434, 546, 582, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 1509, 462, 582, 618, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 1617, 490, 618, 654, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1725, 47, 53, 690, 699, 717, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1761, 95, 105, 699, 735, 765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1821, 105, 115, 717, 765, 795, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1881, 175, 190, 765, 825, 870, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1971, 190, 205, 795, 870, 915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2061, 280, 301, 870, 960, 1023, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2187, 301, 322, 915, 1023, 1086, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 2313, 406, 434, 1023, 1149, 1233, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 2481, 434, 462, 1086, 1233, 1317, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 2649, 546, 582, 1233, 1401, 1509, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 2865, 582, 618, 1317, 1509, 1617, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3081, 735, 765, 1725, 1761, 1821, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3181, 825, 870, 1821, 1881, 1971, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3331, 960, 1023, 1971, 2061, 2187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 3541, 1149, 1233, 2187, 2313, 2481, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 3821, 1401, 1509, 2481, 2649, 2865, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 25, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 46, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 76, pbuffer, 825, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 121, pbuffer, 960, 63, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {95, 105});

                pbuffer.scale(2.0 * a_exp, {175, 190});

                pbuffer.scale(2.0 * a_exp, {280, 301});

                pbuffer.scale(2.0 * a_exp, {735, 765});

                pbuffer.scale(2.0 * a_exp, {825, 870});

                pbuffer.scale(2.0 * a_exp, {960, 1023});

                pbuffer.scale(2.0 * a_exp, {1761, 1821});

                pbuffer.scale(2.0 * a_exp, {1881, 1971});

                pbuffer.scale(2.0 * a_exp, {2061, 2187});

                pbuffer.scale(2.0 * a_exp, {3081, 3181});

                pbuffer.scale(2.0 * a_exp, {3181, 3331});

                pbuffer.scale(2.0 * a_exp, {3331, 3541});

                t2cfunc::reduce(cbuffer, 624, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 634, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 649, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 670, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 700, pbuffer, 825, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 745, pbuffer, 960, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 808, pbuffer, 1761, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 868, pbuffer, 1881, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 958, pbuffer, 2061, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1084, pbuffer, 3081, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1184, pbuffer, 3181, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1334, pbuffer, 3331, 210, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {95, 105});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {175, 190});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {280, 301});

                pbuffer.scale(pfactors, 0, 2.0, {406, 434});

                pbuffer.scale(pfactors, 0, 2.0, {546, 582});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {735, 765});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {825, 870});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {960, 1023});

                pbuffer.scale(pfactors, 0, 2.0, {1149, 1233});

                pbuffer.scale(pfactors, 0, 2.0, {1401, 1509});

                t2cfunc::reduce(cbuffer, 184, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 194, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 209, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 230, pbuffer, 406, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 258, pbuffer, 546, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 294, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 324, pbuffer, 825, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 369, pbuffer, 960, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 432, pbuffer, 1149, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 516, pbuffer, 1401, 108, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {95, 105});

                pbuffer.scale(2.0 * a_exp, {175, 190});

                pbuffer.scale(2.0 * a_exp, {280, 301});

                pbuffer.scale(2.0 * a_exp, {406, 434});

                pbuffer.scale(2.0 * a_exp, {546, 582});

                pbuffer.scale(2.0 * a_exp, {735, 765});

                pbuffer.scale(2.0 * a_exp, {825, 870});

                pbuffer.scale(2.0 * a_exp, {960, 1023});

                pbuffer.scale(2.0 * a_exp, {1149, 1233});

                pbuffer.scale(2.0 * a_exp, {1401, 1509});

                pbuffer.scale(pfactors, 0, 2.0, {1761, 1821});

                pbuffer.scale(pfactors, 0, 2.0, {1881, 1971});

                pbuffer.scale(pfactors, 0, 2.0, {2061, 2187});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2313, 2481});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2649, 2865});

                pbuffer.scale(pfactors, 0, 2.0, {3081, 3181});

                pbuffer.scale(pfactors, 0, 2.0, {3181, 3331});

                pbuffer.scale(pfactors, 0, 2.0, {3331, 3541});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3541, 3821});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3821, 4181});

                t2cfunc::reduce(cbuffer, 1544, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1554, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1569, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1590, pbuffer, 406, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1618, pbuffer, 546, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1654, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1684, pbuffer, 825, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1729, pbuffer, 960, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1792, pbuffer, 1149, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1876, pbuffer, 1401, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1984, pbuffer, 1761, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2044, pbuffer, 1881, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2134, pbuffer, 2061, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2260, pbuffer, 2313, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2428, pbuffer, 2649, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2644, pbuffer, 3081, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2744, pbuffer, 3181, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2894, pbuffer, 3331, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3104, pbuffer, 3541, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3384, pbuffer, 3821, 360, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 222, cbuffer, 0, 10, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 342, cbuffer, 10, 25, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 711, 222, 342, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2187, cbuffer, 46, 76, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 2547, cbuffer, 76, 121, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 3654, 2187, 2547, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 6306, cbuffer, 624, 634, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 6426, cbuffer, 634, 649, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 6795, 6306, 6426, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 8271, cbuffer, 670, 700, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 8631, cbuffer, 700, 745, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 9738, 8271, 8631, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 13500, cbuffer, 808, 868, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 14220, cbuffer, 868, 958, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 16434, 13500, 14220, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 23514, cbuffer, 1084, 1184, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 24714, cbuffer, 1184, 1334, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 28404, 23514, 24714, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 0, cbuffer, 184, 194, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 30, cbuffer, 194, 209, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 75, cbuffer, 209, 230, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 138, cbuffer, 230, 258, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 252, cbuffer, 0, 0, 30, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 387, cbuffer, 10, 30, 75, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 522, cbuffer, 25, 75, 138, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 771, 222, 252, 387, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 951, 342, 387, 522, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 1221, 711, 771, 951, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1521, cbuffer, 294, 324, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 1611, cbuffer, 324, 369, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 1746, cbuffer, 369, 432, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 1935, cbuffer, 432, 516, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2277, cbuffer, 46, 1521, 1611, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 2682, cbuffer, 76, 1611, 1746, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 3087, cbuffer, 121, 1746, 1935, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 3834, 2187, 2277, 2682, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 4374, 2547, 2682, 3087, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 5184, 3654, 3834, 4374, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6084, cbuffer, 1544, 1554, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 6114, cbuffer, 1554, 1569, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 6159, cbuffer, 1569, 1590, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 6222, cbuffer, 1590, 1618, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 6336, cbuffer, 624, 6084, 6114, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 6471, cbuffer, 634, 6114, 6159, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 6606, cbuffer, 649, 6159, 6222, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 6855, 6306, 6336, 6471, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 7035, 6426, 6471, 6606, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 7305, 6795, 6855, 7035, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 7605, cbuffer, 1654, 1684, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 7695, cbuffer, 1684, 1729, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 7830, cbuffer, 1729, 1792, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 8019, cbuffer, 1792, 1876, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 8361, cbuffer, 670, 7605, 7695, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 8766, cbuffer, 700, 7695, 7830, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 9171, cbuffer, 745, 7830, 8019, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 9918, 8271, 8361, 8766, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 10458, 8631, 8766, 9171, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 11268, 9738, 9918, 10458, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 12168, cbuffer, 1984, 2044, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 12348, cbuffer, 2044, 2134, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 12618, cbuffer, 2134, 2260, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 12996, cbuffer, 2260, 2428, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 13680, cbuffer, 808, 12168, 12348, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 14490, cbuffer, 868, 12348, 12618, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 15300, cbuffer, 958, 12618, 12996, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 16794, 13500, 13680, 14490, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 17874, 14220, 14490, 15300, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 19494, 16434, 16794, 17874, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 21294, cbuffer, 2644, 2744, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 21594, cbuffer, 2744, 2894, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 22044, cbuffer, 2894, 3104, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 22674, cbuffer, 3104, 3384, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 23814, cbuffer, 1084, 21294, 21594, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 25164, cbuffer, 1184, 21594, 22044, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 26514, cbuffer, 1334, 22044, 22674, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 29004, 23514, 23814, 25164, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 30804, 24714, 25164, 26514, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 33504, 28404, 29004, 30804, cfactors, 6, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 1221, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 49, ckbuffer, 1321, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 98, ckbuffer, 1421, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 588, ckbuffer, 5184, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 735, ckbuffer, 5484, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 882, ckbuffer, 5784, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 13377, ckbuffer, 7305, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 13426, ckbuffer, 7405, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 13475, ckbuffer, 7505, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 13524, ckbuffer, 11268, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 13671, ckbuffer, 11568, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 13818, ckbuffer, 11868, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 13965, ckbuffer, 19494, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 14259, ckbuffer, 20094, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 14553, ckbuffer, 20694, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 14847, ckbuffer, 33504, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 15337, ckbuffer, 34504, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 15827, ckbuffer, 35504, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 4998, 0, 588, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 5145, 49, 735, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 5292, 98, 882, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 147, 13377, 13524, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 1029, 13524, 13965, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 2352, 13965, 14847, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 5439, 0, 147, 1029, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 6762, 588, 1029, 2352, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 10731, 4998, 5439, 6762, r_ab, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 0, skbuffer, 10731, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 245, skbuffer, 11025, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 490, skbuffer, 11319, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 735, skbuffer, 11613, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 980, skbuffer, 11907, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 1225, skbuffer, 12201, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 1470, skbuffer, 12495, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 1715, skbuffer, 12789, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 1960, skbuffer, 13083, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 0, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDSFF_hpp */
