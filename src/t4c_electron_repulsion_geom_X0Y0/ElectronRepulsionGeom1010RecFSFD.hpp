#ifndef ElectronRepulsionGeom1010RecFSFD_hpp
#define ElectronRepulsionGeom1010RecFSFD_hpp

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
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSSXX.hpp"
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
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSI.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FS|1/|r-r'||FD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fsfd(T& distributor,
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

    CSimdArray<double> pbuffer(5851, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(4995, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(44415, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(33075, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 546, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 549, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 552, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 561, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 570, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 579, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 597, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 615, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 633, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 651, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 681, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 711, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 741, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 771, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 816, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 861, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 906, 135, 220, 235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 951, 190, 280, 301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1014, 205, 301, 322, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1077, 220, 322, 343, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1140, 235, 343, 364, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1203, 301, 406, 434, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1287, 322, 434, 462, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1371, 343, 462, 490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1455, 364, 490, 518, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1539, 2, 3, 546, 549, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1545, 14, 17, 546, 552, 561, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1563, 17, 20, 549, 561, 570, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1581, 41, 47, 552, 579, 597, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1617, 47, 53, 561, 597, 615, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1653, 53, 59, 570, 615, 633, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1689, 95, 105, 597, 651, 681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1749, 105, 115, 615, 681, 711, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1809, 115, 125, 633, 711, 741, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1869, 175, 190, 681, 771, 816, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1959, 190, 205, 711, 816, 861, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2049, 205, 220, 741, 861, 906, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2139, 280, 301, 816, 951, 1014, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2265, 301, 322, 861, 1014, 1077, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2391, 322, 343, 906, 1077, 1140, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 2517, 406, 434, 1014, 1203, 1287, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 2685, 434, 462, 1077, 1287, 1371, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 2853, 462, 490, 1140, 1371, 1455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3021, 552, 561, 1539, 1545, 1563, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3051, 579, 597, 1545, 1581, 1617, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3111, 597, 615, 1563, 1617, 1653, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3171, 651, 681, 1617, 1689, 1749, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3271, 681, 711, 1653, 1749, 1809, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3371, 771, 816, 1749, 1869, 1959, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3521, 816, 861, 1809, 1959, 2049, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3671, 951, 1014, 1959, 2139, 2265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3881, 1014, 1077, 2049, 2265, 2391, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 4091, 1203, 1287, 2265, 2517, 2685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 4371, 1287, 1371, 2391, 2685, 2853, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4651, 1581, 1617, 3021, 3051, 3111, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4741, 1689, 1749, 3111, 3171, 3271, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4891, 1869, 1959, 3271, 3371, 3521, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 5116, 2139, 2265, 3521, 3671, 3881, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 5431, 2517, 2685, 3881, 4091, 4371, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 41, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 31, pbuffer, 579, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 49, pbuffer, 651, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 79, pbuffer, 771, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 124, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 1689, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 220, pbuffer, 1869, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {41, 47});

                pbuffer.scale(2.0 * a_exp, {95, 105});

                pbuffer.scale(2.0 * a_exp, {175, 190});

                pbuffer.scale(2.0 * a_exp, {579, 597});

                pbuffer.scale(2.0 * a_exp, {651, 681});

                pbuffer.scale(2.0 * a_exp, {771, 816});

                pbuffer.scale(2.0 * a_exp, {1581, 1617});

                pbuffer.scale(2.0 * a_exp, {1689, 1749});

                pbuffer.scale(2.0 * a_exp, {1869, 1959});

                pbuffer.scale(2.0 * a_exp, {3051, 3111});

                pbuffer.scale(2.0 * a_exp, {3171, 3271});

                pbuffer.scale(2.0 * a_exp, {3371, 3521});

                pbuffer.scale(2.0 * a_exp, {4651, 4741});

                pbuffer.scale(2.0 * a_exp, {4741, 4891});

                pbuffer.scale(2.0 * a_exp, {4891, 5116});

                t2cfunc::reduce(cbuffer, 1110, pbuffer, 41, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1116, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1126, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1141, pbuffer, 579, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1159, pbuffer, 651, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1189, pbuffer, 771, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1234, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1270, pbuffer, 1689, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1330, pbuffer, 1869, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1420, pbuffer, 3051, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1480, pbuffer, 3171, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1580, pbuffer, 3371, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1730, pbuffer, 4651, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1820, pbuffer, 4741, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1970, pbuffer, 4891, 225, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {41, 47});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {95, 105});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {175, 190});

                pbuffer.scale(pfactors, 0, 2.0, {280, 301});

                pbuffer.scale(pfactors, 0, 2.0, {406, 434});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {579, 597});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {651, 681});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {771, 816});

                pbuffer.scale(pfactors, 0, 2.0, {951, 1014});

                pbuffer.scale(pfactors, 0, 2.0, {1203, 1287});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1581, 1617});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1689, 1749});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1869, 1959});

                pbuffer.scale(pfactors, 0, 2.0, {2139, 2265});

                pbuffer.scale(pfactors, 0, 2.0, {2517, 2685});

                t2cfunc::reduce(cbuffer, 310, pbuffer, 41, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 316, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 326, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 341, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 362, pbuffer, 406, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 390, pbuffer, 579, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 408, pbuffer, 651, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 438, pbuffer, 771, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 483, pbuffer, 951, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 546, pbuffer, 1203, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 630, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 666, pbuffer, 1689, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 726, pbuffer, 1869, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 816, pbuffer, 2139, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 942, pbuffer, 2517, 168, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {41, 47});

                pbuffer.scale(2.0 * a_exp, {95, 105});

                pbuffer.scale(2.0 * a_exp, {175, 190});

                pbuffer.scale(2.0 * a_exp, {280, 301});

                pbuffer.scale(2.0 * a_exp, {406, 434});

                pbuffer.scale(2.0 * a_exp, {579, 597});

                pbuffer.scale(2.0 * a_exp, {651, 681});

                pbuffer.scale(2.0 * a_exp, {771, 816});

                pbuffer.scale(2.0 * a_exp, {951, 1014});

                pbuffer.scale(2.0 * a_exp, {1203, 1287});

                pbuffer.scale(2.0 * a_exp, {1581, 1617});

                pbuffer.scale(2.0 * a_exp, {1689, 1749});

                pbuffer.scale(2.0 * a_exp, {1869, 1959});

                pbuffer.scale(2.0 * a_exp, {2139, 2265});

                pbuffer.scale(2.0 * a_exp, {2517, 2685});

                pbuffer.scale(pfactors, 0, 2.0, {3051, 3111});

                pbuffer.scale(pfactors, 0, 2.0, {3171, 3271});

                pbuffer.scale(pfactors, 0, 2.0, {3371, 3521});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3671, 3881});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {4091, 4371});

                pbuffer.scale(pfactors, 0, 2.0, {4651, 4741});

                pbuffer.scale(pfactors, 0, 2.0, {4741, 4891});

                pbuffer.scale(pfactors, 0, 2.0, {4891, 5116});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5116, 5431});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5431, 5851});

                t2cfunc::reduce(cbuffer, 2195, pbuffer, 41, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2201, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2211, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2226, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2247, pbuffer, 406, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2275, pbuffer, 579, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2293, pbuffer, 651, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2323, pbuffer, 771, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2368, pbuffer, 951, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2431, pbuffer, 1203, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2515, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2551, pbuffer, 1689, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2611, pbuffer, 1869, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2701, pbuffer, 2139, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2827, pbuffer, 2517, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2995, pbuffer, 3051, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3055, pbuffer, 3171, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3155, pbuffer, 3371, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3305, pbuffer, 3671, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3515, pbuffer, 4091, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3795, pbuffer, 4651, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3885, pbuffer, 4741, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4035, pbuffer, 4891, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4260, pbuffer, 5116, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4575, pbuffer, 5431, 420, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 156, cbuffer, 0, 6, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 228, cbuffer, 6, 16, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 483, 156, 228, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1455, cbuffer, 31, 49, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1671, cbuffer, 49, 79, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 2436, 1455, 1671, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4884, cbuffer, 124, 160, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5316, cbuffer, 160, 220, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 6846, 4884, 5316, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 10026, cbuffer, 1110, 1116, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 10098, cbuffer, 1116, 1126, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 10353, 10026, 10098, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 11325, cbuffer, 1141, 1159, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 11541, cbuffer, 1159, 1189, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 12306, 11325, 11541, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 14754, cbuffer, 1234, 1270, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 15186, cbuffer, 1270, 1330, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 16716, 14754, 15186, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 21300, cbuffer, 1420, 1480, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 22020, cbuffer, 1480, 1580, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 24570, 21300, 22020, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 31950, cbuffer, 1730, 1820, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 33030, cbuffer, 1820, 1970, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 36855, 31950, 33030, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 0, cbuffer, 310, 316, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 18, cbuffer, 316, 326, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 48, cbuffer, 326, 341, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 93, cbuffer, 341, 362, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 174, cbuffer, 0, 0, 18, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 258, cbuffer, 6, 18, 48, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 348, cbuffer, 16, 48, 93, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 519, 156, 174, 258, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 627, 228, 258, 348, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 807, 483, 519, 627, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 987, cbuffer, 390, 408, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1041, cbuffer, 408, 438, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 1131, cbuffer, 438, 483, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 1266, cbuffer, 483, 546, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1509, cbuffer, 31, 987, 1041, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 1761, cbuffer, 49, 1041, 1131, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 2031, cbuffer, 79, 1131, 1266, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 2544, 1455, 1509, 1761, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 2868, 1671, 1761, 2031, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 3408, 2436, 2544, 2868, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3948, cbuffer, 630, 666, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 4056, cbuffer, 666, 726, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 4236, cbuffer, 726, 816, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 4506, cbuffer, 816, 942, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4992, cbuffer, 124, 3948, 4056, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 5496, cbuffer, 160, 4056, 4236, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 6036, cbuffer, 220, 4236, 4506, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 7062, 4884, 4992, 5496, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 7710, 5316, 5496, 6036, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 8790, 6846, 7062, 7710, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9870, cbuffer, 2195, 2201, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 9888, cbuffer, 2201, 2211, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 9918, cbuffer, 2211, 2226, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 9963, cbuffer, 2226, 2247, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 10044, cbuffer, 1110, 9870, 9888, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 10128, cbuffer, 1116, 9888, 9918, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 10218, cbuffer, 1126, 9918, 9963, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 10389, 10026, 10044, 10128, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 10497, 10098, 10128, 10218, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 10677, 10353, 10389, 10497, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 10857, cbuffer, 2275, 2293, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 10911, cbuffer, 2293, 2323, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 11001, cbuffer, 2323, 2368, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 11136, cbuffer, 2368, 2431, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 11379, cbuffer, 1141, 10857, 10911, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 11631, cbuffer, 1159, 10911, 11001, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 11901, cbuffer, 1189, 11001, 11136, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 12414, 11325, 11379, 11631, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 12738, 11541, 11631, 11901, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 13278, 12306, 12414, 12738, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 13818, cbuffer, 2515, 2551, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 13926, cbuffer, 2551, 2611, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 14106, cbuffer, 2611, 2701, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 14376, cbuffer, 2701, 2827, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 14862, cbuffer, 1234, 13818, 13926, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 15366, cbuffer, 1270, 13926, 14106, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 15906, cbuffer, 1330, 14106, 14376, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 16932, 14754, 14862, 15366, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 17580, 15186, 15366, 15906, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 18660, 16716, 16932, 17580, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 19740, cbuffer, 2995, 3055, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 19920, cbuffer, 3055, 3155, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 20220, cbuffer, 3155, 3305, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 20670, cbuffer, 3305, 3515, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 21480, cbuffer, 1420, 19740, 19920, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 22320, cbuffer, 1480, 19920, 20220, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 23220, cbuffer, 1580, 20220, 20670, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 24930, 21300, 21480, 22320, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 26010, 22020, 22320, 23220, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 27810, 24570, 24930, 26010, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 29610, cbuffer, 3795, 3885, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 29880, cbuffer, 3885, 4035, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 30330, cbuffer, 4035, 4260, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 31005, cbuffer, 4260, 4575, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 32220, cbuffer, 1730, 29610, 29880, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 33480, cbuffer, 1820, 29880, 30330, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 34830, cbuffer, 1970, 30330, 31005, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 37395, 31950, 32220, 33480, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 39015, 33030, 33480, 34830, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 41715, 36855, 37395, 39015, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 0, ckbuffer, 807, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 35, ckbuffer, 867, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 70, ckbuffer, 927, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 420, ckbuffer, 3408, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 525, ckbuffer, 3588, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 630, ckbuffer, 3768, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 1680, ckbuffer, 8790, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 1890, ckbuffer, 9150, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 2100, ckbuffer, 9510, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 29400, ckbuffer, 10677, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 29435, ckbuffer, 10737, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 29470, ckbuffer, 10797, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 29505, ckbuffer, 13278, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 29610, ckbuffer, 13458, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 29715, ckbuffer, 13638, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 29820, ckbuffer, 18660, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 30030, ckbuffer, 19020, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 30240, ckbuffer, 19380, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 30450, ckbuffer, 27810, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 30800, ckbuffer, 28410, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 31150, ckbuffer, 29010, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 31500, ckbuffer, 41715, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 32025, ckbuffer, 42615, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 32550, ckbuffer, 43515, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 7350, 0, 420, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 7455, 35, 525, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 7560, 70, 630, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 8610, 420, 1680, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 8925, 525, 1890, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 9240, 630, 2100, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 18060, 7350, 8610, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 18270, 7455, 8925, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 18480, 7560, 9240, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 105, 29400, 29505, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 735, 29505, 29820, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 2310, 29820, 30450, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 4200, 30450, 31500, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 7665, 0, 105, 735, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 9555, 420, 735, 2310, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 12390, 1680, 2310, 4200, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 18690, 7350, 7665, 9555, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 20580, 8610, 9555, 12390, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fsxx(skbuffer, 26250, 18060, 18690, 20580, r_ab, 3, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 26250, 3, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 245, skbuffer, 26600, 3, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 490, skbuffer, 26950, 3, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 735, skbuffer, 27300, 3, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 980, skbuffer, 27650, 3, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1225, skbuffer, 28000, 3, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1470, skbuffer, 28350, 3, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1715, skbuffer, 28700, 3, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1960, skbuffer, 29050, 3, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 3, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFSFD_hpp */
