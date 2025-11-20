#ifndef ElectronRepulsionGeom1010RecDPFD_hpp
#define ElectronRepulsionGeom1010RecDPFD_hpp

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
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DP|1/|r-r'||FD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dpfd(T& distributor,
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

    CSimdArray<double> cbuffer(4773, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(42441, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(25620, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(4725, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 579, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 651, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 48, pbuffer, 771, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 93, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 129, pbuffer, 1689, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 189, pbuffer, 1869, 90, ket_width, ket_npgtos);

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

                t2cfunc::reduce(cbuffer, 999, pbuffer, 579, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1017, pbuffer, 651, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1047, pbuffer, 771, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1092, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1128, pbuffer, 1689, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1188, pbuffer, 1869, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1278, pbuffer, 3051, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1338, pbuffer, 3171, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1438, pbuffer, 3371, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1588, pbuffer, 4651, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1678, pbuffer, 4741, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1828, pbuffer, 4891, 225, ket_width, ket_npgtos);

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

                t2cfunc::reduce(cbuffer, 279, pbuffer, 579, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 297, pbuffer, 651, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 327, pbuffer, 771, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 372, pbuffer, 951, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 435, pbuffer, 1203, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 519, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 555, pbuffer, 1689, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 615, pbuffer, 1869, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 705, pbuffer, 2139, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 831, pbuffer, 2517, 168, ket_width, ket_npgtos);

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

                t2cfunc::reduce(cbuffer, 2053, pbuffer, 579, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2071, pbuffer, 651, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2101, pbuffer, 771, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2146, pbuffer, 951, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2209, pbuffer, 1203, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2293, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2329, pbuffer, 1689, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2389, pbuffer, 1869, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2479, pbuffer, 2139, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2605, pbuffer, 2517, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2773, pbuffer, 3051, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2833, pbuffer, 3171, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2933, pbuffer, 3371, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3083, pbuffer, 3671, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3293, pbuffer, 4091, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3573, pbuffer, 4651, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3663, pbuffer, 4741, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3813, pbuffer, 4891, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4038, pbuffer, 5116, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4353, pbuffer, 5431, 420, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 468, cbuffer, 0, 18, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 684, cbuffer, 18, 48, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 1449, 468, 684, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3897, cbuffer, 93, 129, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4329, cbuffer, 129, 189, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 5859, 3897, 4329, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 9351, cbuffer, 999, 1017, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9567, cbuffer, 1017, 1047, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 10332, 9351, 9567, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 12780, cbuffer, 1092, 1128, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 13212, cbuffer, 1128, 1188, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 14742, 12780, 13212, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 19326, cbuffer, 1278, 1338, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 20046, cbuffer, 1338, 1438, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 22596, 19326, 20046, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 29976, cbuffer, 1588, 1678, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 31056, cbuffer, 1678, 1828, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 34881, 29976, 31056, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 0, cbuffer, 279, 297, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 54, cbuffer, 297, 327, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 144, cbuffer, 327, 372, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 279, cbuffer, 372, 435, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 522, cbuffer, 0, 0, 54, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 774, cbuffer, 18, 54, 144, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 1044, cbuffer, 48, 144, 279, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 1557, 468, 522, 774, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 1881, 684, 774, 1044, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 2421, 1449, 1557, 1881, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2961, cbuffer, 519, 555, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3069, cbuffer, 555, 615, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 3249, cbuffer, 615, 705, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 3519, cbuffer, 705, 831, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4005, cbuffer, 93, 2961, 3069, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 4509, cbuffer, 129, 3069, 3249, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 5049, cbuffer, 189, 3249, 3519, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 6075, 3897, 4005, 4509, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 6723, 4329, 4509, 5049, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 7803, 5859, 6075, 6723, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 8883, cbuffer, 2053, 2071, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 8937, cbuffer, 2071, 2101, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 9027, cbuffer, 2101, 2146, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 9162, cbuffer, 2146, 2209, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 9405, cbuffer, 999, 8883, 8937, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 9657, cbuffer, 1017, 8937, 9027, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 9927, cbuffer, 1047, 9027, 9162, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 10440, 9351, 9405, 9657, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 10764, 9567, 9657, 9927, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 11304, 10332, 10440, 10764, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 11844, cbuffer, 2293, 2329, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 11952, cbuffer, 2329, 2389, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 12132, cbuffer, 2389, 2479, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 12402, cbuffer, 2479, 2605, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 12888, cbuffer, 1092, 11844, 11952, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 13392, cbuffer, 1128, 11952, 12132, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 13932, cbuffer, 1188, 12132, 12402, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 14958, 12780, 12888, 13392, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 15606, 13212, 13392, 13932, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 16686, 14742, 14958, 15606, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 17766, cbuffer, 2773, 2833, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 17946, cbuffer, 2833, 2933, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 18246, cbuffer, 2933, 3083, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 18696, cbuffer, 3083, 3293, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 19506, cbuffer, 1278, 17766, 17946, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 20346, cbuffer, 1338, 17946, 18246, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 21246, cbuffer, 1438, 18246, 18696, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 22956, 19326, 19506, 20346, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 24036, 20046, 20346, 21246, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 25836, 22596, 22956, 24036, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 27636, cbuffer, 3573, 3663, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 27906, cbuffer, 3663, 3813, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 28356, cbuffer, 3813, 4038, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 29031, cbuffer, 4038, 4353, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 30246, cbuffer, 1588, 27636, 27906, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 31506, cbuffer, 1678, 27906, 28356, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 32856, cbuffer, 1828, 28356, 29031, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 35421, 29976, 30246, 31506, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 37041, 31056, 31506, 32856, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 39741, 34881, 35421, 37041, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 0, ckbuffer, 2421, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 105, ckbuffer, 2601, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 210, ckbuffer, 2781, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 1260, ckbuffer, 7803, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 1470, ckbuffer, 8163, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 1680, ckbuffer, 8523, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 22050, ckbuffer, 11304, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 22155, ckbuffer, 11484, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 22260, ckbuffer, 11664, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 22365, ckbuffer, 16686, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 22575, ckbuffer, 17046, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 22785, ckbuffer, 17406, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 22995, ckbuffer, 25836, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 23345, ckbuffer, 26436, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 23695, ckbuffer, 27036, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 24045, ckbuffer, 39741, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 24570, ckbuffer, 40641, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 25095, ckbuffer, 41541, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 6930, 0, 1260, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 7245, 105, 1470, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 7560, 210, 1680, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 315, 22050, 22365, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 1890, 22365, 22995, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 3780, 22995, 24045, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 7875, 0, 315, 1890, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 10710, 1260, 1890, 3780, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 16380, 6930, 7875, 10710, r_ab, 3, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 16380, 3, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 525, skbuffer, 17010, 3, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1050, skbuffer, 17640, 3, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1575, skbuffer, 18270, 3, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2100, skbuffer, 18900, 3, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2625, skbuffer, 19530, 3, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 3150, skbuffer, 20160, 3, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 3675, skbuffer, 20790, 3, 2);

            t4cfunc::bra_transform<2, 1>(sbuffer, 4200, skbuffer, 21420, 3, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 3, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDPFD_hpp */
