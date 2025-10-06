#ifndef ElectronRepulsionGeom2000RecDSFF_hpp
#define ElectronRepulsionGeom2000RecDSFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDDXX.hpp"
#include "ElectronRepulsionContrRecDPXX.hpp"
#include "ElectronRepulsionContrRecDSXX.hpp"
#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPSXX.hpp"
#include "ElectronRepulsionContrRecXXDF.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXFF.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
#include "ElectronRepulsionGeom1000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSI.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSI.hpp"
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

/// @brief Computes d^(2)/dA^(2)(DS|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_dsff(T& distributor,
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

    CSimdArray<double> pbuffer(5581, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(3404, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(17848, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(17395, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1470, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 546, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 549, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 558, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 567, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 585, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 603, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 621, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 651, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 681, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 711, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 741, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 786, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 831, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 876, 135, 220, 235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 921, 190, 280, 301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 984, 205, 301, 322, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1047, 220, 322, 343, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1110, 235, 343, 364, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1173, 301, 406, 434, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1257, 322, 434, 462, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1341, 343, 462, 490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1425, 364, 490, 518, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1509, 17, 20, 546, 549, 558, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1527, 47, 53, 549, 567, 585, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1563, 53, 59, 558, 585, 603, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1599, 95, 105, 567, 621, 651, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1659, 105, 115, 585, 651, 681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1719, 115, 125, 603, 681, 711, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1779, 175, 190, 651, 741, 786, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1869, 190, 205, 681, 786, 831, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1959, 205, 220, 711, 831, 876, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2049, 280, 301, 786, 921, 984, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2175, 301, 322, 831, 984, 1047, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2301, 322, 343, 876, 1047, 1110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 2427, 406, 434, 984, 1173, 1257, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 2595, 434, 462, 1047, 1257, 1341, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 2763, 462, 490, 1110, 1341, 1425, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2931, 567, 585, 1509, 1527, 1563, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2991, 621, 651, 1527, 1599, 1659, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3091, 651, 681, 1563, 1659, 1719, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3191, 741, 786, 1659, 1779, 1869, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3341, 786, 831, 1719, 1869, 1959, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3491, 921, 984, 1869, 2049, 2175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3701, 984, 1047, 1959, 2175, 2301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 3911, 1173, 1257, 2175, 2427, 2595, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 4191, 1257, 1341, 2301, 2595, 2763, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4471, 1599, 1659, 2931, 2991, 3091, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4621, 1779, 1869, 3091, 3191, 3341, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 4846, 2049, 2175, 3341, 3491, 3701, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 5161, 2427, 2595, 3701, 3911, 4191, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 25, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 46, pbuffer, 406, 28, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {95, 105});

                pbuffer.scale(2.0 * a_exp, {175, 190});

                pbuffer.scale(2.0 * a_exp, {280, 301});

                pbuffer.scale(2.0 * a_exp, {406, 434});

                pbuffer.scale(2.0 * a_exp, {621, 651});

                pbuffer.scale(2.0 * a_exp, {741, 786});

                pbuffer.scale(2.0 * a_exp, {921, 984});

                pbuffer.scale(2.0 * a_exp, {1173, 1257});

                pbuffer.scale(2.0 * a_exp, {1599, 1659});

                pbuffer.scale(2.0 * a_exp, {1779, 1869});

                pbuffer.scale(2.0 * a_exp, {2049, 2175});

                pbuffer.scale(2.0 * a_exp, {2427, 2595});

                t2cfunc::reduce(cbuffer, 74, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 84, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 99, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 120, pbuffer, 406, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 148, pbuffer, 621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 178, pbuffer, 741, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 223, pbuffer, 921, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 286, pbuffer, 1173, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 370, pbuffer, 1599, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 430, pbuffer, 1779, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 520, pbuffer, 2049, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 646, pbuffer, 2427, 168, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {95, 105});

                pbuffer.scale(2.0 * a_exp, {175, 190});

                pbuffer.scale(2.0 * a_exp, {280, 301});

                pbuffer.scale(2.0 * a_exp, {406, 434});

                pbuffer.scale(2.0 * a_exp, {621, 651});

                pbuffer.scale(2.0 * a_exp, {741, 786});

                pbuffer.scale(2.0 * a_exp, {921, 984});

                pbuffer.scale(2.0 * a_exp, {1173, 1257});

                pbuffer.scale(2.0 * a_exp, {1599, 1659});

                pbuffer.scale(2.0 * a_exp, {1779, 1869});

                pbuffer.scale(2.0 * a_exp, {2049, 2175});

                pbuffer.scale(2.0 * a_exp, {2427, 2595});

                pbuffer.scale(4.0 * a_exp * a_exp, {2991, 3091});

                pbuffer.scale(4.0 * a_exp * a_exp, {3191, 3341});

                pbuffer.scale(4.0 * a_exp * a_exp, {3491, 3701});

                pbuffer.scale(4.0 * a_exp * a_exp, {3911, 4191});

                pbuffer.scale(4.0 * a_exp * a_exp, {4471, 4621});

                pbuffer.scale(4.0 * a_exp * a_exp, {4621, 4846});

                pbuffer.scale(4.0 * a_exp * a_exp, {4846, 5161});

                pbuffer.scale(4.0 * a_exp * a_exp, {5161, 5581});

                t2cfunc::reduce(cbuffer, 814, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 824, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 839, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 860, pbuffer, 406, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 888, pbuffer, 621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 918, pbuffer, 741, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 963, pbuffer, 921, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1026, pbuffer, 1173, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1110, pbuffer, 1599, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1170, pbuffer, 1779, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1260, pbuffer, 2049, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1386, pbuffer, 2427, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1554, pbuffer, 2991, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1654, pbuffer, 3191, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1804, pbuffer, 3491, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2014, pbuffer, 3911, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2294, pbuffer, 4471, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2444, pbuffer, 4621, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2669, pbuffer, 4846, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2984, pbuffer, 5161, 420, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 10, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 30, cbuffer, 10, 25, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 75, cbuffer, 25, 46, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 138, 0, 30, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 198, 30, 75, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 288, 138, 198, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 388, cbuffer, 74, 84, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 418, cbuffer, 84, 99, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 463, cbuffer, 99, 120, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 526, 388, 418, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 586, 418, 463, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 676, 526, 586, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 776, cbuffer, 148, 178, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 866, cbuffer, 178, 223, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 1001, cbuffer, 223, 286, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 1190, 776, 866, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 1370, 866, 1001, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 1640, 1190, 1370, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1940, cbuffer, 370, 430, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 2120, cbuffer, 430, 520, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 2390, cbuffer, 520, 646, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 2768, 1940, 2120, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 3128, 2120, 2390, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 3668, 2768, 3128, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4268, cbuffer, 814, 824, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 4298, cbuffer, 824, 839, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 4343, cbuffer, 839, 860, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 4406, 4268, 4298, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 4466, 4298, 4343, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 4556, 4406, 4466, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4656, cbuffer, 888, 918, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 4746, cbuffer, 918, 963, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 4881, cbuffer, 963, 1026, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 5070, 4656, 4746, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 5250, 4746, 4881, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 5520, 5070, 5250, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5820, cbuffer, 1110, 1170, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 6000, cbuffer, 1170, 1260, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 6270, cbuffer, 1260, 1386, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 6648, 5820, 6000, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 7008, 6000, 6270, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 7548, 6648, 7008, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 8148, cbuffer, 1554, 1654, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 8448, cbuffer, 1654, 1804, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 8898, cbuffer, 1804, 2014, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 9528, 8148, 8448, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 10128, 8448, 8898, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 11028, 9528, 10128, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 12028, cbuffer, 2294, 2444, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 12478, cbuffer, 2444, 2669, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 13153, cbuffer, 2669, 2984, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 14098, 12028, 12478, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 14998, 12478, 13153, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 16348, 14098, 14998, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 288, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 8722, ckbuffer, 676, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 8771, ckbuffer, 1640, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 8918, ckbuffer, 3668, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 9800, ckbuffer, 4556, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 9849, ckbuffer, 5520, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 9996, ckbuffer, 7548, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 10290, ckbuffer, 11028, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 10780, ckbuffer, 16348, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 9212, 8722, 8771, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 9359, 8771, 8918, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 11515, 9800, 9849, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 11662, 9849, 9996, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 12103, 9996, 10290, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 12985, 10290, 10780, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 14455, 11515, 11662, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 14749, 11662, 12103, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 15631, 12103, 12985, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_psxx(skbuffer, 2989, 0, 9212, 9359, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 49, 8722, 14455, 0, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 343, 8771, 14749, 1, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 1225, 8918, 15631, 2, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_psxx(skbuffer, 3430, 9212, 49, 343, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_ppxx(skbuffer, 4312, 9359, 343, 1225, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dsxx(skbuffer, 6958, 2989, 3430, 4312, r_ab, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 0, skbuffer, 6958, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 245, skbuffer, 7252, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 490, skbuffer, 7546, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 735, skbuffer, 7840, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 980, skbuffer, 8134, 3, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 1225, skbuffer, 8428, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 0, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecDSFF_hpp */
