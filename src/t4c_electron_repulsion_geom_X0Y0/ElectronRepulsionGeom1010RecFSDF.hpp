#ifndef ElectronRepulsionGeom1010RecFSDF_hpp
#define ElectronRepulsionGeom1010RecFSDF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPG.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FS|1/|r-r'||DF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fsdf(T& distributor,
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

    CSimdArray<double> cbuffer(4455, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(25785, 1);

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

                t2cfunc::reduce(cbuffer, 25, pbuffer, 621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 55, pbuffer, 741, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 1599, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 1779, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {95, 105});

                pbuffer.scale(2.0 * a_exp, {175, 190});

                pbuffer.scale(2.0 * a_exp, {621, 651});

                pbuffer.scale(2.0 * a_exp, {741, 786});

                pbuffer.scale(2.0 * a_exp, {1599, 1659});

                pbuffer.scale(2.0 * a_exp, {1779, 1869});

                pbuffer.scale(2.0 * a_exp, {2991, 3091});

                pbuffer.scale(2.0 * a_exp, {3191, 3341});

                pbuffer.scale(2.0 * a_exp, {4471, 4621});

                pbuffer.scale(2.0 * a_exp, {4621, 4846});

                t2cfunc::reduce(cbuffer, 990, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1000, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1015, pbuffer, 621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1045, pbuffer, 741, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1090, pbuffer, 1599, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1150, pbuffer, 1779, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1240, pbuffer, 2991, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1340, pbuffer, 3191, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1490, pbuffer, 4471, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1640, pbuffer, 4621, 225, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {95, 105});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {175, 190});

                pbuffer.scale(pfactors, 0, 2.0, {280, 301});

                pbuffer.scale(pfactors, 0, 2.0, {406, 434});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {621, 651});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {741, 786});

                pbuffer.scale(pfactors, 0, 2.0, {921, 984});

                pbuffer.scale(pfactors, 0, 2.0, {1173, 1257});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1599, 1659});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1779, 1869});

                pbuffer.scale(pfactors, 0, 2.0, {2049, 2175});

                pbuffer.scale(pfactors, 0, 2.0, {2427, 2595});

                t2cfunc::reduce(cbuffer, 250, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 260, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 275, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 296, pbuffer, 406, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 324, pbuffer, 621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 354, pbuffer, 741, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 399, pbuffer, 921, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 462, pbuffer, 1173, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 546, pbuffer, 1599, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 606, pbuffer, 1779, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 696, pbuffer, 2049, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 822, pbuffer, 2427, 168, ket_width, ket_npgtos);

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

                pbuffer.scale(pfactors, 0, 2.0, {2991, 3091});

                pbuffer.scale(pfactors, 0, 2.0, {3191, 3341});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3491, 3701});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3911, 4191});

                pbuffer.scale(pfactors, 0, 2.0, {4471, 4621});

                pbuffer.scale(pfactors, 0, 2.0, {4621, 4846});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {4846, 5161});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5161, 5581});

                t2cfunc::reduce(cbuffer, 1865, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1875, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1890, pbuffer, 280, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1911, pbuffer, 406, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1939, pbuffer, 621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1969, pbuffer, 741, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2014, pbuffer, 921, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2077, pbuffer, 1173, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2161, pbuffer, 1599, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2221, pbuffer, 1779, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2311, pbuffer, 2049, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2437, pbuffer, 2427, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2605, pbuffer, 2991, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2705, pbuffer, 3191, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2855, pbuffer, 3491, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3065, pbuffer, 3911, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3345, pbuffer, 4471, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3495, pbuffer, 4621, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3720, pbuffer, 4846, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4035, pbuffer, 5161, 420, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 138, cbuffer, 0, 10, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 987, cbuffer, 25, 55, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 3120, cbuffer, 100, 160, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5868, cbuffer, 990, 1000, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 6717, cbuffer, 1015, 1045, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 8850, cbuffer, 1090, 1150, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 12840, cbuffer, 1240, 1340, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 19260, cbuffer, 1490, 1640, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 0, cbuffer, 250, 260, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 30, cbuffer, 260, 275, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 75, cbuffer, 275, 296, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 168, cbuffer, 0, 0, 30, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 258, cbuffer, 10, 30, 75, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 393, 138, 168, 258, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 573, cbuffer, 324, 354, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 663, cbuffer, 354, 399, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 798, cbuffer, 399, 462, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 1077, cbuffer, 25, 573, 663, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 1347, cbuffer, 55, 663, 798, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 1752, 987, 1077, 1347, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 2292, cbuffer, 546, 606, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 2472, cbuffer, 606, 696, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 2742, cbuffer, 696, 822, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 3300, cbuffer, 100, 2292, 2472, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 3840, cbuffer, 160, 2472, 2742, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 4650, 3120, 3300, 3840, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 5730, cbuffer, 1865, 1875, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 5760, cbuffer, 1875, 1890, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 5805, cbuffer, 1890, 1911, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 5898, cbuffer, 990, 5730, 5760, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 5988, cbuffer, 1000, 5760, 5805, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 6123, 5868, 5898, 5988, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6303, cbuffer, 1939, 1969, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 6393, cbuffer, 1969, 2014, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 6528, cbuffer, 2014, 2077, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 6807, cbuffer, 1015, 6303, 6393, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 7077, cbuffer, 1045, 6393, 6528, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 7482, 6717, 6807, 7077, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 8022, cbuffer, 2161, 2221, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 8202, cbuffer, 2221, 2311, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 8472, cbuffer, 2311, 2437, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 9030, cbuffer, 1090, 8022, 8202, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 9570, cbuffer, 1150, 8202, 8472, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 10380, 8850, 9030, 9570, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 11460, cbuffer, 2605, 2705, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 11760, cbuffer, 2705, 2855, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 12210, cbuffer, 2855, 3065, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 13140, cbuffer, 1240, 11460, 11760, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 14040, cbuffer, 1340, 11760, 12210, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 15390, 12840, 13140, 14040, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 17190, cbuffer, 3345, 3495, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 17640, cbuffer, 3495, 3720, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 18315, cbuffer, 3720, 4035, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 19710, cbuffer, 1490, 17190, 17640, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 21060, cbuffer, 1640, 17640, 18315, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 23085, 19260, 19710, 21060, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<2, 3>(skbuffer, 0, ckbuffer, 393, 0, 0);

            t4cfunc::ket_transform<2, 3>(skbuffer, 35, ckbuffer, 453, 0, 0);

            t4cfunc::ket_transform<2, 3>(skbuffer, 70, ckbuffer, 513, 0, 0);

            t4cfunc::ket_transform<2, 3>(skbuffer, 420, ckbuffer, 1752, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 525, ckbuffer, 1932, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 630, ckbuffer, 2112, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 1680, ckbuffer, 4650, 0, 2);

            t4cfunc::ket_transform<2, 3>(skbuffer, 1890, ckbuffer, 5010, 0, 2);

            t4cfunc::ket_transform<2, 3>(skbuffer, 2100, ckbuffer, 5370, 0, 2);

            t4cfunc::ket_transform<2, 3>(skbuffer, 29400, ckbuffer, 6123, 0, 0);

            t4cfunc::ket_transform<2, 3>(skbuffer, 29435, ckbuffer, 6183, 0, 0);

            t4cfunc::ket_transform<2, 3>(skbuffer, 29470, ckbuffer, 6243, 0, 0);

            t4cfunc::ket_transform<2, 3>(skbuffer, 29505, ckbuffer, 7482, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 29610, ckbuffer, 7662, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 29715, ckbuffer, 7842, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 29820, ckbuffer, 10380, 0, 2);

            t4cfunc::ket_transform<2, 3>(skbuffer, 30030, ckbuffer, 10740, 0, 2);

            t4cfunc::ket_transform<2, 3>(skbuffer, 30240, ckbuffer, 11100, 0, 2);

            t4cfunc::ket_transform<2, 3>(skbuffer, 30450, ckbuffer, 15390, 0, 3);

            t4cfunc::ket_transform<2, 3>(skbuffer, 30800, ckbuffer, 15990, 0, 3);

            t4cfunc::ket_transform<2, 3>(skbuffer, 31150, ckbuffer, 16590, 0, 3);

            t4cfunc::ket_transform<2, 3>(skbuffer, 31500, ckbuffer, 23085, 0, 4);

            t4cfunc::ket_transform<2, 3>(skbuffer, 32025, ckbuffer, 23985, 0, 4);

            t4cfunc::ket_transform<2, 3>(skbuffer, 32550, ckbuffer, 24885, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 7350, 0, 420, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 7455, 35, 525, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 7560, 70, 630, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 8610, 420, 1680, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 8925, 525, 1890, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 9240, 630, 2100, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 18060, 7350, 8610, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 18270, 7455, 8925, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 18480, 7560, 9240, r_ab, 2, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 105, 29400, 29505, r_ab, 2, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 735, 29505, 29820, r_ab, 2, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 2310, 29820, 30450, r_ab, 2, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 4200, 30450, 31500, r_ab, 2, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 7665, 0, 105, 735, r_ab, 2, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 9555, 420, 735, 2310, r_ab, 2, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 12390, 1680, 2310, 4200, r_ab, 2, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 18690, 7350, 7665, 9555, r_ab, 2, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 20580, 8610, 9555, 12390, r_ab, 2, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fsxx(skbuffer, 26250, 18060, 18690, 20580, r_ab, 2, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 26250, 2, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 245, skbuffer, 26600, 2, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 490, skbuffer, 26950, 2, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 735, skbuffer, 27300, 2, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 980, skbuffer, 27650, 2, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1225, skbuffer, 28000, 2, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1470, skbuffer, 28350, 2, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1715, skbuffer, 28700, 2, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1960, skbuffer, 29050, 2, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 2, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFSDF_hpp */
