#ifndef ElectronRepulsionGeom2000RecDSDF_hpp
#define ElectronRepulsionGeom2000RecDSDF_hpp

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
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionGeom1000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
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

/// @brief Computes d^(2)/dA^(2)(DS|1/|r-r'||DF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_dsdf(T& distributor,
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

    CSimdArray<double> pbuffer(3565, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2116, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(6210, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(12425, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1050, 1);

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

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 245, 85, 95, 155, 170, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 266, 95, 105, 170, 185, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 287, 105, 115, 185, 200, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 308, 115, 125, 200, 215, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 329, 125, 135, 215, 230, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 350, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 353, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 362, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 371, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 389, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 407, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 425, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 455, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 485, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 515, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 545, 95, 155, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 590, 105, 170, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 635, 115, 185, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 680, 125, 200, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 725, 170, 245, 266, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 788, 185, 266, 287, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 851, 200, 287, 308, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 914, 215, 308, 329, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 977, 16, 19, 350, 353, 362, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 995, 43, 49, 353, 371, 389, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1031, 49, 55, 362, 389, 407, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1067, 85, 95, 371, 425, 455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1127, 95, 105, 389, 455, 485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1187, 105, 115, 407, 485, 515, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1247, 155, 170, 455, 545, 590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1337, 170, 185, 485, 590, 635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1427, 185, 200, 515, 635, 680, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1517, 245, 266, 590, 725, 788, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1643, 266, 287, 635, 788, 851, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1769, 287, 308, 680, 851, 914, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1895, 371, 389, 977, 995, 1031, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1955, 425, 455, 995, 1067, 1127, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2055, 455, 485, 1031, 1127, 1187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2155, 545, 590, 1127, 1247, 1337, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2305, 590, 635, 1187, 1337, 1427, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 2455, 725, 788, 1337, 1517, 1643, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 2665, 788, 851, 1427, 1643, 1769, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2875, 1067, 1127, 1895, 1955, 2055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3025, 1247, 1337, 2055, 2155, 2305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 3250, 1517, 1643, 2305, 2455, 2665, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 25, pbuffer, 245, 21, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {85, 95});

                pbuffer.scale(2.0 * a_exp, {155, 170});

                pbuffer.scale(2.0 * a_exp, {245, 266});

                pbuffer.scale(2.0 * a_exp, {425, 455});

                pbuffer.scale(2.0 * a_exp, {545, 590});

                pbuffer.scale(2.0 * a_exp, {725, 788});

                pbuffer.scale(2.0 * a_exp, {1067, 1127});

                pbuffer.scale(2.0 * a_exp, {1247, 1337});

                pbuffer.scale(2.0 * a_exp, {1517, 1643});

                t2cfunc::reduce(cbuffer, 46, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 56, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 71, pbuffer, 245, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 92, pbuffer, 425, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 122, pbuffer, 545, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 167, pbuffer, 725, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 230, pbuffer, 1067, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 290, pbuffer, 1247, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 380, pbuffer, 1517, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {85, 95});

                pbuffer.scale(2.0 * a_exp, {155, 170});

                pbuffer.scale(2.0 * a_exp, {245, 266});

                pbuffer.scale(2.0 * a_exp, {425, 455});

                pbuffer.scale(2.0 * a_exp, {545, 590});

                pbuffer.scale(2.0 * a_exp, {725, 788});

                pbuffer.scale(2.0 * a_exp, {1067, 1127});

                pbuffer.scale(2.0 * a_exp, {1247, 1337});

                pbuffer.scale(2.0 * a_exp, {1517, 1643});

                pbuffer.scale(4.0 * a_exp * a_exp, {1955, 2055});

                pbuffer.scale(4.0 * a_exp * a_exp, {2155, 2305});

                pbuffer.scale(4.0 * a_exp * a_exp, {2455, 2665});

                pbuffer.scale(4.0 * a_exp * a_exp, {2875, 3025});

                pbuffer.scale(4.0 * a_exp * a_exp, {3025, 3250});

                pbuffer.scale(4.0 * a_exp * a_exp, {3250, 3565});

                t2cfunc::reduce(cbuffer, 506, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 516, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 531, pbuffer, 245, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 552, pbuffer, 425, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 582, pbuffer, 545, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 627, pbuffer, 725, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 690, pbuffer, 1067, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 750, pbuffer, 1247, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 840, pbuffer, 1517, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 966, pbuffer, 1955, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1066, pbuffer, 2155, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1216, pbuffer, 2455, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1426, pbuffer, 2875, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1576, pbuffer, 3025, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1801, pbuffer, 3250, 315, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 10, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 30, cbuffer, 10, 25, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 75, 0, 30, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 135, cbuffer, 46, 56, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 165, cbuffer, 56, 71, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 210, 135, 165, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 270, cbuffer, 92, 122, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 360, cbuffer, 122, 167, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 495, 270, 360, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 675, cbuffer, 230, 290, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 855, cbuffer, 290, 380, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 1125, 675, 855, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1485, cbuffer, 506, 516, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 1515, cbuffer, 516, 531, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 1560, 1485, 1515, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1620, cbuffer, 552, 582, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 1710, cbuffer, 582, 627, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 1845, 1620, 1710, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2025, cbuffer, 690, 750, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 2205, cbuffer, 750, 840, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 2475, 2025, 2205, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2835, cbuffer, 966, 1066, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 3135, cbuffer, 1066, 1216, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 3585, 2835, 3135, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4185, cbuffer, 1426, 1576, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 4635, cbuffer, 1576, 1801, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 5310, 4185, 4635, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<2, 3>(skbuffer, 0, ckbuffer, 75, 0, 0);

            t4cfunc::ket_transform<2, 3>(skbuffer, 6230, ckbuffer, 210, 0, 0);

            t4cfunc::ket_transform<2, 3>(skbuffer, 6265, ckbuffer, 495, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 6370, ckbuffer, 1125, 0, 2);

            t4cfunc::ket_transform<2, 3>(skbuffer, 7000, ckbuffer, 1560, 0, 0);

            t4cfunc::ket_transform<2, 3>(skbuffer, 7035, ckbuffer, 1845, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 7140, ckbuffer, 2475, 0, 2);

            t4cfunc::ket_transform<2, 3>(skbuffer, 7350, ckbuffer, 3585, 0, 3);

            t4cfunc::ket_transform<2, 3>(skbuffer, 7700, ckbuffer, 5310, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 6580, 6230, 6265, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 6685, 6265, 6370, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 8225, 7000, 7035, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 8330, 7035, 7140, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 8645, 7140, 7350, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 9275, 7350, 7700, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 10325, 8225, 8330, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 10535, 8330, 8645, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 11165, 8645, 9275, r_ab, 2, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_psxx(skbuffer, 2135, 0, 6580, 6685, r_ab, 2, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 35, 6230, 10325, 0, 2, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 245, 6265, 10535, 1, 2, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 875, 6370, 11165, 2, 2, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_psxx(skbuffer, 2450, 6580, 35, 245, r_ab, 2, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_ppxx(skbuffer, 3080, 6685, 245, 875, r_ab, 2, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dsxx(skbuffer, 4970, 2135, 2450, 3080, r_ab, 2, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 0, skbuffer, 4970, 2, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 175, skbuffer, 5180, 2, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 350, skbuffer, 5390, 2, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 525, skbuffer, 5600, 2, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 700, skbuffer, 5810, 2, 3);

            t4cfunc::bra_transform<2, 0>(sbuffer, 875, skbuffer, 6020, 2, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 0, 2, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecDSDF_hpp */
