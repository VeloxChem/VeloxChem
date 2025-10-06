#ifndef ElectronRepulsionGeom1100RecDDDF_hpp
#define ElectronRepulsionGeom1100RecDDDF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecXXDF.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionGeom0100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
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
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DD|1/|r-r'||DF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dddf(T& distributor,
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

    CSimdArray<double> pbuffer(10746, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(6532, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(31410, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(53795, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(7875, 1);

    // setup Boys fuction data

    const CBoysFunc<11> bf_table;

    CSimdArray<double> bf_data(13, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 12, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 12, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 12, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 12);
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

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 11, pfactors, 16, bf_data, 11);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 12, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 21, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 24, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 27, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 30, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 33, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 36, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 39, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 42, 10, 11, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 45, 0, 1, 12, 15, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 51, 1, 2, 15, 18, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 57, 2, 3, 18, 21, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 63, 3, 4, 21, 24, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 69, 4, 5, 24, 27, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 75, 5, 6, 27, 30, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 81, 6, 7, 30, 33, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 87, 7, 8, 33, 36, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 93, 8, 9, 36, 39, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 99, 9, 10, 39, 42, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 12, 15, 45, 51, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 15, 18, 51, 57, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 18, 21, 57, 63, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 21, 24, 63, 69, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 24, 27, 69, 75, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 27, 30, 75, 81, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 30, 33, 81, 87, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 175, 33, 36, 87, 93, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 185, 36, 39, 93, 99, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 195, 45, 51, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 210, 51, 57, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 225, 57, 63, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 240, 63, 69, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 255, 69, 75, 145, 155, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 270, 75, 81, 155, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 285, 81, 87, 165, 175, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 300, 87, 93, 175, 185, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 315, 105, 115, 195, 210, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 336, 115, 125, 210, 225, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 357, 125, 135, 225, 240, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 378, 135, 145, 240, 255, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 399, 145, 155, 255, 270, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 420, 155, 165, 270, 285, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 441, 165, 175, 285, 300, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 462, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 465, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 468, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 471, 3, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 480, 4, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 489, 5, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 498, 6, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 507, 18, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 525, 21, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 543, 24, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 561, 27, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 579, 30, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 597, 51, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 627, 57, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 657, 63, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 687, 69, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 717, 75, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 747, 81, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 777, 115, 195, 210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 822, 125, 210, 225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 867, 135, 225, 240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 912, 145, 240, 255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 957, 155, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1002, 165, 270, 285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1047, 210, 315, 336, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1110, 225, 336, 357, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1173, 240, 357, 378, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1236, 255, 378, 399, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1299, 270, 399, 420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1362, 285, 420, 441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1425, 3, 4, 462, 465, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1431, 4, 5, 465, 468, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1437, 18, 21, 462, 471, 480, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1455, 21, 24, 465, 480, 489, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1473, 24, 27, 468, 489, 498, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1491, 51, 57, 471, 507, 525, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1527, 57, 63, 480, 525, 543, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1563, 63, 69, 489, 543, 561, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1599, 69, 75, 498, 561, 579, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1635, 105, 115, 507, 597, 627, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1695, 115, 125, 525, 627, 657, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1755, 125, 135, 543, 657, 687, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1815, 135, 145, 561, 687, 717, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1875, 145, 155, 579, 717, 747, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1935, 195, 210, 627, 777, 822, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2025, 210, 225, 657, 822, 867, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2115, 225, 240, 687, 867, 912, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2205, 240, 255, 717, 912, 957, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2295, 255, 270, 747, 957, 1002, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2385, 315, 336, 822, 1047, 1110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2511, 336, 357, 867, 1110, 1173, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2637, 357, 378, 912, 1173, 1236, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2763, 378, 399, 957, 1236, 1299, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2889, 399, 420, 1002, 1299, 1362, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3015, 462, 465, 1425, 1431, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3025, 471, 480, 1425, 1437, 1455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3055, 480, 489, 1431, 1455, 1473, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3085, 507, 525, 1437, 1491, 1527, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3145, 525, 543, 1455, 1527, 1563, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3205, 543, 561, 1473, 1563, 1599, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3265, 597, 627, 1491, 1635, 1695, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3365, 627, 657, 1527, 1695, 1755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3465, 657, 687, 1563, 1755, 1815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3565, 687, 717, 1599, 1815, 1875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3665, 777, 822, 1695, 1935, 2025, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3815, 822, 867, 1755, 2025, 2115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3965, 867, 912, 1815, 2115, 2205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4115, 912, 957, 1875, 2205, 2295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4265, 1047, 1110, 2025, 2385, 2511, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4475, 1110, 1173, 2115, 2511, 2637, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4685, 1173, 1236, 2205, 2637, 2763, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4895, 1236, 1299, 2295, 2763, 2889, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 5105, 1437, 1455, 3015, 3025, 3055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 5150, 1491, 1527, 3025, 3085, 3145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 5240, 1527, 1563, 3055, 3145, 3205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5330, 1635, 1695, 3085, 3265, 3365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5480, 1695, 1755, 3145, 3365, 3465, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5630, 1755, 1815, 3205, 3465, 3565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5780, 1935, 2025, 3365, 3665, 3815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6005, 2025, 2115, 3465, 3815, 3965, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6230, 2115, 2205, 3565, 3965, 4115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 6455, 2385, 2511, 3815, 4265, 4475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 6770, 2511, 2637, 3965, 4475, 4685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 7085, 2637, 2763, 4115, 4685, 4895, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 7400, 3085, 3145, 5105, 5150, 5240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 7526, 3265, 3365, 5150, 5330, 5480, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 7736, 3365, 3465, 5240, 5480, 5630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 7946, 3665, 3815, 5480, 5780, 6005, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 8261, 3815, 3965, 5630, 6005, 6230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 8576, 4265, 4475, 6005, 6455, 6770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 9017, 4475, 4685, 6230, 6770, 7085, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 9458, 5330, 5480, 7400, 7526, 7736, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 9738, 5780, 6005, 7736, 7946, 8261, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 10158, 6455, 6770, 8261, 8576, 9017, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 597, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 777, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 75, pbuffer, 1047, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 138, pbuffer, 1635, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 198, pbuffer, 1935, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 288, pbuffer, 2385, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {3265, 3365});

                pbuffer.scale(2.0 * b_exp, {3665, 3815});

                pbuffer.scale(2.0 * b_exp, {4265, 4475});

                pbuffer.scale(2.0 * b_exp, {5330, 5480});

                pbuffer.scale(2.0 * b_exp, {5780, 6005});

                pbuffer.scale(2.0 * b_exp, {6455, 6770});

                t2cfunc::reduce(cbuffer, 414, pbuffer, 3265, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 514, pbuffer, 3665, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 664, pbuffer, 4265, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 874, pbuffer, 5330, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1024, pbuffer, 5780, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1249, pbuffer, 6455, 315, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {597, 627});

                pbuffer.scale(2.0 * a_exp, {777, 822});

                pbuffer.scale(2.0 * a_exp, {1047, 1110});

                pbuffer.scale(2.0 * a_exp, {1635, 1695});

                pbuffer.scale(2.0 * a_exp, {1935, 2025});

                pbuffer.scale(2.0 * a_exp, {2385, 2511});

                pbuffer.scale(a_exp / b_exp, {3265, 3365});

                pbuffer.scale(a_exp / b_exp, {3665, 3815});

                pbuffer.scale(a_exp / b_exp, {4265, 4475});

                pbuffer.scale(a_exp / b_exp, {5330, 5480});

                pbuffer.scale(a_exp / b_exp, {5780, 6005});

                pbuffer.scale(a_exp / b_exp, {6455, 6770});

                t2cfunc::reduce(cbuffer, 1564, pbuffer, 597, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1594, pbuffer, 777, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1639, pbuffer, 1047, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1702, pbuffer, 1635, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1762, pbuffer, 1935, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1852, pbuffer, 2385, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1978, pbuffer, 3265, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2078, pbuffer, 3665, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2228, pbuffer, 4265, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2438, pbuffer, 5330, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2588, pbuffer, 5780, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2813, pbuffer, 6455, 315, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {3265, 3365});

                pbuffer.scale(2.0 * b_exp, {3665, 3815});

                pbuffer.scale(2.0 * b_exp, {4265, 4475});

                pbuffer.scale(2.0 * b_exp, {5330, 5480});

                pbuffer.scale(2.0 * b_exp, {5780, 6005});

                pbuffer.scale(2.0 * b_exp, {6455, 6770});

                pbuffer.scale(4.0 * a_exp * b_exp, {7526, 7736});

                pbuffer.scale(4.0 * a_exp * b_exp, {7946, 8261});

                pbuffer.scale(4.0 * a_exp * b_exp, {8576, 9017});

                pbuffer.scale(4.0 * a_exp * b_exp, {9458, 9738});

                pbuffer.scale(4.0 * a_exp * b_exp, {9738, 10158});

                pbuffer.scale(4.0 * a_exp * b_exp, {10158, 10746});

                t2cfunc::reduce(cbuffer, 3128, pbuffer, 3265, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3228, pbuffer, 3665, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3378, pbuffer, 4265, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3588, pbuffer, 5330, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3738, pbuffer, 5780, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3963, pbuffer, 6455, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4278, pbuffer, 7526, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4488, pbuffer, 7946, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4803, pbuffer, 8576, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5244, pbuffer, 9458, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5524, pbuffer, 9738, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5944, pbuffer, 10158, 588, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 30, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 90, cbuffer, 30, 75, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 225, 0, 90, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 405, cbuffer, 138, 198, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 585, cbuffer, 198, 288, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 855, 405, 585, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4095, cbuffer, 414, 514, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 4395, cbuffer, 514, 664, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 4845, 4095, 4395, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5445, cbuffer, 874, 1024, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 5895, cbuffer, 1024, 1249, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 6570, 5445, 5895, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 7470, cbuffer, 1564, 1594, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 7560, cbuffer, 1594, 1639, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 7695, 7470, 7560, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 7875, cbuffer, 1702, 1762, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 8055, cbuffer, 1762, 1852, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 8325, 7875, 8055, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9765, cbuffer, 1978, 2078, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 10065, cbuffer, 2078, 2228, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 10515, 9765, 10065, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 12915, cbuffer, 2438, 2588, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 13365, cbuffer, 2588, 2813, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 14040, 12915, 13365, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 21420, cbuffer, 3128, 3228, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 21720, cbuffer, 3228, 3378, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 22170, 21420, 21720, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 22770, cbuffer, 3588, 3738, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 23220, cbuffer, 3738, 3963, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 23895, 22770, 23220, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 24795, cbuffer, 4278, 4488, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 25425, cbuffer, 4488, 4803, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 26370, 24795, 25425, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 27630, cbuffer, 5244, 5524, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 28470, cbuffer, 5524, 5944, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 29730, 27630, 28470, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<2, 3>(skbuffer, 0, ckbuffer, 225, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 105, ckbuffer, 855, 0, 2);

            t4cfunc::ket_transform<2, 3>(skbuffer, 42000, ckbuffer, 4845, 0, 3);

            t4cfunc::ket_transform<2, 3>(skbuffer, 42350, ckbuffer, 6570, 0, 4);

            t4cfunc::ket_transform<2, 3>(skbuffer, 42875, ckbuffer, 7695, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 42980, ckbuffer, 8325, 0, 2);

            t4cfunc::ket_transform<2, 3>(skbuffer, 43820, ckbuffer, 10515, 0, 3);

            t4cfunc::ket_transform<2, 3>(skbuffer, 45220, ckbuffer, 14040, 0, 4);

            t4cfunc::ket_transform<2, 3>(skbuffer, 51205, ckbuffer, 22170, 0, 3);

            t4cfunc::ket_transform<2, 3>(skbuffer, 51555, ckbuffer, 23895, 0, 4);

            t4cfunc::ket_transform<2, 3>(skbuffer, 52080, ckbuffer, 26370, 0, 5);

            t4cfunc::ket_transform<2, 3>(skbuffer, 52815, ckbuffer, 29730, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 49525, 42980, 43820, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 50155, 43820, 45220, r_ab, 2, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 13650, 105, 49525, 50155, r_ab, 2, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 315, 0, 42000, 2, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 2835, 105, 42350, 2, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 11760, 105, 315, 2835, r_ab, 2, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 43190, 42875, 51205, 2, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 44170, 42980, 51555, 2, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 45745, 43820, 52080, 2, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 47320, 45220, 52815, 2, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 945, 42980, 43190, 44170, r_ab, 2, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 3885, 43820, 44170, 45745, r_ab, 2, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 7035, 45220, 45745, 47320, r_ab, 2, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 15540, 315, 49525, 945, 3885, r_ab, 2, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 21210, 2835, 50155, 3885, 7035, r_ab, 2, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 30660, 11760, 13650, 15540, 21210, r_ab, 2, 3);

            t4cfunc::bra_transform<2, 2>(sbuffer, 0, skbuffer, 30660, 2, 3);

            t4cfunc::bra_transform<2, 2>(sbuffer, 875, skbuffer, 31920, 2, 3);

            t4cfunc::bra_transform<2, 2>(sbuffer, 1750, skbuffer, 33180, 2, 3);

            t4cfunc::bra_transform<2, 2>(sbuffer, 2625, skbuffer, 34440, 2, 3);

            t4cfunc::bra_transform<2, 2>(sbuffer, 3500, skbuffer, 35700, 2, 3);

            t4cfunc::bra_transform<2, 2>(sbuffer, 4375, skbuffer, 36960, 2, 3);

            t4cfunc::bra_transform<2, 2>(sbuffer, 5250, skbuffer, 38220, 2, 3);

            t4cfunc::bra_transform<2, 2>(sbuffer, 6125, skbuffer, 39480, 2, 3);

            t4cfunc::bra_transform<2, 2>(sbuffer, 7000, skbuffer, 40740, 2, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 2, 2, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDDDF_hpp */
