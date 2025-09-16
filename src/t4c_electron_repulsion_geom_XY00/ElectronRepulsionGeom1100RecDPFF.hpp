#ifndef ElectronRepulsionGeom1100RecDPFF_hpp
#define ElectronRepulsionGeom1100RecDPFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecXXDF.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXFF.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
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
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSI.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DP|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dpff(T& distributor,
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

    CSimdArray<double> pbuffer(10005, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(6808, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(48596, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(43022, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(6615, 1);

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

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 462, 195, 210, 315, 336, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 490, 210, 225, 336, 357, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 518, 225, 240, 357, 378, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 546, 240, 255, 378, 399, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 574, 255, 270, 399, 420, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 602, 270, 285, 420, 441, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 630, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 633, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 636, 3, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 645, 4, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 654, 5, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 663, 18, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 681, 21, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 699, 24, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 717, 27, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 735, 51, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 765, 57, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 795, 63, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 825, 69, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 855, 75, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 885, 115, 195, 210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 930, 125, 210, 225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 975, 135, 225, 240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1020, 145, 240, 255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1065, 155, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1110, 210, 315, 336, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1173, 225, 336, 357, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1236, 240, 357, 378, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1299, 255, 378, 399, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1362, 270, 399, 420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1425, 336, 462, 490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1509, 357, 490, 518, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1593, 378, 518, 546, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1677, 399, 546, 574, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1761, 420, 574, 602, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1845, 3, 4, 630, 633, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1851, 18, 21, 630, 636, 645, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1869, 21, 24, 633, 645, 654, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1887, 51, 57, 636, 663, 681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1923, 57, 63, 645, 681, 699, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1959, 63, 69, 654, 699, 717, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1995, 105, 115, 663, 735, 765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2055, 115, 125, 681, 765, 795, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2115, 125, 135, 699, 795, 825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2175, 135, 145, 717, 825, 855, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2235, 195, 210, 765, 885, 930, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2325, 210, 225, 795, 930, 975, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2415, 225, 240, 825, 975, 1020, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2505, 240, 255, 855, 1020, 1065, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2595, 315, 336, 930, 1110, 1173, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2721, 336, 357, 975, 1173, 1236, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2847, 357, 378, 1020, 1236, 1299, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2973, 378, 399, 1065, 1299, 1362, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3099, 462, 490, 1173, 1425, 1509, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3267, 490, 518, 1236, 1509, 1593, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3435, 518, 546, 1299, 1593, 1677, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3603, 546, 574, 1362, 1677, 1761, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3771, 636, 645, 1845, 1851, 1869, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3801, 663, 681, 1851, 1887, 1923, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3861, 681, 699, 1869, 1923, 1959, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3921, 735, 765, 1887, 1995, 2055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4021, 765, 795, 1923, 2055, 2115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4121, 795, 825, 1959, 2115, 2175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4221, 885, 930, 2055, 2235, 2325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4371, 930, 975, 2115, 2325, 2415, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4521, 975, 1020, 2175, 2415, 2505, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4671, 1110, 1173, 2325, 2595, 2721, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4881, 1173, 1236, 2415, 2721, 2847, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5091, 1236, 1299, 2505, 2847, 2973, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 5301, 1425, 1509, 2721, 3099, 3267, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 5581, 1509, 1593, 2847, 3267, 3435, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 5861, 1593, 1677, 2973, 3435, 3603, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 6141, 1887, 1923, 3771, 3801, 3861, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 6231, 1995, 2055, 3801, 3921, 4021, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 6381, 2055, 2115, 3861, 4021, 4121, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6531, 2235, 2325, 4021, 4221, 4371, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6756, 2325, 2415, 4121, 4371, 4521, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 6981, 2595, 2721, 4371, 4671, 4881, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 7296, 2721, 2847, 4521, 4881, 5091, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 7611, 3099, 3267, 4881, 5301, 5581, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 8031, 3267, 3435, 5091, 5581, 5861, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 8451, 3921, 4021, 6141, 6231, 6381, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 8661, 4221, 4371, 6381, 6531, 6756, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 8976, 4671, 4881, 6756, 6981, 7296, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 9417, 5301, 5581, 7296, 7611, 8031, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 105, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 195, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 25, pbuffer, 315, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 46, pbuffer, 462, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 74, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 104, pbuffer, 885, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 149, pbuffer, 1110, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 212, pbuffer, 1425, 84, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {1995, 2055});

                pbuffer.scale(2.0 * b_exp, {2235, 2325});

                pbuffer.scale(2.0 * b_exp, {2595, 2721});

                pbuffer.scale(2.0 * b_exp, {3099, 3267});

                pbuffer.scale(2.0 * b_exp, {3921, 4021});

                pbuffer.scale(2.0 * b_exp, {4221, 4371});

                pbuffer.scale(2.0 * b_exp, {4671, 4881});

                pbuffer.scale(2.0 * b_exp, {5301, 5581});

                t2cfunc::reduce(cbuffer, 296, pbuffer, 1995, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 356, pbuffer, 2235, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 446, pbuffer, 2595, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 572, pbuffer, 3099, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 740, pbuffer, 3921, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 840, pbuffer, 4221, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 990, pbuffer, 4671, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1200, pbuffer, 5301, 280, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {105, 115});

                pbuffer.scale(2.0 * a_exp, {195, 210});

                pbuffer.scale(2.0 * a_exp, {315, 336});

                pbuffer.scale(2.0 * a_exp, {462, 490});

                pbuffer.scale(2.0 * a_exp, {735, 765});

                pbuffer.scale(2.0 * a_exp, {885, 930});

                pbuffer.scale(2.0 * a_exp, {1110, 1173});

                pbuffer.scale(2.0 * a_exp, {1425, 1509});

                pbuffer.scale(a_exp / b_exp, {1995, 2055});

                pbuffer.scale(a_exp / b_exp, {2235, 2325});

                pbuffer.scale(a_exp / b_exp, {2595, 2721});

                pbuffer.scale(a_exp / b_exp, {3099, 3267});

                pbuffer.scale(a_exp / b_exp, {3921, 4021});

                pbuffer.scale(a_exp / b_exp, {4221, 4371});

                pbuffer.scale(a_exp / b_exp, {4671, 4881});

                pbuffer.scale(a_exp / b_exp, {5301, 5581});

                t2cfunc::reduce(cbuffer, 1480, pbuffer, 105, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1490, pbuffer, 195, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1505, pbuffer, 315, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1526, pbuffer, 462, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1554, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1584, pbuffer, 885, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1629, pbuffer, 1110, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1692, pbuffer, 1425, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1776, pbuffer, 1995, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1836, pbuffer, 2235, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1926, pbuffer, 2595, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2052, pbuffer, 3099, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2220, pbuffer, 3921, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2320, pbuffer, 4221, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2470, pbuffer, 4671, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2680, pbuffer, 5301, 280, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {1995, 2055});

                pbuffer.scale(2.0 * b_exp, {2235, 2325});

                pbuffer.scale(2.0 * b_exp, {2595, 2721});

                pbuffer.scale(2.0 * b_exp, {3099, 3267});

                pbuffer.scale(2.0 * b_exp, {3921, 4021});

                pbuffer.scale(2.0 * b_exp, {4221, 4371});

                pbuffer.scale(2.0 * b_exp, {4671, 4881});

                pbuffer.scale(2.0 * b_exp, {5301, 5581});

                pbuffer.scale(4.0 * a_exp * b_exp, {6231, 6381});

                pbuffer.scale(4.0 * a_exp * b_exp, {6531, 6756});

                pbuffer.scale(4.0 * a_exp * b_exp, {6981, 7296});

                pbuffer.scale(4.0 * a_exp * b_exp, {7611, 8031});

                pbuffer.scale(4.0 * a_exp * b_exp, {8451, 8661});

                pbuffer.scale(4.0 * a_exp * b_exp, {8661, 8976});

                pbuffer.scale(4.0 * a_exp * b_exp, {8976, 9417});

                pbuffer.scale(4.0 * a_exp * b_exp, {9417, 10005});

                t2cfunc::reduce(cbuffer, 2960, pbuffer, 1995, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3020, pbuffer, 2235, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3110, pbuffer, 2595, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3236, pbuffer, 3099, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3404, pbuffer, 3921, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3504, pbuffer, 4221, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3654, pbuffer, 4671, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3864, pbuffer, 5301, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4144, pbuffer, 6231, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4294, pbuffer, 6531, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4519, pbuffer, 6981, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4834, pbuffer, 7611, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5254, pbuffer, 8451, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5464, pbuffer, 8661, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5779, pbuffer, 8976, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6220, pbuffer, 9417, 588, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 10, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 30, cbuffer, 10, 25, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 75, cbuffer, 25, 46, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 138, 0, 30, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 198, 30, 75, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 288, 138, 198, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 388, cbuffer, 74, 104, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 478, cbuffer, 104, 149, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 613, cbuffer, 149, 212, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 802, 388, 478, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 982, 478, 613, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 1252, 802, 982, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4252, cbuffer, 296, 356, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 4432, cbuffer, 356, 446, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 4702, cbuffer, 446, 572, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 5080, 4252, 4432, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 5440, 4432, 4702, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 5980, 5080, 5440, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 6580, cbuffer, 740, 840, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 6880, cbuffer, 840, 990, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 7330, cbuffer, 990, 1200, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 7960, 6580, 6880, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 8560, 6880, 7330, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 9460, 7960, 8560, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 10460, cbuffer, 1480, 1490, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 10490, cbuffer, 1490, 1505, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 10535, cbuffer, 1505, 1526, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 10598, 10460, 10490, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 10658, 10490, 10535, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 10748, 10598, 10658, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 10848, cbuffer, 1554, 1584, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 10938, cbuffer, 1584, 1629, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 11073, cbuffer, 1629, 1692, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 11262, 10848, 10938, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 11442, 10938, 11073, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 11712, 11262, 11442, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 12912, cbuffer, 1776, 1836, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 13092, cbuffer, 1836, 1926, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 13362, cbuffer, 1926, 2052, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 13740, 12912, 13092, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 14100, 13092, 13362, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 14640, 13740, 14100, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 17040, cbuffer, 2220, 2320, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 17340, cbuffer, 2320, 2470, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 17790, cbuffer, 2470, 2680, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 18420, 17040, 17340, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 19020, 17340, 17790, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 19920, 18420, 19020, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 28420, cbuffer, 2960, 3020, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 28600, cbuffer, 3020, 3110, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 28870, cbuffer, 3110, 3236, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 29248, 28420, 28600, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 29608, 28600, 28870, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 30148, 29248, 29608, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 30748, cbuffer, 3404, 3504, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 31048, cbuffer, 3504, 3654, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 31498, cbuffer, 3654, 3864, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 32128, 30748, 31048, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 32728, 31048, 31498, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 33628, 32128, 32728, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 34628, cbuffer, 4144, 4294, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 35078, cbuffer, 4294, 4519, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 35753, cbuffer, 4519, 4834, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 36698, 34628, 35078, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 37598, 35078, 35753, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 38948, 36698, 37598, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 40448, cbuffer, 5254, 5464, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 41078, cbuffer, 5464, 5779, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 42023, cbuffer, 5779, 6220, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 43346, 40448, 41078, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 44606, 41078, 42023, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 46496, 43346, 44606, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 288, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 49, ckbuffer, 1252, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 32389, ckbuffer, 5980, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 32683, ckbuffer, 9460, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 33173, ckbuffer, 10748, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 33222, ckbuffer, 11712, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 33810, ckbuffer, 14640, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 34986, ckbuffer, 19920, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 40474, ckbuffer, 30148, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 40768, ckbuffer, 33628, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 41258, ckbuffer, 38948, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 41993, ckbuffer, 46496, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 39151, 33222, 33810, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 39592, 33810, 34986, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 11221, 49, 39151, 39592, r_ab, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 196, 0, 32389, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 1960, 49, 32683, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 9898, 49, 196, 1960, r_ab, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 33369, 33173, 40474, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 34104, 33222, 40768, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 35476, 33810, 41258, 3, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 36946, 34986, 41993, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 637, 33222, 33369, 34104, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 2842, 33810, 34104, 35476, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 5488, 34986, 35476, 36946, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 12544, 196, 39151, 637, 2842, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 16513, 1960, 39592, 2842, 5488, r_ab, 3, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 24451, 9898, 11221, 12544, 16513, r_ab, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 24451, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 735, skbuffer, 25333, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1470, skbuffer, 26215, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2205, skbuffer, 27097, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2940, skbuffer, 27979, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 3675, skbuffer, 28861, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 4410, skbuffer, 29743, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 5145, skbuffer, 30625, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 5880, skbuffer, 31507, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDPFF_hpp */
