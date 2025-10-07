#ifndef ElectronRepulsionGeom1010RecDPFF_hpp
#define ElectronRepulsionGeom1010RecDPFF_hpp

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
#include "ElectronRepulsionPrimRecSDSK.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSI.hpp"
#include "ElectronRepulsionPrimRecSFSK.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSI.hpp"
#include "ElectronRepulsionPrimRecSGSK.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
#include "ElectronRepulsionPrimRecSPSI.hpp"
#include "ElectronRepulsionPrimRecSPSK.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DP|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dpff(T& distributor,
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

    CSimdArray<double> pbuffer(8185, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(6708, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(65403, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(35868, 1);

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

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 630, 315, 336, 462, 490, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 666, 336, 357, 490, 518, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 702, 357, 378, 518, 546, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 738, 378, 399, 546, 574, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssk(pbuffer, 774, 399, 420, 574, 602, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 810, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 813, 3, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 822, 4, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 831, 18, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 849, 21, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 867, 24, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 885, 51, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 915, 57, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 945, 63, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 975, 69, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1005, 115, 195, 210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1050, 125, 210, 225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1095, 135, 225, 240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1140, 145, 240, 255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1185, 210, 315, 336, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1248, 225, 336, 357, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1311, 240, 357, 378, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1374, 255, 378, 399, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1437, 336, 462, 490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1521, 357, 490, 518, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1605, 378, 518, 546, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1689, 399, 546, 574, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 1773, 490, 630, 666, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 1881, 518, 666, 702, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 1989, 546, 702, 738, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsk(pbuffer, 2097, 574, 738, 774, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 2205, 18, 21, 810, 813, 822, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2223, 51, 57, 813, 831, 849, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2259, 57, 63, 822, 849, 867, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2295, 105, 115, 831, 885, 915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2355, 115, 125, 849, 915, 945, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2415, 125, 135, 867, 945, 975, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2475, 195, 210, 915, 1005, 1050, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2565, 210, 225, 945, 1050, 1095, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2655, 225, 240, 975, 1095, 1140, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2745, 315, 336, 1050, 1185, 1248, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2871, 336, 357, 1095, 1248, 1311, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2997, 357, 378, 1140, 1311, 1374, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3123, 462, 490, 1248, 1437, 1521, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3291, 490, 518, 1311, 1521, 1605, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3459, 518, 546, 1374, 1605, 1689, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 3627, 630, 666, 1521, 1773, 1881, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 3843, 666, 702, 1605, 1881, 1989, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsk(pbuffer, 4059, 702, 738, 1689, 1989, 2097, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4275, 831, 849, 2205, 2223, 2259, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4335, 885, 915, 2223, 2295, 2355, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4435, 915, 945, 2259, 2355, 2415, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4535, 1005, 1050, 2355, 2475, 2565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4685, 1050, 1095, 2415, 2565, 2655, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4835, 1185, 1248, 2565, 2745, 2871, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5045, 1248, 1311, 2655, 2871, 2997, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 5255, 1437, 1521, 2871, 3123, 3291, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 5535, 1521, 1605, 2997, 3291, 3459, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 5815, 1773, 1881, 3291, 3627, 3843, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsk(pbuffer, 6175, 1881, 1989, 3459, 3843, 4059, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 6535, 2295, 2355, 4275, 4335, 4435, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6685, 2475, 2565, 4435, 4535, 4685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 6910, 2745, 2871, 4685, 4835, 5045, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 7225, 3123, 3291, 5045, 5255, 5535, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsk(pbuffer, 7645, 3627, 3843, 5535, 5815, 6175, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 885, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 1005, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 75, pbuffer, 1185, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 138, pbuffer, 2295, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 198, pbuffer, 2475, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 288, pbuffer, 2745, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {885, 915});

                pbuffer.scale(2.0 * a_exp, {1005, 1050});

                pbuffer.scale(2.0 * a_exp, {1185, 1248});

                pbuffer.scale(2.0 * a_exp, {2295, 2355});

                pbuffer.scale(2.0 * a_exp, {2475, 2565});

                pbuffer.scale(2.0 * a_exp, {2745, 2871});

                pbuffer.scale(2.0 * a_exp, {4335, 4435});

                pbuffer.scale(2.0 * a_exp, {4535, 4685});

                pbuffer.scale(2.0 * a_exp, {4835, 5045});

                pbuffer.scale(2.0 * a_exp, {6535, 6685});

                pbuffer.scale(2.0 * a_exp, {6685, 6910});

                pbuffer.scale(2.0 * a_exp, {6910, 7225});

                t2cfunc::reduce(cbuffer, 1404, pbuffer, 885, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1434, pbuffer, 1005, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1479, pbuffer, 1185, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1542, pbuffer, 2295, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1602, pbuffer, 2475, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1692, pbuffer, 2745, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1818, pbuffer, 4335, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1918, pbuffer, 4535, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2068, pbuffer, 4835, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2278, pbuffer, 6535, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2428, pbuffer, 6685, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2653, pbuffer, 6910, 315, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {885, 915});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1005, 1050});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1185, 1248});

                pbuffer.scale(pfactors, 0, 2.0, {1437, 1521});

                pbuffer.scale(pfactors, 0, 2.0, {1773, 1881});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2295, 2355});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2475, 2565});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2745, 2871});

                pbuffer.scale(pfactors, 0, 2.0, {3123, 3291});

                pbuffer.scale(pfactors, 0, 2.0, {3627, 3843});

                t2cfunc::reduce(cbuffer, 414, pbuffer, 885, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 444, pbuffer, 1005, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 489, pbuffer, 1185, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 552, pbuffer, 1437, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 636, pbuffer, 1773, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 744, pbuffer, 2295, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 804, pbuffer, 2475, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 894, pbuffer, 2745, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1020, pbuffer, 3123, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1188, pbuffer, 3627, 216, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {885, 915});

                pbuffer.scale(2.0 * a_exp, {1005, 1050});

                pbuffer.scale(2.0 * a_exp, {1185, 1248});

                pbuffer.scale(2.0 * a_exp, {1437, 1521});

                pbuffer.scale(2.0 * a_exp, {1773, 1881});

                pbuffer.scale(2.0 * a_exp, {2295, 2355});

                pbuffer.scale(2.0 * a_exp, {2475, 2565});

                pbuffer.scale(2.0 * a_exp, {2745, 2871});

                pbuffer.scale(2.0 * a_exp, {3123, 3291});

                pbuffer.scale(2.0 * a_exp, {3627, 3843});

                pbuffer.scale(pfactors, 0, 2.0, {4335, 4435});

                pbuffer.scale(pfactors, 0, 2.0, {4535, 4685});

                pbuffer.scale(pfactors, 0, 2.0, {4835, 5045});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5255, 5535});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5815, 6175});

                pbuffer.scale(pfactors, 0, 2.0, {6535, 6685});

                pbuffer.scale(pfactors, 0, 2.0, {6685, 6910});

                pbuffer.scale(pfactors, 0, 2.0, {6910, 7225});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {7225, 7645});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {7645, 8185});

                t2cfunc::reduce(cbuffer, 2968, pbuffer, 885, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2998, pbuffer, 1005, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3043, pbuffer, 1185, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3106, pbuffer, 1437, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3190, pbuffer, 1773, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3298, pbuffer, 2295, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3358, pbuffer, 2475, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3448, pbuffer, 2745, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3574, pbuffer, 3123, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3742, pbuffer, 3627, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3958, pbuffer, 4335, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4058, pbuffer, 4535, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4208, pbuffer, 4835, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4418, pbuffer, 5255, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4698, pbuffer, 5815, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5058, pbuffer, 6535, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5208, pbuffer, 6685, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5433, pbuffer, 6910, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5748, pbuffer, 7225, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6168, pbuffer, 7645, 540, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 666, cbuffer, 0, 30, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 1026, cbuffer, 30, 75, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 2133, 666, 1026, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5895, cbuffer, 138, 198, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 6615, cbuffer, 198, 288, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 8829, 5895, 6615, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 14355, cbuffer, 1404, 1434, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 14715, cbuffer, 1434, 1479, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 15822, 14355, 14715, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 19584, cbuffer, 1542, 1602, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 20304, cbuffer, 1602, 1692, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 22518, 19584, 20304, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 29598, cbuffer, 1818, 1918, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 30798, cbuffer, 1918, 2068, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 34488, 29598, 30798, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 45918, cbuffer, 2278, 2428, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 47718, cbuffer, 2428, 2653, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 53253, 45918, 47718, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 0, cbuffer, 414, 444, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 90, cbuffer, 444, 489, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 225, cbuffer, 489, 552, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 414, cbuffer, 552, 636, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 756, cbuffer, 0, 0, 90, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 1161, cbuffer, 30, 90, 225, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 1566, cbuffer, 75, 225, 414, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 2313, 666, 756, 1161, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 2853, 1026, 1161, 1566, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 3663, 2133, 2313, 2853, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 4563, cbuffer, 744, 804, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 4743, cbuffer, 804, 894, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 5013, cbuffer, 894, 1020, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 5391, cbuffer, 1020, 1188, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 6075, cbuffer, 138, 4563, 4743, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 6885, cbuffer, 198, 4743, 5013, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 7695, cbuffer, 288, 5013, 5391, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 9189, 5895, 6075, 6885, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 10269, 6615, 6885, 7695, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 11889, 8829, 9189, 10269, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 13689, cbuffer, 2968, 2998, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 13779, cbuffer, 2998, 3043, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 13914, cbuffer, 3043, 3106, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 14103, cbuffer, 3106, 3190, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 14445, cbuffer, 1404, 13689, 13779, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 14850, cbuffer, 1434, 13779, 13914, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 15255, cbuffer, 1479, 13914, 14103, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 16002, 14355, 14445, 14850, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 16542, 14715, 14850, 15255, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 17352, 15822, 16002, 16542, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 18252, cbuffer, 3298, 3358, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 18432, cbuffer, 3358, 3448, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 18702, cbuffer, 3448, 3574, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 19080, cbuffer, 3574, 3742, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 19764, cbuffer, 1542, 18252, 18432, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 20574, cbuffer, 1602, 18432, 18702, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 21384, cbuffer, 1692, 18702, 19080, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 22878, 19584, 19764, 20574, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 23958, 20304, 20574, 21384, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 25578, 22518, 22878, 23958, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 27378, cbuffer, 3958, 4058, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 27678, cbuffer, 4058, 4208, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 28128, cbuffer, 4208, 4418, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 28758, cbuffer, 4418, 4698, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 29898, cbuffer, 1818, 27378, 27678, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 31248, cbuffer, 1918, 27678, 28128, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 32598, cbuffer, 2068, 28128, 28758, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 35088, 29598, 29898, 31248, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 36888, 30798, 31248, 32598, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 39588, 34488, 35088, 36888, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 42588, cbuffer, 5058, 5208, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 43038, cbuffer, 5208, 5433, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 43713, cbuffer, 5433, 5748, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 44658, cbuffer, 5748, 6168, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 46368, cbuffer, 2278, 42588, 43038, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 48393, cbuffer, 2428, 43038, 43713, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 50418, cbuffer, 2653, 43713, 44658, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 54153, 45918, 46368, 48393, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 56853, 47718, 48393, 50418, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 60903, 53253, 54153, 56853, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 3663, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 147, ckbuffer, 3963, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 294, ckbuffer, 4263, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 1764, ckbuffer, 11889, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 2058, ckbuffer, 12489, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 2352, ckbuffer, 13089, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 30870, ckbuffer, 17352, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 31017, ckbuffer, 17652, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 31164, ckbuffer, 17952, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 31311, ckbuffer, 25578, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 31605, ckbuffer, 26178, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 31899, ckbuffer, 26778, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 32193, ckbuffer, 39588, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 32683, ckbuffer, 40588, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 33173, ckbuffer, 41588, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 33663, ckbuffer, 60903, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 34398, ckbuffer, 62403, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 35133, ckbuffer, 63903, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 9702, 0, 1764, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 10143, 147, 2058, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 10584, 294, 2352, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 441, 30870, 31311, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 2646, 31311, 32193, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 5292, 32193, 33663, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 11025, 0, 441, 2646, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 14994, 1764, 2646, 5292, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 22932, 9702, 11025, 14994, r_ab, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 22932, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 735, skbuffer, 23814, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1470, skbuffer, 24696, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2205, skbuffer, 25578, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2940, skbuffer, 26460, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 3675, skbuffer, 27342, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 4410, skbuffer, 28224, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 5145, skbuffer, 29106, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 5880, skbuffer, 29988, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDPFF_hpp */
