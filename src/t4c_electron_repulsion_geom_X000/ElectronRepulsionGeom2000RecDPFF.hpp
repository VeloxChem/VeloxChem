#ifndef ElectronRepulsionGeom2000RecDPFF_hpp
#define ElectronRepulsionGeom2000RecDPFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDDXX.hpp"
#include "ElectronRepulsionContrRecDFXX.hpp"
#include "ElectronRepulsionContrRecDPXX.hpp"
#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecXXDF.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXFF.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"
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

/// @brief Computes d^(2)/dA^(2)(DP|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_dpff(T& distributor,
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

    CSimdArray<double> cbuffer(5698, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(29876, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(35819, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(4410, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 885, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 75, pbuffer, 1110, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 138, pbuffer, 1425, 84, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {735, 765});

                pbuffer.scale(2.0 * a_exp, {885, 930});

                pbuffer.scale(2.0 * a_exp, {1110, 1173});

                pbuffer.scale(2.0 * a_exp, {1425, 1509});

                pbuffer.scale(2.0 * a_exp, {1995, 2055});

                pbuffer.scale(2.0 * a_exp, {2235, 2325});

                pbuffer.scale(2.0 * a_exp, {2595, 2721});

                pbuffer.scale(2.0 * a_exp, {3099, 3267});

                pbuffer.scale(2.0 * a_exp, {3921, 4021});

                pbuffer.scale(2.0 * a_exp, {4221, 4371});

                pbuffer.scale(2.0 * a_exp, {4671, 4881});

                pbuffer.scale(2.0 * a_exp, {5301, 5581});

                t2cfunc::reduce(cbuffer, 222, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 252, pbuffer, 885, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 297, pbuffer, 1110, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 360, pbuffer, 1425, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 444, pbuffer, 1995, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 504, pbuffer, 2235, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 594, pbuffer, 2595, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 720, pbuffer, 3099, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 888, pbuffer, 3921, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 988, pbuffer, 4221, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1138, pbuffer, 4671, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1348, pbuffer, 5301, 280, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {735, 765});

                pbuffer.scale(2.0 * a_exp, {885, 930});

                pbuffer.scale(2.0 * a_exp, {1110, 1173});

                pbuffer.scale(2.0 * a_exp, {1425, 1509});

                pbuffer.scale(2.0 * a_exp, {1995, 2055});

                pbuffer.scale(2.0 * a_exp, {2235, 2325});

                pbuffer.scale(2.0 * a_exp, {2595, 2721});

                pbuffer.scale(2.0 * a_exp, {3099, 3267});

                pbuffer.scale(2.0 * a_exp, {3921, 4021});

                pbuffer.scale(2.0 * a_exp, {4221, 4371});

                pbuffer.scale(2.0 * a_exp, {4671, 4881});

                pbuffer.scale(2.0 * a_exp, {5301, 5581});

                pbuffer.scale(4.0 * a_exp * a_exp, {6231, 6381});

                pbuffer.scale(4.0 * a_exp * a_exp, {6531, 6756});

                pbuffer.scale(4.0 * a_exp * a_exp, {6981, 7296});

                pbuffer.scale(4.0 * a_exp * a_exp, {7611, 8031});

                pbuffer.scale(4.0 * a_exp * a_exp, {8451, 8661});

                pbuffer.scale(4.0 * a_exp * a_exp, {8661, 8976});

                pbuffer.scale(4.0 * a_exp * a_exp, {8976, 9417});

                pbuffer.scale(4.0 * a_exp * a_exp, {9417, 10005});

                t2cfunc::reduce(cbuffer, 1628, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1658, pbuffer, 885, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1703, pbuffer, 1110, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1766, pbuffer, 1425, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1850, pbuffer, 1995, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1910, pbuffer, 2235, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2000, pbuffer, 2595, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2126, pbuffer, 3099, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2294, pbuffer, 3921, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2394, pbuffer, 4221, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2544, pbuffer, 4671, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2754, pbuffer, 5301, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3034, pbuffer, 6231, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3184, pbuffer, 6531, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3409, pbuffer, 6981, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3724, pbuffer, 7611, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4144, pbuffer, 8451, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4354, pbuffer, 8661, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4669, pbuffer, 8976, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5110, pbuffer, 9417, 588, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 30, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 90, cbuffer, 30, 75, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 225, cbuffer, 75, 138, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 414, 0, 90, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 594, 90, 225, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 864, 414, 594, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1164, cbuffer, 222, 252, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 1254, cbuffer, 252, 297, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 1389, cbuffer, 297, 360, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 1578, 1164, 1254, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 1758, 1254, 1389, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 2028, 1578, 1758, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2328, cbuffer, 444, 504, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 2508, cbuffer, 504, 594, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 2778, cbuffer, 594, 720, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 3156, 2328, 2508, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 3516, 2508, 2778, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 4056, 3156, 3516, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4656, cbuffer, 888, 988, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 4956, cbuffer, 988, 1138, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 5406, cbuffer, 1138, 1348, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 6036, 4656, 4956, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 6636, 4956, 5406, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 7536, 6036, 6636, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 8536, cbuffer, 1628, 1658, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 8626, cbuffer, 1658, 1703, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 8761, cbuffer, 1703, 1766, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 8950, 8536, 8626, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 9130, 8626, 8761, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 9400, 8950, 9130, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9700, cbuffer, 1850, 1910, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 9880, cbuffer, 1910, 2000, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 10150, cbuffer, 2000, 2126, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 10528, 9700, 9880, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 10888, 9880, 10150, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 11428, 10528, 10888, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 12028, cbuffer, 2294, 2394, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 12328, cbuffer, 2394, 2544, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 12778, cbuffer, 2544, 2754, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 13408, 12028, 12328, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 14008, 12328, 12778, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 14908, 13408, 14008, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 15908, cbuffer, 3034, 3184, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 16358, cbuffer, 3184, 3409, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 17033, cbuffer, 3409, 3724, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 17978, 15908, 16358, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 18878, 16358, 17033, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 20228, 17978, 18878, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 21728, cbuffer, 4144, 4354, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 22358, cbuffer, 4354, 4669, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 23303, cbuffer, 4669, 5110, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 24626, 21728, 22358, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 25886, 22358, 23303, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 27776, 24626, 25886, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 864, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 20286, ckbuffer, 2028, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 20433, ckbuffer, 4056, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 20727, ckbuffer, 7536, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 22540, ckbuffer, 9400, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 22687, ckbuffer, 11428, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 22981, ckbuffer, 14908, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 23471, ckbuffer, 20228, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 24206, ckbuffer, 27776, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 21217, 20286, 20433, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 21658, 20433, 20727, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 25235, 22540, 22687, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 25676, 22687, 22981, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 26558, 22981, 23471, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 28028, 23471, 24206, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 30233, 25235, 25676, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 31115, 25676, 26558, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 32879, 26558, 28028, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 5733, 0, 21217, 21658, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 147, 20286, 30233, 1, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 1029, 20433, 31115, 2, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 2793, 20727, 32879, 3, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_ppxx(skbuffer, 7056, 21217, 147, 1029, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_pdxx(skbuffer, 9702, 21658, 1029, 2793, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dpxx(skbuffer, 14994, 5733, 7056, 9702, r_ab, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 14994, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 735, skbuffer, 15876, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1470, skbuffer, 16758, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2205, skbuffer, 17640, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2940, skbuffer, 18522, 3, 3);

            t4cfunc::bra_transform<2, 1>(sbuffer, 3675, skbuffer, 19404, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecDPFF_hpp */
