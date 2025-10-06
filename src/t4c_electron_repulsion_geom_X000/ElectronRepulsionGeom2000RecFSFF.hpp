#ifndef ElectronRepulsionGeom2000RecFSFF_hpp
#define ElectronRepulsionGeom2000RecFSFF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDDXX.hpp"
#include "ElectronRepulsionContrRecDFXX.hpp"
#include "ElectronRepulsionContrRecDPXX.hpp"
#include "ElectronRepulsionContrRecDSXX.hpp"
#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPSXX.hpp"
#include "ElectronRepulsionContrRecXXDF.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXFF.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
#include "ElectronRepulsionGeom1000ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecFSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPSXX.hpp"
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

/// @brief Computes d^(2)/dA^(2)(FS|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_fsff(T& distributor,
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

    CSimdArray<double> cbuffer(5920, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(31040, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(43904, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2058, 1);

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

                pbuffer.scale(2.0 * a_exp, {105, 115});

                pbuffer.scale(2.0 * a_exp, {195, 210});

                pbuffer.scale(2.0 * a_exp, {315, 336});

                pbuffer.scale(2.0 * a_exp, {462, 490});

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

                t2cfunc::reduce(cbuffer, 296, pbuffer, 105, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 306, pbuffer, 195, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 321, pbuffer, 315, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 342, pbuffer, 462, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 370, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 400, pbuffer, 885, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 445, pbuffer, 1110, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 508, pbuffer, 1425, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 592, pbuffer, 1995, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 652, pbuffer, 2235, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 742, pbuffer, 2595, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 868, pbuffer, 3099, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1036, pbuffer, 3921, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1136, pbuffer, 4221, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1286, pbuffer, 4671, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1496, pbuffer, 5301, 280, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {105, 115});

                pbuffer.scale(2.0 * a_exp, {195, 210});

                pbuffer.scale(2.0 * a_exp, {315, 336});

                pbuffer.scale(2.0 * a_exp, {462, 490});

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

                t2cfunc::reduce(cbuffer, 1776, pbuffer, 105, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1786, pbuffer, 195, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1801, pbuffer, 315, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1822, pbuffer, 462, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1850, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1880, pbuffer, 885, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1925, pbuffer, 1110, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1988, pbuffer, 1425, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2072, pbuffer, 1995, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2132, pbuffer, 2235, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2222, pbuffer, 2595, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2348, pbuffer, 3099, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2516, pbuffer, 3921, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2616, pbuffer, 4221, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2766, pbuffer, 4671, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2976, pbuffer, 5301, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3256, pbuffer, 6231, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3406, pbuffer, 6531, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3631, pbuffer, 6981, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3946, pbuffer, 7611, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4366, pbuffer, 8451, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4576, pbuffer, 8661, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4891, pbuffer, 8976, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5332, pbuffer, 9417, 588, ket_width, ket_npgtos);

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

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1552, cbuffer, 296, 306, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 1582, cbuffer, 306, 321, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 1627, cbuffer, 321, 342, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 1690, 1552, 1582, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 1750, 1582, 1627, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 1840, 1690, 1750, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1940, cbuffer, 370, 400, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 2030, cbuffer, 400, 445, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 2165, cbuffer, 445, 508, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 2354, 1940, 2030, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 2534, 2030, 2165, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 2804, 2354, 2534, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 3104, cbuffer, 592, 652, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 3284, cbuffer, 652, 742, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 3554, cbuffer, 742, 868, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 3932, 3104, 3284, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 4292, 3284, 3554, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 4832, 3932, 4292, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5432, cbuffer, 1036, 1136, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 5732, cbuffer, 1136, 1286, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 6182, cbuffer, 1286, 1496, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 6812, 5432, 5732, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 7412, 5732, 6182, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 8312, 6812, 7412, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9312, cbuffer, 1776, 1786, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 9342, cbuffer, 1786, 1801, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 9387, cbuffer, 1801, 1822, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 9450, 9312, 9342, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 9510, 9342, 9387, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 9600, 9450, 9510, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9700, cbuffer, 1850, 1880, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 9790, cbuffer, 1880, 1925, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 9925, cbuffer, 1925, 1988, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 10114, 9700, 9790, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 10294, 9790, 9925, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 10564, 10114, 10294, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 10864, cbuffer, 2072, 2132, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 11044, cbuffer, 2132, 2222, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 11314, cbuffer, 2222, 2348, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 11692, 10864, 11044, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 12052, 11044, 11314, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 12592, 11692, 12052, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 13192, cbuffer, 2516, 2616, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 13492, cbuffer, 2616, 2766, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 13942, cbuffer, 2766, 2976, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 14572, 13192, 13492, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 15172, 13492, 13942, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 16072, 14572, 15172, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 17072, cbuffer, 3256, 3406, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 17522, cbuffer, 3406, 3631, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 18197, cbuffer, 3631, 3946, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 19142, 17072, 17522, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 20042, 17522, 18197, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 21392, 19142, 20042, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 22892, cbuffer, 4366, 4576, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 23522, cbuffer, 4576, 4891, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxph(ckbuffer, 24467, cbuffer, 4891, 5332, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 25790, 22892, 23522, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdg(ckbuffer, 27050, 23522, 24467, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxff(ckbuffer, 28940, 25790, 27050, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 288, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 343, ckbuffer, 1252, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 27685, ckbuffer, 1840, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 27734, ckbuffer, 2804, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 27881, ckbuffer, 4832, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 28175, ckbuffer, 8312, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 30135, ckbuffer, 9600, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 30184, ckbuffer, 10564, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 30331, ckbuffer, 12592, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 30625, ckbuffer, 16072, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 31115, ckbuffer, 21392, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 31850, ckbuffer, 28940, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 6076, 0, 343, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 28665, 27685, 27734, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 28812, 27734, 27881, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 29253, 27881, 28175, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 32879, 30135, 30184, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 33026, 30184, 30331, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 33467, 30331, 30625, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 34349, 30625, 31115, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 35819, 31115, 31850, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 38024, 32879, 33026, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 38318, 33026, 33467, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 39200, 33467, 34349, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 40964, 34349, 35819, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_psxx(skbuffer, 6223, 0, 28665, 28812, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 7546, 343, 28812, 29253, r_ab, 3, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dsxx(skbuffer, 16807, 6076, 6223, 7546, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 49, 27685, 38024, 0, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 490, 27734, 38318, 1, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 1372, 27881, 39200, 2, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 3136, 28175, 40964, 3, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_psxx(skbuffer, 6664, 28665, 49, 490, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_ppxx(skbuffer, 8869, 28812, 490, 1372, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_pdxx(skbuffer, 11515, 29253, 1372, 3136, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dsxx(skbuffer, 17689, 6223, 6664, 8869, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dpxx(skbuffer, 19453, 7546, 8869, 11515, r_ab, 3, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_fsxx(skbuffer, 24745, 16807, 17689, 19453, r_ab, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 24745, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 343, skbuffer, 25235, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 686, skbuffer, 25725, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1029, skbuffer, 26215, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1372, skbuffer, 26705, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1715, skbuffer, 27195, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecFSFF_hpp */
