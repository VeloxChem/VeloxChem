#ifndef ElectronRepulsionGeom1010RecFSFF_hpp
#define ElectronRepulsionGeom1010RecFSFF_hpp

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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FS|1/|r-r'||FF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fsff(T& distributor,
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

    CSimdArray<double> cbuffer(7020, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(68445, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(46305, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3087, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 105, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 195, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 25, pbuffer, 315, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 46, pbuffer, 885, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 76, pbuffer, 1005, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 121, pbuffer, 1185, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 184, pbuffer, 2295, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 244, pbuffer, 2475, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 334, pbuffer, 2745, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {105, 115});

                pbuffer.scale(2.0 * a_exp, {195, 210});

                pbuffer.scale(2.0 * a_exp, {315, 336});

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

                t2cfunc::reduce(cbuffer, 1560, pbuffer, 105, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1570, pbuffer, 195, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1585, pbuffer, 315, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1606, pbuffer, 885, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1636, pbuffer, 1005, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1681, pbuffer, 1185, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1744, pbuffer, 2295, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1804, pbuffer, 2475, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1894, pbuffer, 2745, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2020, pbuffer, 4335, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2120, pbuffer, 4535, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2270, pbuffer, 4835, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2480, pbuffer, 6535, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2630, pbuffer, 6685, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2855, pbuffer, 6910, 315, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {105, 115});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {195, 210});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {315, 336});

                pbuffer.scale(pfactors, 0, 2.0, {462, 490});

                pbuffer.scale(pfactors, 0, 2.0, {630, 666});

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

                t2cfunc::reduce(cbuffer, 460, pbuffer, 105, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 470, pbuffer, 195, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 485, pbuffer, 315, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 506, pbuffer, 462, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 534, pbuffer, 630, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 570, pbuffer, 885, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 600, pbuffer, 1005, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 645, pbuffer, 1185, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 708, pbuffer, 1437, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 792, pbuffer, 1773, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 900, pbuffer, 2295, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 960, pbuffer, 2475, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1050, pbuffer, 2745, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1176, pbuffer, 3123, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1344, pbuffer, 3627, 216, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {105, 115});

                pbuffer.scale(2.0 * a_exp, {195, 210});

                pbuffer.scale(2.0 * a_exp, {315, 336});

                pbuffer.scale(2.0 * a_exp, {462, 490});

                pbuffer.scale(2.0 * a_exp, {630, 666});

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

                t2cfunc::reduce(cbuffer, 3170, pbuffer, 105, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3180, pbuffer, 195, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3195, pbuffer, 315, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3216, pbuffer, 462, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3244, pbuffer, 630, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3280, pbuffer, 885, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3310, pbuffer, 1005, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3355, pbuffer, 1185, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3418, pbuffer, 1437, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3502, pbuffer, 1773, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3610, pbuffer, 2295, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3670, pbuffer, 2475, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3760, pbuffer, 2745, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3886, pbuffer, 3123, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4054, pbuffer, 3627, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4270, pbuffer, 4335, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4370, pbuffer, 4535, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4520, pbuffer, 4835, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4730, pbuffer, 5255, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5010, pbuffer, 5815, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5370, pbuffer, 6535, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5520, pbuffer, 6685, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5745, pbuffer, 6910, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6060, pbuffer, 7225, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6480, pbuffer, 7645, 540, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 222, cbuffer, 0, 10, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 342, cbuffer, 10, 25, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 711, 222, 342, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2187, cbuffer, 46, 76, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 2547, cbuffer, 76, 121, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 3654, 2187, 2547, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 7416, cbuffer, 184, 244, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 8136, cbuffer, 244, 334, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 10350, 7416, 8136, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 15432, cbuffer, 1560, 1570, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 15552, cbuffer, 1570, 1585, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 15921, 15432, 15552, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 17397, cbuffer, 1606, 1636, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 17757, cbuffer, 1636, 1681, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 18864, 17397, 17757, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 22626, cbuffer, 1744, 1804, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 23346, cbuffer, 1804, 1894, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 25560, 22626, 23346, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 32640, cbuffer, 2020, 2120, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 33840, cbuffer, 2120, 2270, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 37530, 32640, 33840, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 48960, cbuffer, 2480, 2630, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 50760, cbuffer, 2630, 2855, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 56295, 48960, 50760, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 0, cbuffer, 460, 470, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 30, cbuffer, 470, 485, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 75, cbuffer, 485, 506, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 138, cbuffer, 506, 534, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 252, cbuffer, 0, 0, 30, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 387, cbuffer, 10, 30, 75, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 522, cbuffer, 25, 75, 138, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 771, 222, 252, 387, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 951, 342, 387, 522, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 1221, 711, 771, 951, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1521, cbuffer, 570, 600, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 1611, cbuffer, 600, 645, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 1746, cbuffer, 645, 708, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 1935, cbuffer, 708, 792, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2277, cbuffer, 46, 1521, 1611, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 2682, cbuffer, 76, 1611, 1746, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 3087, cbuffer, 121, 1746, 1935, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 3834, 2187, 2277, 2682, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 4374, 2547, 2682, 3087, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 5184, 3654, 3834, 4374, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6084, cbuffer, 900, 960, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 6264, cbuffer, 960, 1050, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 6534, cbuffer, 1050, 1176, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 6912, cbuffer, 1176, 1344, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 7596, cbuffer, 184, 6084, 6264, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 8406, cbuffer, 244, 6264, 6534, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 9216, cbuffer, 334, 6534, 6912, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 10710, 7416, 7596, 8406, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 11790, 8136, 8406, 9216, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 13410, 10350, 10710, 11790, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 15210, cbuffer, 3170, 3180, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 15240, cbuffer, 3180, 3195, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 15285, cbuffer, 3195, 3216, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 15348, cbuffer, 3216, 3244, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 15462, cbuffer, 1560, 15210, 15240, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 15597, cbuffer, 1570, 15240, 15285, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 15732, cbuffer, 1585, 15285, 15348, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 15981, 15432, 15462, 15597, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 16161, 15552, 15597, 15732, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 16431, 15921, 15981, 16161, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 16731, cbuffer, 3280, 3310, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 16821, cbuffer, 3310, 3355, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 16956, cbuffer, 3355, 3418, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 17145, cbuffer, 3418, 3502, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 17487, cbuffer, 1606, 16731, 16821, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 17892, cbuffer, 1636, 16821, 16956, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 18297, cbuffer, 1681, 16956, 17145, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 19044, 17397, 17487, 17892, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 19584, 17757, 17892, 18297, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 20394, 18864, 19044, 19584, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 21294, cbuffer, 3610, 3670, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 21474, cbuffer, 3670, 3760, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 21744, cbuffer, 3760, 3886, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 22122, cbuffer, 3886, 4054, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 22806, cbuffer, 1744, 21294, 21474, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 23616, cbuffer, 1804, 21474, 21744, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 24426, cbuffer, 1894, 21744, 22122, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 25920, 22626, 22806, 23616, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 27000, 23346, 23616, 24426, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 28620, 25560, 25920, 27000, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 30420, cbuffer, 4270, 4370, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 30720, cbuffer, 4370, 4520, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 31170, cbuffer, 4520, 4730, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 31800, cbuffer, 4730, 5010, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 32940, cbuffer, 2020, 30420, 30720, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 34290, cbuffer, 2120, 30720, 31170, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 35640, cbuffer, 2270, 31170, 31800, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 38130, 32640, 32940, 34290, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 39930, 33840, 34290, 35640, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 42630, 37530, 38130, 39930, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 45630, cbuffer, 5370, 5520, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 46080, cbuffer, 5520, 5745, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 46755, cbuffer, 5745, 6060, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsi(ckbuffer, 47700, cbuffer, 6060, 6480, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 49410, cbuffer, 2480, 45630, 46080, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 51435, cbuffer, 2630, 46080, 46755, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxph(ckbuffer, 53460, cbuffer, 2855, 46755, 47700, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 57195, 48960, 49410, 51435, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdg(ckbuffer, 59895, 50760, 51435, 53460, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxff(ckbuffer, 63945, 56295, 57195, 59895, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 0, ckbuffer, 1221, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 49, ckbuffer, 1321, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 98, ckbuffer, 1421, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 588, ckbuffer, 5184, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 735, ckbuffer, 5484, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 882, ckbuffer, 5784, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 2352, ckbuffer, 13410, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 2646, ckbuffer, 14010, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 2940, ckbuffer, 14610, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 41160, ckbuffer, 16431, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 41209, ckbuffer, 16531, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 41258, ckbuffer, 16631, 0, 0);

            t4cfunc::ket_transform<3, 3>(skbuffer, 41307, ckbuffer, 20394, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 41454, ckbuffer, 20694, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 41601, ckbuffer, 20994, 0, 1);

            t4cfunc::ket_transform<3, 3>(skbuffer, 41748, ckbuffer, 28620, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 42042, ckbuffer, 29220, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 42336, ckbuffer, 29820, 0, 2);

            t4cfunc::ket_transform<3, 3>(skbuffer, 42630, ckbuffer, 42630, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 43120, ckbuffer, 43630, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 43610, ckbuffer, 44630, 0, 3);

            t4cfunc::ket_transform<3, 3>(skbuffer, 44100, ckbuffer, 63945, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 44835, ckbuffer, 65445, 0, 4);

            t4cfunc::ket_transform<3, 3>(skbuffer, 45570, ckbuffer, 66945, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 10290, 0, 588, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 10437, 49, 735, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 10584, 98, 882, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 12054, 588, 2352, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 12495, 735, 2646, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 12936, 882, 2940, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 25284, 10290, 12054, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 25578, 10437, 12495, r_ab, 3, 3);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 25872, 10584, 12936, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 147, 41160, 41307, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 1029, 41307, 41748, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 3234, 41748, 42630, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 5880, 42630, 44100, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 10731, 0, 147, 1029, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 13377, 588, 1029, 3234, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 17346, 2352, 3234, 5880, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 26166, 10290, 10731, 13377, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 28812, 12054, 13377, 17346, r_ab, 3, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fsxx(skbuffer, 36750, 25284, 26166, 28812, r_ab, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 36750, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 343, skbuffer, 37240, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 686, skbuffer, 37730, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1029, skbuffer, 38220, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1372, skbuffer, 38710, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1715, skbuffer, 39200, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 2058, skbuffer, 39690, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 2401, skbuffer, 40180, 3, 3);

            t4cfunc::bra_transform<3, 0>(sbuffer, 2744, skbuffer, 40670, 3, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 3, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFSFF_hpp */
