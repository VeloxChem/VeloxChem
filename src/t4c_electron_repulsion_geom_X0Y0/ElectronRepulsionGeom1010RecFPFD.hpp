#ifndef ElectronRepulsionGeom1010RecFPFD_hpp
#define ElectronRepulsionGeom1010RecFPFD_hpp

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
#include "ElectronRepulsionGeom1010ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
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
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSI.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FP|1/|r-r'||FD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fpfd(T& distributor,
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

    CSimdArray<double> pbuffer(10456, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(8214, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(73038, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(67620, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 630, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 633, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 636, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 639, 2, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 648, 3, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 657, 4, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 666, 5, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 675, 15, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 693, 18, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 711, 21, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 729, 24, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 747, 27, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 765, 51, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 795, 57, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 825, 63, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 855, 69, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 885, 75, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 915, 115, 195, 210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 960, 125, 210, 225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1005, 135, 225, 240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1050, 145, 240, 255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1095, 155, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1140, 210, 315, 336, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1203, 225, 336, 357, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1266, 240, 357, 378, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1329, 255, 378, 399, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1392, 270, 399, 420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1455, 336, 462, 490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1539, 357, 490, 518, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1623, 378, 518, 546, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1707, 399, 546, 574, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1791, 420, 574, 602, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1875, 2, 3, 630, 633, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1881, 3, 4, 633, 636, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1887, 15, 18, 630, 639, 648, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1905, 18, 21, 633, 648, 657, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1923, 21, 24, 636, 657, 666, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1941, 45, 51, 639, 675, 693, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1977, 51, 57, 648, 693, 711, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2013, 57, 63, 657, 711, 729, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 2049, 63, 69, 666, 729, 747, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2085, 105, 115, 693, 765, 795, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2145, 115, 125, 711, 795, 825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2205, 125, 135, 729, 825, 855, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2265, 135, 145, 747, 855, 885, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2325, 195, 210, 795, 915, 960, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2415, 210, 225, 825, 960, 1005, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2505, 225, 240, 855, 1005, 1050, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2595, 240, 255, 885, 1050, 1095, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2685, 315, 336, 960, 1140, 1203, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2811, 336, 357, 1005, 1203, 1266, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2937, 357, 378, 1050, 1266, 1329, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3063, 378, 399, 1095, 1329, 1392, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3189, 462, 490, 1203, 1455, 1539, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3357, 490, 518, 1266, 1539, 1623, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3525, 518, 546, 1329, 1623, 1707, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 3693, 546, 574, 1392, 1707, 1791, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3861, 630, 633, 1875, 1881, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3871, 639, 648, 1875, 1887, 1905, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3901, 648, 657, 1881, 1905, 1923, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3931, 675, 693, 1887, 1941, 1977, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3991, 693, 711, 1905, 1977, 2013, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 4051, 711, 729, 1923, 2013, 2049, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4111, 765, 795, 1977, 2085, 2145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4211, 795, 825, 2013, 2145, 2205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 4311, 825, 855, 2049, 2205, 2265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4411, 915, 960, 2145, 2325, 2415, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4561, 960, 1005, 2205, 2415, 2505, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4711, 1005, 1050, 2265, 2505, 2595, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4861, 1140, 1203, 2415, 2685, 2811, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5071, 1203, 1266, 2505, 2811, 2937, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5281, 1266, 1329, 2595, 2937, 3063, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 5491, 1455, 1539, 2811, 3189, 3357, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 5771, 1539, 1623, 2937, 3357, 3525, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 6051, 1623, 1707, 3063, 3525, 3693, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 6331, 1887, 1905, 3861, 3871, 3901, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 6376, 1941, 1977, 3871, 3931, 3991, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 6466, 1977, 2013, 3901, 3991, 4051, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 6556, 2085, 2145, 3991, 4111, 4211, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 6706, 2145, 2205, 4051, 4211, 4311, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6856, 2325, 2415, 4211, 4411, 4561, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 7081, 2415, 2505, 4311, 4561, 4711, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 7306, 2685, 2811, 4561, 4861, 5071, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 7621, 2811, 2937, 4711, 5071, 5281, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 7936, 3189, 3357, 5071, 5491, 5771, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsi(pbuffer, 8356, 3357, 3525, 5281, 5771, 6051, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 8776, 3931, 3991, 6331, 6376, 6466, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 8902, 4111, 4211, 6466, 6556, 6706, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 9112, 4411, 4561, 6706, 6856, 7081, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 9427, 4861, 5071, 7081, 7306, 7621, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsi(pbuffer, 9868, 5491, 5771, 7621, 7936, 8356, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 675, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 765, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 48, pbuffer, 915, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 93, pbuffer, 1941, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 129, pbuffer, 2085, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 189, pbuffer, 2325, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 279, pbuffer, 3931, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 339, pbuffer, 4111, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 439, pbuffer, 4411, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {675, 693});

                pbuffer.scale(2.0 * a_exp, {765, 795});

                pbuffer.scale(2.0 * a_exp, {915, 960});

                pbuffer.scale(2.0 * a_exp, {1941, 1977});

                pbuffer.scale(2.0 * a_exp, {2085, 2145});

                pbuffer.scale(2.0 * a_exp, {2325, 2415});

                pbuffer.scale(2.0 * a_exp, {3931, 3991});

                pbuffer.scale(2.0 * a_exp, {4111, 4211});

                pbuffer.scale(2.0 * a_exp, {4411, 4561});

                pbuffer.scale(2.0 * a_exp, {6376, 6466});

                pbuffer.scale(2.0 * a_exp, {6556, 6706});

                pbuffer.scale(2.0 * a_exp, {6856, 7081});

                pbuffer.scale(2.0 * a_exp, {8776, 8902});

                pbuffer.scale(2.0 * a_exp, {8902, 9112});

                pbuffer.scale(2.0 * a_exp, {9112, 9427});

                t2cfunc::reduce(cbuffer, 2109, pbuffer, 675, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2127, pbuffer, 765, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2157, pbuffer, 915, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2202, pbuffer, 1941, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2238, pbuffer, 2085, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2298, pbuffer, 2325, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2388, pbuffer, 3931, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2448, pbuffer, 4111, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2548, pbuffer, 4411, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2698, pbuffer, 6376, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2788, pbuffer, 6556, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2938, pbuffer, 6856, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3163, pbuffer, 8776, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3289, pbuffer, 8902, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3499, pbuffer, 9112, 315, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {675, 693});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {765, 795});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {915, 960});

                pbuffer.scale(pfactors, 0, 2.0, {1140, 1203});

                pbuffer.scale(pfactors, 0, 2.0, {1455, 1539});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1941, 1977});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2085, 2145});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2325, 2415});

                pbuffer.scale(pfactors, 0, 2.0, {2685, 2811});

                pbuffer.scale(pfactors, 0, 2.0, {3189, 3357});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3931, 3991});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {4111, 4211});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {4411, 4561});

                pbuffer.scale(pfactors, 0, 2.0, {4861, 5071});

                pbuffer.scale(pfactors, 0, 2.0, {5491, 5771});

                t2cfunc::reduce(cbuffer, 589, pbuffer, 675, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 607, pbuffer, 765, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 637, pbuffer, 915, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 682, pbuffer, 1140, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 745, pbuffer, 1455, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 829, pbuffer, 1941, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 865, pbuffer, 2085, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 925, pbuffer, 2325, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1015, pbuffer, 2685, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1141, pbuffer, 3189, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1309, pbuffer, 3931, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1369, pbuffer, 4111, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1469, pbuffer, 4411, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1619, pbuffer, 4861, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1829, pbuffer, 5491, 280, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {675, 693});

                pbuffer.scale(2.0 * a_exp, {765, 795});

                pbuffer.scale(2.0 * a_exp, {915, 960});

                pbuffer.scale(2.0 * a_exp, {1140, 1203});

                pbuffer.scale(2.0 * a_exp, {1455, 1539});

                pbuffer.scale(2.0 * a_exp, {1941, 1977});

                pbuffer.scale(2.0 * a_exp, {2085, 2145});

                pbuffer.scale(2.0 * a_exp, {2325, 2415});

                pbuffer.scale(2.0 * a_exp, {2685, 2811});

                pbuffer.scale(2.0 * a_exp, {3189, 3357});

                pbuffer.scale(2.0 * a_exp, {3931, 3991});

                pbuffer.scale(2.0 * a_exp, {4111, 4211});

                pbuffer.scale(2.0 * a_exp, {4411, 4561});

                pbuffer.scale(2.0 * a_exp, {4861, 5071});

                pbuffer.scale(2.0 * a_exp, {5491, 5771});

                pbuffer.scale(pfactors, 0, 2.0, {6376, 6466});

                pbuffer.scale(pfactors, 0, 2.0, {6556, 6706});

                pbuffer.scale(pfactors, 0, 2.0, {6856, 7081});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {7306, 7621});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {7936, 8356});

                pbuffer.scale(pfactors, 0, 2.0, {8776, 8902});

                pbuffer.scale(pfactors, 0, 2.0, {8902, 9112});

                pbuffer.scale(pfactors, 0, 2.0, {9112, 9427});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {9427, 9868});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {9868, 10456});

                t2cfunc::reduce(cbuffer, 3814, pbuffer, 675, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3832, pbuffer, 765, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3862, pbuffer, 915, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3907, pbuffer, 1140, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3970, pbuffer, 1455, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4054, pbuffer, 1941, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4090, pbuffer, 2085, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4150, pbuffer, 2325, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4240, pbuffer, 2685, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4366, pbuffer, 3189, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4534, pbuffer, 3931, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4594, pbuffer, 4111, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4694, pbuffer, 4411, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4844, pbuffer, 4861, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5054, pbuffer, 5491, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5334, pbuffer, 6376, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5424, pbuffer, 6556, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5574, pbuffer, 6856, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5799, pbuffer, 7306, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6114, pbuffer, 7936, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6534, pbuffer, 8776, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6660, pbuffer, 8902, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6870, pbuffer, 9112, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7185, pbuffer, 9427, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7626, pbuffer, 9868, 588, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 468, cbuffer, 0, 18, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 684, cbuffer, 18, 48, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 1449, 468, 684, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3897, cbuffer, 93, 129, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4329, cbuffer, 129, 189, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 5859, 3897, 4329, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 10443, cbuffer, 279, 339, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 11163, cbuffer, 339, 439, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 13713, 10443, 11163, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 19221, cbuffer, 2109, 2127, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 19437, cbuffer, 2127, 2157, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 20202, 19221, 19437, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 22650, cbuffer, 2202, 2238, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 23082, cbuffer, 2238, 2298, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 24612, 22650, 23082, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 29196, cbuffer, 2388, 2448, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 29916, cbuffer, 2448, 2548, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 32466, 29196, 29916, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 39846, cbuffer, 2698, 2788, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 40926, cbuffer, 2788, 2938, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 44751, 39846, 40926, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 55587, cbuffer, 3163, 3289, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 57099, cbuffer, 3289, 3499, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 62454, 55587, 57099, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 0, cbuffer, 589, 607, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 54, cbuffer, 607, 637, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 144, cbuffer, 637, 682, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 279, cbuffer, 682, 745, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 522, cbuffer, 0, 0, 54, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 774, cbuffer, 18, 54, 144, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 1044, cbuffer, 48, 144, 279, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 1557, 468, 522, 774, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 1881, 684, 774, 1044, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 2421, 1449, 1557, 1881, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2961, cbuffer, 829, 865, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3069, cbuffer, 865, 925, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 3249, cbuffer, 925, 1015, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 3519, cbuffer, 1015, 1141, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4005, cbuffer, 93, 2961, 3069, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 4509, cbuffer, 129, 3069, 3249, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 5049, cbuffer, 189, 3249, 3519, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 6075, 3897, 4005, 4509, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 6723, 4329, 4509, 5049, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 7803, 5859, 6075, 6723, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 8883, cbuffer, 1309, 1369, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 9063, cbuffer, 1369, 1469, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 9363, cbuffer, 1469, 1619, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 9813, cbuffer, 1619, 1829, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 10623, cbuffer, 279, 8883, 9063, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 11463, cbuffer, 339, 9063, 9363, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 12363, cbuffer, 439, 9363, 9813, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 14073, 10443, 10623, 11463, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 15153, 11163, 11463, 12363, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 16953, 13713, 14073, 15153, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 18753, cbuffer, 3814, 3832, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 18807, cbuffer, 3832, 3862, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 18897, cbuffer, 3862, 3907, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 19032, cbuffer, 3907, 3970, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 19275, cbuffer, 2109, 18753, 18807, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 19527, cbuffer, 2127, 18807, 18897, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 19797, cbuffer, 2157, 18897, 19032, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 20310, 19221, 19275, 19527, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 20634, 19437, 19527, 19797, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 21174, 20202, 20310, 20634, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 21714, cbuffer, 4054, 4090, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 21822, cbuffer, 4090, 4150, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 22002, cbuffer, 4150, 4240, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 22272, cbuffer, 4240, 4366, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 22758, cbuffer, 2202, 21714, 21822, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 23262, cbuffer, 2238, 21822, 22002, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 23802, cbuffer, 2298, 22002, 22272, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 24828, 22650, 22758, 23262, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 25476, 23082, 23262, 23802, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 26556, 24612, 24828, 25476, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 27636, cbuffer, 4534, 4594, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 27816, cbuffer, 4594, 4694, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 28116, cbuffer, 4694, 4844, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 28566, cbuffer, 4844, 5054, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 29376, cbuffer, 2388, 27636, 27816, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 30216, cbuffer, 2448, 27816, 28116, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 31116, cbuffer, 2548, 28116, 28566, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 32826, 29196, 29376, 30216, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 33906, 29916, 30216, 31116, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 35706, 32466, 32826, 33906, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 37506, cbuffer, 5334, 5424, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 37776, cbuffer, 5424, 5574, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 38226, cbuffer, 5574, 5799, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 38901, cbuffer, 5799, 6114, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 40116, cbuffer, 2698, 37506, 37776, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 41376, cbuffer, 2788, 37776, 38226, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 42726, cbuffer, 2938, 38226, 38901, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 45291, 39846, 40116, 41376, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 46911, 40926, 41376, 42726, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 49611, 44751, 45291, 46911, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 52311, cbuffer, 6534, 6660, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 52689, cbuffer, 6660, 6870, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 53319, cbuffer, 6870, 7185, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 54264, cbuffer, 7185, 7626, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 55965, cbuffer, 3163, 52311, 52689, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 57729, cbuffer, 3289, 52689, 53319, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 59619, cbuffer, 3499, 53319, 54264, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 63210, 55587, 55965, 57729, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 65478, 57099, 57729, 59619, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 69258, 62454, 63210, 65478, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 0, ckbuffer, 2421, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 105, ckbuffer, 2601, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 210, ckbuffer, 2781, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 1260, ckbuffer, 7803, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 1470, ckbuffer, 8163, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 1680, ckbuffer, 8523, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 3780, ckbuffer, 16953, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 4130, ckbuffer, 17553, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 4480, ckbuffer, 18153, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 61845, ckbuffer, 21174, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 61950, ckbuffer, 21354, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 62055, ckbuffer, 21534, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 62160, ckbuffer, 26556, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 62370, ckbuffer, 26916, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 62580, ckbuffer, 27276, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 62790, ckbuffer, 35706, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 63140, ckbuffer, 36306, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 63490, ckbuffer, 36906, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 63840, ckbuffer, 49611, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 64365, ckbuffer, 50511, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 64890, ckbuffer, 51411, 0, 4);

            t4cfunc::ket_transform<3, 2>(skbuffer, 65415, ckbuffer, 69258, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 66150, ckbuffer, 70518, 0, 5);

            t4cfunc::ket_transform<3, 2>(skbuffer, 66885, ckbuffer, 71778, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 12705, 0, 1260, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 13020, 105, 1470, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 13335, 210, 1680, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 16485, 1260, 3780, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 17115, 1470, 4130, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 17745, 1680, 4480, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 33495, 12705, 16485, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 34125, 13020, 17115, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 34755, 13335, 17745, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 315, 61845, 62160, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 1890, 62160, 62790, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 4830, 62790, 63840, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 7980, 63840, 65415, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 13650, 0, 315, 1890, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 18375, 1260, 1890, 4830, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 24045, 3780, 4830, 7980, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 35385, 12705, 13650, 18375, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 41055, 16485, 18375, 24045, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fpxx(skbuffer, 52395, 33495, 35385, 41055, r_ab, 3, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 52395, 3, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 735, skbuffer, 53445, 3, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1470, skbuffer, 54495, 3, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2205, skbuffer, 55545, 3, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2940, skbuffer, 56595, 3, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 3675, skbuffer, 57645, 3, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 4410, skbuffer, 58695, 3, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 5145, skbuffer, 59745, 3, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 5880, skbuffer, 60795, 3, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 3, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFPFD_hpp */
