#ifndef ElectronRepulsionGeom1010RecFDFP_hpp
#define ElectronRepulsionGeom1010RecFDFP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXDP.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXFP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom1010ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
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
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FD|1/|r-r'||FP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fdfp(T& distributor,
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

    CSimdArray<double> pbuffer(11747, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(8214, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(63270, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(69174, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 462, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 465, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 468, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 471, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 474, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 477, 1, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 486, 2, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 495, 3, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 504, 4, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 513, 5, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 522, 6, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 531, 15, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 549, 18, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 567, 21, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 585, 24, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 603, 27, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 621, 30, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 639, 51, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 669, 57, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 699, 63, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 729, 69, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 759, 75, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 789, 81, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 819, 115, 195, 210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 864, 125, 210, 225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 909, 135, 225, 240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 954, 145, 240, 255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 999, 155, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1044, 165, 270, 285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1089, 210, 315, 336, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1152, 225, 336, 357, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1215, 240, 357, 378, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1278, 255, 378, 399, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1341, 270, 399, 420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1404, 285, 420, 441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1467, 1, 2, 462, 465, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1473, 2, 3, 465, 468, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1479, 3, 4, 468, 471, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1485, 4, 5, 471, 474, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1491, 12, 15, 462, 477, 486, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1509, 15, 18, 465, 486, 495, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1527, 18, 21, 468, 495, 504, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1545, 21, 24, 471, 504, 513, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1563, 24, 27, 474, 513, 522, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1581, 45, 51, 486, 531, 549, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1617, 51, 57, 495, 549, 567, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1653, 57, 63, 504, 567, 585, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1689, 63, 69, 513, 585, 603, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1725, 69, 75, 522, 603, 621, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1761, 105, 115, 549, 639, 669, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1821, 115, 125, 567, 669, 699, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1881, 125, 135, 585, 699, 729, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1941, 135, 145, 603, 729, 759, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 2001, 145, 155, 621, 759, 789, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2061, 195, 210, 669, 819, 864, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2151, 210, 225, 699, 864, 909, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2241, 225, 240, 729, 909, 954, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2331, 240, 255, 759, 954, 999, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2421, 255, 270, 789, 999, 1044, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2511, 315, 336, 864, 1089, 1152, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2637, 336, 357, 909, 1152, 1215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2763, 357, 378, 954, 1215, 1278, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2889, 378, 399, 999, 1278, 1341, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 3015, 399, 420, 1044, 1341, 1404, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3141, 462, 465, 1467, 1473, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3151, 465, 468, 1473, 1479, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 3161, 468, 471, 1479, 1485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3171, 477, 486, 1467, 1491, 1509, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3201, 486, 495, 1473, 1509, 1527, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3231, 495, 504, 1479, 1527, 1545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 3261, 504, 513, 1485, 1545, 1563, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3291, 531, 549, 1509, 1581, 1617, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3351, 549, 567, 1527, 1617, 1653, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3411, 567, 585, 1545, 1653, 1689, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3471, 585, 603, 1563, 1689, 1725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3531, 639, 669, 1617, 1761, 1821, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3631, 669, 699, 1653, 1821, 1881, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3731, 699, 729, 1689, 1881, 1941, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3831, 729, 759, 1725, 1941, 2001, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3931, 819, 864, 1821, 2061, 2151, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4081, 864, 909, 1881, 2151, 2241, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4231, 909, 954, 1941, 2241, 2331, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4381, 954, 999, 2001, 2331, 2421, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4531, 1089, 1152, 2151, 2511, 2637, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4741, 1152, 1215, 2241, 2637, 2763, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 4951, 1215, 1278, 2331, 2763, 2889, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 5161, 1278, 1341, 2421, 2889, 3015, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 5371, 1467, 1473, 3141, 3151, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 5386, 1473, 1479, 3151, 3161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 5401, 1491, 1509, 3141, 3171, 3201, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 5446, 1509, 1527, 3151, 3201, 3231, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 5491, 1527, 1545, 3161, 3231, 3261, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 5536, 1581, 1617, 3201, 3291, 3351, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 5626, 1617, 1653, 3231, 3351, 3411, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 5716, 1653, 1689, 3261, 3411, 3471, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5806, 1761, 1821, 3351, 3531, 3631, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5956, 1821, 1881, 3411, 3631, 3731, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 6106, 1881, 1941, 3471, 3731, 3831, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6256, 2061, 2151, 3631, 3931, 4081, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6481, 2151, 2241, 3731, 4081, 4231, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6706, 2241, 2331, 3831, 4231, 4381, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 6931, 2511, 2637, 4081, 4531, 4741, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 7246, 2637, 2763, 4231, 4741, 4951, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 7561, 2763, 2889, 4381, 4951, 5161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 7876, 3141, 3151, 5371, 5386, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 7897, 3171, 3201, 5371, 5401, 5446, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 7960, 3201, 3231, 5386, 5446, 5491, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 8023, 3291, 3351, 5446, 5536, 5626, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 8149, 3351, 3411, 5491, 5626, 5716, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 8275, 3531, 3631, 5626, 5806, 5956, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 8485, 3631, 3731, 5716, 5956, 6106, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 8695, 3931, 4081, 5956, 6256, 6481, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 9010, 4081, 4231, 6106, 6481, 6706, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 9325, 4531, 4741, 6481, 6931, 7246, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 9766, 4741, 4951, 6706, 7246, 7561, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 10207, 5401, 5446, 7876, 7897, 7960, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 10291, 5536, 5626, 7960, 8023, 8149, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 10459, 5806, 5956, 8149, 8275, 8485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 10739, 6256, 6481, 8485, 8695, 9010, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sish(pbuffer, 11159, 6931, 7246, 9010, 9325, 9766, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1491, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 54, pbuffer, 1761, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 114, pbuffer, 3171, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 144, pbuffer, 3291, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 204, pbuffer, 3531, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 304, pbuffer, 5401, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 349, pbuffer, 5536, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 439, pbuffer, 5806, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1491, 1509});

                pbuffer.scale(2.0 * a_exp, {1581, 1617});

                pbuffer.scale(2.0 * a_exp, {1761, 1821});

                pbuffer.scale(2.0 * a_exp, {3171, 3201});

                pbuffer.scale(2.0 * a_exp, {3291, 3351});

                pbuffer.scale(2.0 * a_exp, {3531, 3631});

                pbuffer.scale(2.0 * a_exp, {5401, 5446});

                pbuffer.scale(2.0 * a_exp, {5536, 5626});

                pbuffer.scale(2.0 * a_exp, {5806, 5956});

                pbuffer.scale(2.0 * a_exp, {7897, 7960});

                pbuffer.scale(2.0 * a_exp, {8023, 8149});

                pbuffer.scale(2.0 * a_exp, {8275, 8485});

                pbuffer.scale(2.0 * a_exp, {10207, 10291});

                pbuffer.scale(2.0 * a_exp, {10291, 10459});

                pbuffer.scale(2.0 * a_exp, {10459, 10739});

                t2cfunc::reduce(cbuffer, 2294, pbuffer, 1491, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2312, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2348, pbuffer, 1761, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2408, pbuffer, 3171, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2438, pbuffer, 3291, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2498, pbuffer, 3531, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2598, pbuffer, 5401, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2643, pbuffer, 5536, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2733, pbuffer, 5806, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2883, pbuffer, 7897, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2946, pbuffer, 8023, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3072, pbuffer, 8275, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3282, pbuffer, 10207, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3366, pbuffer, 10291, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3534, pbuffer, 10459, 280, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1491, 1509});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1581, 1617});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1761, 1821});

                pbuffer.scale(pfactors, 0, 2.0, {2061, 2151});

                pbuffer.scale(pfactors, 0, 2.0, {2511, 2637});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3171, 3201});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3291, 3351});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3531, 3631});

                pbuffer.scale(pfactors, 0, 2.0, {3931, 4081});

                pbuffer.scale(pfactors, 0, 2.0, {4531, 4741});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {5401, 5446});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {5536, 5626});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {5806, 5956});

                pbuffer.scale(pfactors, 0, 2.0, {6256, 6481});

                pbuffer.scale(pfactors, 0, 2.0, {6931, 7246});

                t2cfunc::reduce(cbuffer, 589, pbuffer, 1491, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 607, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 643, pbuffer, 1761, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 703, pbuffer, 2061, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 793, pbuffer, 2511, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 919, pbuffer, 3171, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 949, pbuffer, 3291, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1009, pbuffer, 3531, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1109, pbuffer, 3931, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1259, pbuffer, 4531, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1469, pbuffer, 5401, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1514, pbuffer, 5536, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1604, pbuffer, 5806, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1754, pbuffer, 6256, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1979, pbuffer, 6931, 315, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1491, 1509});

                pbuffer.scale(2.0 * a_exp, {1581, 1617});

                pbuffer.scale(2.0 * a_exp, {1761, 1821});

                pbuffer.scale(2.0 * a_exp, {2061, 2151});

                pbuffer.scale(2.0 * a_exp, {2511, 2637});

                pbuffer.scale(2.0 * a_exp, {3171, 3201});

                pbuffer.scale(2.0 * a_exp, {3291, 3351});

                pbuffer.scale(2.0 * a_exp, {3531, 3631});

                pbuffer.scale(2.0 * a_exp, {3931, 4081});

                pbuffer.scale(2.0 * a_exp, {4531, 4741});

                pbuffer.scale(2.0 * a_exp, {5401, 5446});

                pbuffer.scale(2.0 * a_exp, {5536, 5626});

                pbuffer.scale(2.0 * a_exp, {5806, 5956});

                pbuffer.scale(2.0 * a_exp, {6256, 6481});

                pbuffer.scale(2.0 * a_exp, {6931, 7246});

                pbuffer.scale(pfactors, 0, 2.0, {7897, 7960});

                pbuffer.scale(pfactors, 0, 2.0, {8023, 8149});

                pbuffer.scale(pfactors, 0, 2.0, {8275, 8485});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {8695, 9010});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {9325, 9766});

                pbuffer.scale(pfactors, 0, 2.0, {10207, 10291});

                pbuffer.scale(pfactors, 0, 2.0, {10291, 10459});

                pbuffer.scale(pfactors, 0, 2.0, {10459, 10739});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {10739, 11159});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {11159, 11747});

                t2cfunc::reduce(cbuffer, 3814, pbuffer, 1491, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3832, pbuffer, 1581, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3868, pbuffer, 1761, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3928, pbuffer, 2061, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4018, pbuffer, 2511, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4144, pbuffer, 3171, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4174, pbuffer, 3291, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4234, pbuffer, 3531, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4334, pbuffer, 3931, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4484, pbuffer, 4531, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4694, pbuffer, 5401, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4739, pbuffer, 5536, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4829, pbuffer, 5806, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4979, pbuffer, 6256, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5204, pbuffer, 6931, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5519, pbuffer, 7897, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5582, pbuffer, 8023, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5708, pbuffer, 8275, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5918, pbuffer, 8695, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6233, pbuffer, 9325, 441, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6674, pbuffer, 10207, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6758, pbuffer, 10291, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6926, pbuffer, 10459, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7206, pbuffer, 10739, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7626, pbuffer, 11159, 588, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 612, cbuffer, 0, 18, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 828, cbuffer, 18, 54, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 1800, 612, 828, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4440, cbuffer, 114, 144, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4800, cbuffer, 144, 204, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 6420, 4440, 4800, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 10650, cbuffer, 304, 349, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 11190, cbuffer, 349, 439, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 13620, 10650, 11190, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 18282, cbuffer, 2294, 2312, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 18498, cbuffer, 2312, 2348, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 19470, 18282, 18498, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 22110, cbuffer, 2408, 2438, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 22470, cbuffer, 2438, 2498, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 24090, 22110, 22470, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 28320, cbuffer, 2598, 2643, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 28860, cbuffer, 2643, 2733, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 31290, 28320, 28860, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 37482, cbuffer, 2883, 2946, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 38238, cbuffer, 2946, 3072, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 41640, 37482, 38238, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 50166, cbuffer, 3282, 3366, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 51174, cbuffer, 3366, 3534, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 55710, 50166, 51174, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 589, 607, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 54, cbuffer, 607, 643, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 162, cbuffer, 643, 703, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 342, cbuffer, 703, 793, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 666, cbuffer, 0, 0, 54, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 936, cbuffer, 18, 54, 162, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 1260, cbuffer, 54, 162, 342, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1908, 612, 666, 936, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 2232, 828, 936, 1260, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 2880, 1800, 1908, 2232, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3420, cbuffer, 919, 949, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3510, cbuffer, 949, 1009, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3690, cbuffer, 1009, 1109, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 3990, cbuffer, 1109, 1259, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4530, cbuffer, 114, 3420, 3510, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4980, cbuffer, 144, 3510, 3690, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 5520, cbuffer, 204, 3690, 3990, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 6600, 4440, 4530, 4980, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 7140, 4800, 4980, 5520, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 8220, 6420, 6600, 7140, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9120, cbuffer, 1469, 1514, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9255, cbuffer, 1514, 1604, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 9525, cbuffer, 1604, 1754, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 9975, cbuffer, 1754, 1979, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 10785, cbuffer, 304, 9120, 9255, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 11460, cbuffer, 349, 9255, 9525, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 12270, cbuffer, 439, 9525, 9975, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 13890, 10650, 10785, 11460, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 14700, 11190, 11460, 12270, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 16320, 13620, 13890, 14700, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 17670, cbuffer, 3814, 3832, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 17724, cbuffer, 3832, 3868, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 17832, cbuffer, 3868, 3928, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 18012, cbuffer, 3928, 4018, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 18336, cbuffer, 2294, 17670, 17724, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 18606, cbuffer, 2312, 17724, 17832, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 18930, cbuffer, 2348, 17832, 18012, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 19578, 18282, 18336, 18606, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 19902, 18498, 18606, 18930, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 20550, 19470, 19578, 19902, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 21090, cbuffer, 4144, 4174, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 21180, cbuffer, 4174, 4234, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 21360, cbuffer, 4234, 4334, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 21660, cbuffer, 4334, 4484, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 22200, cbuffer, 2408, 21090, 21180, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 22650, cbuffer, 2438, 21180, 21360, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 23190, cbuffer, 2498, 21360, 21660, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 24270, 22110, 22200, 22650, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 24810, 22470, 22650, 23190, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 25890, 24090, 24270, 24810, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 26790, cbuffer, 4694, 4739, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 26925, cbuffer, 4739, 4829, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 27195, cbuffer, 4829, 4979, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 27645, cbuffer, 4979, 5204, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 28455, cbuffer, 2598, 26790, 26925, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 29130, cbuffer, 2643, 26925, 27195, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 29940, cbuffer, 2733, 27195, 27645, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 31560, 28320, 28455, 29130, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 32370, 28860, 29130, 29940, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 33990, 31290, 31560, 32370, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 35340, cbuffer, 5519, 5582, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 35529, cbuffer, 5582, 5708, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 35907, cbuffer, 5708, 5918, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 36537, cbuffer, 5918, 6233, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 37671, cbuffer, 2883, 35340, 35529, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 38616, cbuffer, 2946, 35529, 35907, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 39750, cbuffer, 3072, 35907, 36537, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 42018, 37482, 37671, 38616, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 43152, 38238, 38616, 39750, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 45420, 41640, 42018, 43152, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 47310, cbuffer, 6674, 6758, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 47562, cbuffer, 6758, 6926, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 48066, cbuffer, 6926, 7206, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 48906, cbuffer, 7206, 7626, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 50418, cbuffer, 3282, 47310, 47562, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 51678, cbuffer, 3366, 47562, 48066, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 53190, cbuffer, 3534, 48066, 48906, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 56214, 50166, 50418, 51678, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 57726, 51174, 51678, 53190, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 60750, 55710, 56214, 57726, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<3, 1>(skbuffer, 0, ckbuffer, 2880, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 126, ckbuffer, 3060, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 252, ckbuffer, 3240, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 1512, ckbuffer, 8220, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 1722, ckbuffer, 8520, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 1932, ckbuffer, 8820, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 4032, ckbuffer, 16320, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 4347, ckbuffer, 16770, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 4662, ckbuffer, 17220, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 64134, ckbuffer, 20550, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 64260, ckbuffer, 20730, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 64386, ckbuffer, 20910, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 64512, ckbuffer, 25890, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 64722, ckbuffer, 26190, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 64932, ckbuffer, 26490, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 65142, ckbuffer, 33990, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 65457, ckbuffer, 34440, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 65772, ckbuffer, 34890, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 66087, ckbuffer, 45420, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 66528, ckbuffer, 46050, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 66969, ckbuffer, 46680, 0, 5);

            t4cfunc::ket_transform<3, 1>(skbuffer, 67410, ckbuffer, 60750, 0, 6);

            t4cfunc::ket_transform<3, 1>(skbuffer, 67998, ckbuffer, 61590, 0, 6);

            t4cfunc::ket_transform<3, 1>(skbuffer, 68586, ckbuffer, 62430, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 11781, 0, 1512, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 12159, 126, 1722, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 12537, 252, 1932, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 16317, 1512, 4032, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 16947, 1722, 4347, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 17577, 1932, 4662, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 32382, 11781, 16317, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 33138, 12159, 16947, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 33894, 12537, 17577, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 378, 64134, 64512, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 2142, 64512, 65142, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 4977, 65142, 66087, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 7812, 66087, 67410, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 12915, 0, 378, 2142, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 18207, 1512, 2142, 4977, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 23877, 4032, 4977, 7812, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 34650, 11781, 12915, 18207, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 41454, 16317, 18207, 23877, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fdxx(skbuffer, 52794, 32382, 34650, 41454, r_ab, 3, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 52794, 3, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 735, skbuffer, 54054, 3, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1470, skbuffer, 55314, 3, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2205, skbuffer, 56574, 3, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2940, skbuffer, 57834, 3, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 3675, skbuffer, 59094, 3, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 4410, skbuffer, 60354, 3, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 5145, skbuffer, 61614, 3, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 5880, skbuffer, 62874, 3, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 3, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFDFP_hpp */
