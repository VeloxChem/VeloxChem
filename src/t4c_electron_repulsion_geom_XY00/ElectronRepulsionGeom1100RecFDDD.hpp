#ifndef ElectronRepulsionGeom1100RecFDDD_hpp
#define ElectronRepulsionGeom1100RecFDDD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecXXDD.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0100ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSIXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecFDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSHXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSF.hpp"
#include "ElectronRepulsionPrimRecSKSG.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FD|1/|r-r'||DD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_fddd(T& distributor,
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

    CSimdArray<double> pbuffer(11026, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(7130, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(31308, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(97175, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 315, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 318, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 321, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 324, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 327, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 330, 2, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 339, 3, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 348, 4, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 357, 5, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 366, 6, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 375, 7, 30, 33, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 384, 15, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 402, 18, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 420, 21, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 438, 24, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 456, 27, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 474, 30, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 492, 33, 81, 87, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 510, 51, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 540, 57, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 570, 63, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 600, 69, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 630, 75, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 660, 81, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 690, 87, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 720, 115, 195, 210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 765, 125, 210, 225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 810, 135, 225, 240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 855, 145, 240, 255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 900, 155, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 945, 165, 270, 285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 990, 175, 285, 300, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1035, 2, 3, 315, 318, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1041, 3, 4, 318, 321, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1047, 4, 5, 321, 324, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1053, 5, 6, 324, 327, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1059, 15, 18, 315, 330, 339, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1077, 18, 21, 318, 339, 348, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1095, 21, 24, 321, 348, 357, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1113, 24, 27, 324, 357, 366, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1131, 27, 30, 327, 366, 375, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1149, 45, 51, 330, 384, 402, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1185, 51, 57, 339, 402, 420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1221, 57, 63, 348, 420, 438, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1257, 63, 69, 357, 438, 456, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1293, 69, 75, 366, 456, 474, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1329, 75, 81, 375, 474, 492, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1365, 105, 115, 402, 510, 540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1425, 115, 125, 420, 540, 570, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1485, 125, 135, 438, 570, 600, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1545, 135, 145, 456, 600, 630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1605, 145, 155, 474, 630, 660, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1665, 155, 165, 492, 660, 690, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1725, 195, 210, 540, 720, 765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1815, 210, 225, 570, 765, 810, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1905, 225, 240, 600, 810, 855, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1995, 240, 255, 630, 855, 900, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2085, 255, 270, 660, 900, 945, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2175, 270, 285, 690, 945, 990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2265, 315, 318, 1035, 1041, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2275, 318, 321, 1041, 1047, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2285, 321, 324, 1047, 1053, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2295, 330, 339, 1035, 1059, 1077, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2325, 339, 348, 1041, 1077, 1095, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2355, 348, 357, 1047, 1095, 1113, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2385, 357, 366, 1053, 1113, 1131, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2415, 384, 402, 1059, 1149, 1185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2475, 402, 420, 1077, 1185, 1221, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2535, 420, 438, 1095, 1221, 1257, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2595, 438, 456, 1113, 1257, 1293, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2655, 456, 474, 1131, 1293, 1329, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2715, 510, 540, 1185, 1365, 1425, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2815, 540, 570, 1221, 1425, 1485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2915, 570, 600, 1257, 1485, 1545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3015, 600, 630, 1293, 1545, 1605, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3115, 630, 660, 1329, 1605, 1665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3215, 720, 765, 1425, 1725, 1815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3365, 765, 810, 1485, 1815, 1905, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3515, 810, 855, 1545, 1905, 1995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3665, 855, 900, 1605, 1995, 2085, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3815, 900, 945, 1665, 2085, 2175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 3965, 1035, 1041, 2265, 2275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 3980, 1041, 1047, 2275, 2285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3995, 1059, 1077, 2265, 2295, 2325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4040, 1077, 1095, 2275, 2325, 2355, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4085, 1095, 1113, 2285, 2355, 2385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4130, 1149, 1185, 2295, 2415, 2475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4220, 1185, 1221, 2325, 2475, 2535, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4310, 1221, 1257, 2355, 2535, 2595, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4400, 1257, 1293, 2385, 2595, 2655, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4490, 1365, 1425, 2475, 2715, 2815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4640, 1425, 1485, 2535, 2815, 2915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4790, 1485, 1545, 2595, 2915, 3015, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4940, 1545, 1605, 2655, 3015, 3115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5090, 1725, 1815, 2815, 3215, 3365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5315, 1815, 1905, 2915, 3365, 3515, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5540, 1905, 1995, 3015, 3515, 3665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5765, 1995, 2085, 3115, 3665, 3815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 5990, 2265, 2275, 3965, 3980, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 6011, 2295, 2325, 3965, 3995, 4040, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 6074, 2325, 2355, 3980, 4040, 4085, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 6137, 2415, 2475, 3995, 4130, 4220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 6263, 2475, 2535, 4040, 4220, 4310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 6389, 2535, 2595, 4085, 4310, 4400, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 6515, 2715, 2815, 4220, 4490, 4640, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 6725, 2815, 2915, 4310, 4640, 4790, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 6935, 2915, 3015, 4400, 4790, 4940, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 7145, 3215, 3365, 4640, 5090, 5315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 7460, 3365, 3515, 4790, 5315, 5540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 7775, 3515, 3665, 4940, 5540, 5765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 8090, 3995, 4040, 5990, 6011, 6074, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 8174, 4130, 4220, 6011, 6137, 6263, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 8342, 4220, 4310, 6074, 6263, 6389, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 8510, 4490, 4640, 6263, 6515, 6725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 8790, 4640, 4790, 6389, 6725, 6935, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 9070, 5090, 5315, 6725, 7145, 7460, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 9490, 5315, 5540, 6935, 7460, 7775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 9910, 6137, 6263, 8090, 8174, 8342, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 10126, 6515, 6725, 8342, 8510, 8790, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 10486, 7145, 7460, 8790, 9070, 9490, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 384, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 510, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 48, pbuffer, 720, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 93, pbuffer, 1149, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 129, pbuffer, 1365, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 189, pbuffer, 1725, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 279, pbuffer, 2415, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 339, pbuffer, 2715, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 439, pbuffer, 3215, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {2415, 2475});

                pbuffer.scale(2.0 * b_exp, {2715, 2815});

                pbuffer.scale(2.0 * b_exp, {3215, 3365});

                pbuffer.scale(2.0 * b_exp, {4130, 4220});

                pbuffer.scale(2.0 * b_exp, {4490, 4640});

                pbuffer.scale(2.0 * b_exp, {5090, 5315});

                pbuffer.scale(2.0 * b_exp, {6137, 6263});

                pbuffer.scale(2.0 * b_exp, {6515, 6725});

                pbuffer.scale(2.0 * b_exp, {7145, 7460});

                t2cfunc::reduce(cbuffer, 589, pbuffer, 2415, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 649, pbuffer, 2715, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 749, pbuffer, 3215, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 899, pbuffer, 4130, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 989, pbuffer, 4490, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1139, pbuffer, 5090, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1364, pbuffer, 6137, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1490, pbuffer, 6515, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1700, pbuffer, 7145, 315, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {384, 402});

                pbuffer.scale(2.0 * a_exp, {510, 540});

                pbuffer.scale(2.0 * a_exp, {720, 765});

                pbuffer.scale(2.0 * a_exp, {1149, 1185});

                pbuffer.scale(2.0 * a_exp, {1365, 1425});

                pbuffer.scale(2.0 * a_exp, {1725, 1815});

                pbuffer.scale(a_exp / b_exp, {2415, 2475});

                pbuffer.scale(a_exp / b_exp, {2715, 2815});

                pbuffer.scale(a_exp / b_exp, {3215, 3365});

                pbuffer.scale(a_exp / b_exp, {4130, 4220});

                pbuffer.scale(a_exp / b_exp, {4490, 4640});

                pbuffer.scale(a_exp / b_exp, {5090, 5315});

                pbuffer.scale(a_exp / b_exp, {6137, 6263});

                pbuffer.scale(a_exp / b_exp, {6515, 6725});

                pbuffer.scale(a_exp / b_exp, {7145, 7460});

                t2cfunc::reduce(cbuffer, 2015, pbuffer, 384, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2033, pbuffer, 510, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2063, pbuffer, 720, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2108, pbuffer, 1149, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2144, pbuffer, 1365, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2204, pbuffer, 1725, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2294, pbuffer, 2415, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2354, pbuffer, 2715, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2454, pbuffer, 3215, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2604, pbuffer, 4130, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2694, pbuffer, 4490, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2844, pbuffer, 5090, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3069, pbuffer, 6137, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3195, pbuffer, 6515, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3405, pbuffer, 7145, 315, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {2415, 2475});

                pbuffer.scale(2.0 * b_exp, {2715, 2815});

                pbuffer.scale(2.0 * b_exp, {3215, 3365});

                pbuffer.scale(2.0 * b_exp, {4130, 4220});

                pbuffer.scale(2.0 * b_exp, {4490, 4640});

                pbuffer.scale(2.0 * b_exp, {5090, 5315});

                pbuffer.scale(2.0 * b_exp, {6137, 6263});

                pbuffer.scale(2.0 * b_exp, {6515, 6725});

                pbuffer.scale(2.0 * b_exp, {7145, 7460});

                pbuffer.scale(4.0 * a_exp * b_exp, {8174, 8342});

                pbuffer.scale(4.0 * a_exp * b_exp, {8510, 8790});

                pbuffer.scale(4.0 * a_exp * b_exp, {9070, 9490});

                pbuffer.scale(4.0 * a_exp * b_exp, {9910, 10126});

                pbuffer.scale(4.0 * a_exp * b_exp, {10126, 10486});

                pbuffer.scale(4.0 * a_exp * b_exp, {10486, 11026});

                t2cfunc::reduce(cbuffer, 3720, pbuffer, 2415, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3780, pbuffer, 2715, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3880, pbuffer, 3215, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4030, pbuffer, 4130, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4120, pbuffer, 4490, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4270, pbuffer, 5090, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4495, pbuffer, 6137, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4621, pbuffer, 6515, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4831, pbuffer, 7145, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5146, pbuffer, 8174, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5314, pbuffer, 8510, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5594, pbuffer, 9070, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6014, pbuffer, 9910, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6230, pbuffer, 10126, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6590, pbuffer, 10486, 540, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 18, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 54, cbuffer, 18, 48, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 144, 0, 54, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 252, cbuffer, 93, 129, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 360, cbuffer, 129, 189, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 540, 252, 360, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1404, cbuffer, 279, 339, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1584, cbuffer, 339, 439, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 1884, 1404, 1584, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4944, cbuffer, 589, 649, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5124, cbuffer, 649, 749, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 5424, 4944, 5124, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5784, cbuffer, 899, 989, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 6054, cbuffer, 989, 1139, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 6504, 5784, 6054, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 7044, cbuffer, 1364, 1490, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 7422, cbuffer, 1490, 1700, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 8052, 7044, 7422, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8808, cbuffer, 2015, 2033, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 8862, cbuffer, 2033, 2063, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 8952, 8808, 8862, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 9060, cbuffer, 2108, 2144, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9168, cbuffer, 2144, 2204, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 9348, 9060, 9168, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 10212, cbuffer, 2294, 2354, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 10392, cbuffer, 2354, 2454, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 10692, 10212, 10392, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 12132, cbuffer, 2604, 2694, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 12402, cbuffer, 2694, 2844, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 12852, 12132, 12402, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 15012, cbuffer, 3069, 3195, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 15390, cbuffer, 3195, 3405, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 16020, 15012, 15390, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 22068, cbuffer, 3720, 3780, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 22248, cbuffer, 3780, 3880, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 22548, 22068, 22248, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 22908, cbuffer, 4030, 4120, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 23178, cbuffer, 4120, 4270, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 23628, 22908, 23178, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 24168, cbuffer, 4495, 4621, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 24546, cbuffer, 4621, 4831, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 25176, 24168, 24546, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 25932, cbuffer, 5146, 5314, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 26436, cbuffer, 5314, 5594, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 27276, 25932, 26436, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 28284, cbuffer, 6014, 6230, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 28932, cbuffer, 6230, 6590, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 30012, 28284, 28932, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<2, 2>(skbuffer, 0, ckbuffer, 144, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 75, ckbuffer, 540, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 2025, ckbuffer, 1884, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 83575, ckbuffer, 5424, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 83825, ckbuffer, 6504, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 84200, ckbuffer, 8052, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 84725, ckbuffer, 8952, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 84800, ckbuffer, 9348, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 85400, ckbuffer, 10692, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 86400, ckbuffer, 12852, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 87900, ckbuffer, 16020, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 94425, ckbuffer, 22548, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 94675, ckbuffer, 23628, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 95050, ckbuffer, 25176, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 95575, ckbuffer, 27276, 0, 6);

            t4cfunc::ket_transform<2, 2>(skbuffer, 96275, ckbuffer, 30012, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 14500, 75, 2025, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 92100, 84800, 85400, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 92550, 85400, 86400, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 93300, 86400, 87900, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 16300, 75, 92100, 92550, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 23950, 2025, 92550, 93300, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ddxx(skbuffer, 45775, 14500, 16300, 23950, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 225, 0, 83575, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 2275, 75, 83825, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 5275, 2025, 84200, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 14950, 75, 225, 2275, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 21700, 2025, 2275, 5275, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ddxx(skbuffer, 43075, 14500, 14950, 21700, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 84950, 84725, 94425, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 85650, 84800, 94675, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 86775, 85400, 95050, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 88425, 86400, 95575, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 90000, 87900, 96275, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 675, 84800, 84950, 85650, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 3025, 85400, 85650, 86775, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 6400, 86400, 86775, 88425, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 9775, 87900, 88425, 90000, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 17650, 225, 92100, 675, 3025, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 26200, 2275, 92550, 3025, 6400, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 32950, 5275, 93300, 6400, 9775, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 48475, 14950, 16300, 17650, 26200, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 56575, 21700, 23950, 26200, 32950, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_fdxx(skbuffer, 70075, 43075, 45775, 48475, 56575, r_ab, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 70075, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 875, skbuffer, 71575, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1750, skbuffer, 73075, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2625, skbuffer, 74575, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 3500, skbuffer, 76075, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 4375, skbuffer, 77575, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 5250, skbuffer, 79075, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 6125, skbuffer, 80575, 2, 2);

            t4cfunc::bra_transform<3, 2>(sbuffer, 7000, skbuffer, 82075, 2, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 2, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFDDD_hpp */
