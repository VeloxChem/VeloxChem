#ifndef ElectronRepulsionGeom1100RecFFPF_hpp
#define ElectronRepulsionGeom1100RecFFPF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0100ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSIXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSKXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecFFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPHXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSIXX.hpp"
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
#include "ElectronRepulsionPrimRecSLSF.hpp"
#include "ElectronRepulsionPrimRecSLSG.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FF|1/|r-r'||PF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_ffpf(T& distributor,
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

    CSimdArray<double> pbuffer(15161, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(8000, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(23640, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(124404, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(9261, 1);

    // setup Boys fuction data

    const CBoysFunc<12> bf_table;

    CSimdArray<double> bf_data(14, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 13, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 13, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 13, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 13);
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

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 12, pfactors, 16, bf_data, 12);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 25, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 28, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 31, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 34, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 37, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 40, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 43, 10, 11, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 46, 11, 12, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 49, 0, 1, 13, 16, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 55, 1, 2, 16, 19, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 61, 2, 3, 19, 22, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 67, 3, 4, 22, 25, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 73, 4, 5, 25, 28, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 79, 5, 6, 28, 31, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 85, 6, 7, 31, 34, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 91, 7, 8, 34, 37, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 97, 8, 9, 37, 40, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 103, 9, 10, 40, 43, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 109, 10, 11, 43, 46, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 13, 16, 49, 55, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 16, 19, 55, 61, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 19, 22, 61, 67, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 22, 25, 67, 73, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 25, 28, 73, 79, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 28, 31, 79, 85, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 175, 31, 34, 85, 91, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 185, 34, 37, 91, 97, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 195, 37, 40, 97, 103, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 205, 40, 43, 103, 109, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 215, 49, 55, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 230, 55, 61, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 245, 61, 67, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 260, 67, 73, 145, 155, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 275, 73, 79, 155, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 290, 79, 85, 165, 175, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 305, 85, 91, 175, 185, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 320, 91, 97, 185, 195, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 335, 97, 103, 195, 205, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 350, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 353, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 356, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 359, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 362, 7, 8, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 365, 3, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 374, 4, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 383, 5, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 392, 6, 28, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 401, 7, 31, 34, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 410, 8, 34, 37, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 419, 19, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 437, 22, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 455, 25, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 473, 28, 73, 79, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 491, 31, 79, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 509, 34, 85, 91, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 527, 37, 91, 97, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 545, 55, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 575, 61, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 605, 67, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 635, 73, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 665, 79, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 695, 85, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 725, 91, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 755, 97, 185, 195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 785, 125, 215, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 830, 135, 230, 245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 875, 145, 245, 260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 920, 155, 260, 275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 965, 165, 275, 290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1010, 175, 290, 305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1055, 185, 305, 320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1100, 195, 320, 335, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1145, 3, 4, 350, 353, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1151, 4, 5, 353, 356, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1157, 5, 6, 356, 359, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1163, 6, 7, 359, 362, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1169, 19, 22, 350, 365, 374, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1187, 22, 25, 353, 374, 383, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1205, 25, 28, 356, 383, 392, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1223, 28, 31, 359, 392, 401, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1241, 31, 34, 362, 401, 410, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1259, 55, 61, 365, 419, 437, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1295, 61, 67, 374, 437, 455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1331, 67, 73, 383, 455, 473, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1367, 73, 79, 392, 473, 491, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1403, 79, 85, 401, 491, 509, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1439, 85, 91, 410, 509, 527, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1475, 115, 125, 419, 545, 575, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1535, 125, 135, 437, 575, 605, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1595, 135, 145, 455, 605, 635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1655, 145, 155, 473, 635, 665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1715, 155, 165, 491, 665, 695, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1775, 165, 175, 509, 695, 725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1835, 175, 185, 527, 725, 755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1895, 215, 230, 575, 785, 830, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1985, 230, 245, 605, 830, 875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2075, 245, 260, 635, 875, 920, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2165, 260, 275, 665, 920, 965, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2255, 275, 290, 695, 965, 1010, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2345, 290, 305, 725, 1010, 1055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2435, 305, 320, 755, 1055, 1100, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2525, 350, 353, 1145, 1151, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2535, 353, 356, 1151, 1157, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2545, 356, 359, 1157, 1163, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2555, 365, 374, 1145, 1169, 1187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2585, 374, 383, 1151, 1187, 1205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2615, 383, 392, 1157, 1205, 1223, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2645, 392, 401, 1163, 1223, 1241, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2675, 419, 437, 1169, 1259, 1295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2735, 437, 455, 1187, 1295, 1331, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2795, 455, 473, 1205, 1331, 1367, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2855, 473, 491, 1223, 1367, 1403, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2915, 491, 509, 1241, 1403, 1439, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2975, 545, 575, 1259, 1475, 1535, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3075, 575, 605, 1295, 1535, 1595, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3175, 605, 635, 1331, 1595, 1655, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3275, 635, 665, 1367, 1655, 1715, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3375, 665, 695, 1403, 1715, 1775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3475, 695, 725, 1439, 1775, 1835, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3575, 785, 830, 1535, 1895, 1985, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3725, 830, 875, 1595, 1985, 2075, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3875, 875, 920, 1655, 2075, 2165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4025, 920, 965, 1715, 2165, 2255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4175, 965, 1010, 1775, 2255, 2345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4325, 1010, 1055, 1835, 2345, 2435, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4475, 1145, 1151, 2525, 2535, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4490, 1151, 1157, 2535, 2545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4505, 1169, 1187, 2525, 2555, 2585, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4550, 1187, 1205, 2535, 2585, 2615, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4595, 1205, 1223, 2545, 2615, 2645, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4640, 1259, 1295, 2555, 2675, 2735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4730, 1295, 1331, 2585, 2735, 2795, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4820, 1331, 1367, 2615, 2795, 2855, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4910, 1367, 1403, 2645, 2855, 2915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5000, 1475, 1535, 2675, 2975, 3075, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5150, 1535, 1595, 2735, 3075, 3175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5300, 1595, 1655, 2795, 3175, 3275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5450, 1655, 1715, 2855, 3275, 3375, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5600, 1715, 1775, 2915, 3375, 3475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5750, 1895, 1985, 3075, 3575, 3725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5975, 1985, 2075, 3175, 3725, 3875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6200, 2075, 2165, 3275, 3875, 4025, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6425, 2165, 2255, 3375, 4025, 4175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6650, 2255, 2345, 3475, 4175, 4325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 6875, 2525, 2535, 4475, 4490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 6896, 2555, 2585, 4475, 4505, 4550, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 6959, 2585, 2615, 4490, 4550, 4595, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 7022, 2675, 2735, 4505, 4640, 4730, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 7148, 2735, 2795, 4550, 4730, 4820, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 7274, 2795, 2855, 4595, 4820, 4910, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 7400, 2975, 3075, 4640, 5000, 5150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 7610, 3075, 3175, 4730, 5150, 5300, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 7820, 3175, 3275, 4820, 5300, 5450, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 8030, 3275, 3375, 4910, 5450, 5600, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 8240, 3575, 3725, 5150, 5750, 5975, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 8555, 3725, 3875, 5300, 5975, 6200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 8870, 3875, 4025, 5450, 6200, 6425, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 9185, 4025, 4175, 5600, 6425, 6650, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 9500, 4505, 4550, 6875, 6896, 6959, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 9584, 4640, 4730, 6896, 7022, 7148, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 9752, 4730, 4820, 6959, 7148, 7274, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 9920, 5000, 5150, 7022, 7400, 7610, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 10200, 5150, 5300, 7148, 7610, 7820, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 10480, 5300, 5450, 7274, 7820, 8030, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 10760, 5750, 5975, 7610, 8240, 8555, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 11180, 5975, 6200, 7820, 8555, 8870, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 11600, 6200, 6425, 8030, 8870, 9185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 12020, 7022, 7148, 9500, 9584, 9752, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 12236, 7400, 7610, 9584, 9920, 10200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 12596, 7610, 7820, 9752, 10200, 10480, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 12956, 8240, 8555, 10200, 10760, 11180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 13496, 8555, 8870, 10480, 11180, 11600, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsf(pbuffer, 14036, 9920, 10200, 12020, 12236, 12596, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsg(pbuffer, 14486, 10760, 11180, 12596, 12956, 13496, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1475, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 1895, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 150, pbuffer, 2975, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 250, pbuffer, 3575, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 400, pbuffer, 5000, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 550, pbuffer, 5750, 225, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {5000, 5150});

                pbuffer.scale(2.0 * b_exp, {5750, 5975});

                pbuffer.scale(2.0 * b_exp, {7400, 7610});

                pbuffer.scale(2.0 * b_exp, {8240, 8555});

                pbuffer.scale(2.0 * b_exp, {9920, 10200});

                pbuffer.scale(2.0 * b_exp, {10760, 11180});

                t2cfunc::reduce(cbuffer, 775, pbuffer, 5000, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 925, pbuffer, 5750, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1150, pbuffer, 7400, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1360, pbuffer, 8240, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1675, pbuffer, 9920, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1955, pbuffer, 10760, 420, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1475, 1535});

                pbuffer.scale(2.0 * a_exp, {1895, 1985});

                pbuffer.scale(2.0 * a_exp, {2975, 3075});

                pbuffer.scale(2.0 * a_exp, {3575, 3725});

                pbuffer.scale(a_exp / b_exp, {5000, 5150});

                pbuffer.scale(a_exp / b_exp, {5750, 5975});

                pbuffer.scale(a_exp / b_exp, {7400, 7610});

                pbuffer.scale(a_exp / b_exp, {8240, 8555});

                pbuffer.scale(a_exp / b_exp, {9920, 10200});

                pbuffer.scale(a_exp / b_exp, {10760, 11180});

                t2cfunc::reduce(cbuffer, 2375, pbuffer, 1475, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2435, pbuffer, 1895, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2525, pbuffer, 2975, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2625, pbuffer, 3575, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2775, pbuffer, 5000, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2925, pbuffer, 5750, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3150, pbuffer, 7400, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3360, pbuffer, 8240, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3675, pbuffer, 9920, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3955, pbuffer, 10760, 420, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {5000, 5150});

                pbuffer.scale(2.0 * b_exp, {5750, 5975});

                pbuffer.scale(2.0 * b_exp, {7400, 7610});

                pbuffer.scale(2.0 * b_exp, {8240, 8555});

                pbuffer.scale(2.0 * b_exp, {9920, 10200});

                pbuffer.scale(2.0 * b_exp, {10760, 11180});

                pbuffer.scale(4.0 * a_exp * b_exp, {12236, 12596});

                pbuffer.scale(4.0 * a_exp * b_exp, {12956, 13496});

                pbuffer.scale(4.0 * a_exp * b_exp, {14036, 14486});

                pbuffer.scale(4.0 * a_exp * b_exp, {14486, 15161});

                t2cfunc::reduce(cbuffer, 4375, pbuffer, 5000, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4525, pbuffer, 5750, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4750, pbuffer, 7400, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4960, pbuffer, 8240, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5275, pbuffer, 9920, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5555, pbuffer, 10760, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5975, pbuffer, 12236, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6335, pbuffer, 12956, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6875, pbuffer, 14036, 450, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7325, pbuffer, 14486, 675, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 60, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 180, cbuffer, 150, 250, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1380, cbuffer, 400, 550, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5070, cbuffer, 775, 925, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5520, cbuffer, 1150, 1360, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 6150, cbuffer, 1675, 1955, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 6990, cbuffer, 2375, 2435, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 7170, cbuffer, 2525, 2625, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 8370, cbuffer, 2775, 2925, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 10170, cbuffer, 3150, 3360, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 12690, cbuffer, 3675, 3955, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 19290, cbuffer, 4375, 4525, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 19740, cbuffer, 4750, 4960, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 20370, cbuffer, 5275, 5555, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 21210, cbuffer, 5975, 6335, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 22290, cbuffer, 6875, 7325, cfactors, 6, 0, 8);

            t4cfunc::ket_transform<1, 3>(skbuffer, 0, ckbuffer, 0, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 126, ckbuffer, 180, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 2856, ckbuffer, 1380, 0, 4);

            t4cfunc::ket_transform<1, 3>(skbuffer, 108507, ckbuffer, 5070, 0, 4);

            t4cfunc::ket_transform<1, 3>(skbuffer, 108822, ckbuffer, 5520, 0, 5);

            t4cfunc::ket_transform<1, 3>(skbuffer, 109263, ckbuffer, 6150, 0, 6);

            t4cfunc::ket_transform<1, 3>(skbuffer, 109851, ckbuffer, 6990, 0, 2);

            t4cfunc::ket_transform<1, 3>(skbuffer, 109977, ckbuffer, 7170, 0, 3);

            t4cfunc::ket_transform<1, 3>(skbuffer, 110817, ckbuffer, 8370, 0, 4);

            t4cfunc::ket_transform<1, 3>(skbuffer, 112077, ckbuffer, 10170, 0, 5);

            t4cfunc::ket_transform<1, 3>(skbuffer, 113841, ckbuffer, 12690, 0, 6);

            t4cfunc::ket_transform<1, 3>(skbuffer, 121359, ckbuffer, 19290, 0, 4);

            t4cfunc::ket_transform<1, 3>(skbuffer, 121674, ckbuffer, 19740, 0, 5);

            t4cfunc::ket_transform<1, 3>(skbuffer, 122115, ckbuffer, 20370, 0, 6);

            t4cfunc::ket_transform<1, 3>(skbuffer, 122703, ckbuffer, 21210, 0, 7);

            t4cfunc::ket_transform<1, 3>(skbuffer, 123459, ckbuffer, 22290, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 17535, 126, 2856, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 118461, 109977, 110817, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 119091, 110817, 112077, r_ab, 1, 3);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 120036, 112077, 113841, r_ab, 1, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 20055, 126, 118461, 119091, r_ab, 1, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 30450, 2856, 119091, 120036, r_ab, 1, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 57477, 17535, 20055, 30450, r_ab, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 336, 0, 108507, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 3171, 126, 108822, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 6951, 2856, 109263, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 18165, 126, 336, 3171, r_ab, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pgxx(skbuffer, 27615, 2856, 3171, 6951, r_ab, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dfxx(skbuffer, 53697, 17535, 18165, 27615, r_ab, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 110187, 109851, 121359, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 111132, 109977, 121674, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 112518, 110817, 122115, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 114429, 112077, 122703, 1, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_skxx(skbuffer, 116193, 113841, 123459, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 966, 109977, 110187, 111132, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 4116, 110817, 111132, 112518, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 8274, 112077, 112518, 114429, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sixx(skbuffer, 12243, 113841, 114429, 116193, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 21945, 336, 118461, 966, 4116, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 33285, 3171, 119091, 4116, 8274, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_phxx(skbuffer, 41790, 6951, 120036, 8274, 12243, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 61257, 18165, 20055, 21945, 33285, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dgxx(skbuffer, 72597, 27615, 30450, 33285, 41790, r_ab, 1, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ffxx(skbuffer, 89607, 53697, 57477, 61257, 72597, r_ab, 1, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 89607, 1, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1029, skbuffer, 91707, 1, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2058, skbuffer, 93807, 1, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 3087, skbuffer, 95907, 1, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 4116, skbuffer, 98007, 1, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 5145, skbuffer, 100107, 1, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 6174, skbuffer, 102207, 1, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 7203, skbuffer, 104307, 1, 3);

            t4cfunc::bra_transform<3, 3>(sbuffer, 8232, skbuffer, 106407, 1, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 1, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFFPF_hpp */
