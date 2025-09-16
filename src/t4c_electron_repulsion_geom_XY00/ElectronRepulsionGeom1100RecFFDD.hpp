#ifndef ElectronRepulsionGeom1100RecFFDD_hpp
#define ElectronRepulsionGeom1100RecFFDD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecXXDD.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
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
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSF.hpp"
#include "ElectronRepulsionPrimRecSKSG.hpp"
#include "ElectronRepulsionPrimRecSKSP.hpp"
#include "ElectronRepulsionPrimRecSLSD.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FF|1/|r-r'||DD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_ffdd(T& distributor,
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

    CSimdArray<double> pbuffer(16585, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(9920, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(43728, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(148100, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(11025, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 350, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 353, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 356, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 359, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 362, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 365, 7, 8, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 368, 2, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 377, 3, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 386, 4, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 395, 5, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 404, 6, 28, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 413, 7, 31, 34, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 422, 8, 34, 37, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 431, 16, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 449, 19, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 467, 22, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 485, 25, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 503, 28, 73, 79, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 521, 31, 79, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 539, 34, 85, 91, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 557, 37, 91, 97, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 575, 55, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 605, 61, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 635, 67, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 665, 73, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 695, 79, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 725, 85, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 755, 91, 175, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 785, 97, 185, 195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 815, 125, 215, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 860, 135, 230, 245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 905, 145, 245, 260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 950, 155, 260, 275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 995, 165, 275, 290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1040, 175, 290, 305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1085, 185, 305, 320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1130, 195, 320, 335, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1175, 2, 3, 350, 353, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1181, 3, 4, 353, 356, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1187, 4, 5, 356, 359, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1193, 5, 6, 359, 362, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1199, 6, 7, 362, 365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1205, 16, 19, 350, 368, 377, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1223, 19, 22, 353, 377, 386, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1241, 22, 25, 356, 386, 395, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1259, 25, 28, 359, 395, 404, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1277, 28, 31, 362, 404, 413, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1295, 31, 34, 365, 413, 422, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1313, 49, 55, 368, 431, 449, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1349, 55, 61, 377, 449, 467, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1385, 61, 67, 386, 467, 485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1421, 67, 73, 395, 485, 503, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1457, 73, 79, 404, 503, 521, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1493, 79, 85, 413, 521, 539, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1529, 85, 91, 422, 539, 557, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1565, 115, 125, 449, 575, 605, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1625, 125, 135, 467, 605, 635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1685, 135, 145, 485, 635, 665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1745, 145, 155, 503, 665, 695, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1805, 155, 165, 521, 695, 725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1865, 165, 175, 539, 725, 755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1925, 175, 185, 557, 755, 785, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1985, 215, 230, 605, 815, 860, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2075, 230, 245, 635, 860, 905, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2165, 245, 260, 665, 905, 950, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2255, 260, 275, 695, 950, 995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2345, 275, 290, 725, 995, 1040, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2435, 290, 305, 755, 1040, 1085, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2525, 305, 320, 785, 1085, 1130, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2615, 350, 353, 1175, 1181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2625, 353, 356, 1181, 1187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2635, 356, 359, 1187, 1193, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2645, 359, 362, 1193, 1199, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2655, 368, 377, 1175, 1205, 1223, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2685, 377, 386, 1181, 1223, 1241, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2715, 386, 395, 1187, 1241, 1259, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2745, 395, 404, 1193, 1259, 1277, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2775, 404, 413, 1199, 1277, 1295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2805, 431, 449, 1205, 1313, 1349, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2865, 449, 467, 1223, 1349, 1385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2925, 467, 485, 1241, 1385, 1421, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2985, 485, 503, 1259, 1421, 1457, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3045, 503, 521, 1277, 1457, 1493, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 3105, 521, 539, 1295, 1493, 1529, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3165, 575, 605, 1349, 1565, 1625, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3265, 605, 635, 1385, 1625, 1685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3365, 635, 665, 1421, 1685, 1745, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3465, 665, 695, 1457, 1745, 1805, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3565, 695, 725, 1493, 1805, 1865, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3665, 725, 755, 1529, 1865, 1925, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3765, 815, 860, 1625, 1985, 2075, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3915, 860, 905, 1685, 2075, 2165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4065, 905, 950, 1745, 2165, 2255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4215, 950, 995, 1805, 2255, 2345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4365, 995, 1040, 1865, 2345, 2435, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 4515, 1040, 1085, 1925, 2435, 2525, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4665, 1175, 1181, 2615, 2625, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4680, 1181, 1187, 2625, 2635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4695, 1187, 1193, 2635, 2645, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4710, 1205, 1223, 2615, 2655, 2685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4755, 1223, 1241, 2625, 2685, 2715, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4800, 1241, 1259, 2635, 2715, 2745, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4845, 1259, 1277, 2645, 2745, 2775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4890, 1313, 1349, 2655, 2805, 2865, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4980, 1349, 1385, 2685, 2865, 2925, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 5070, 1385, 1421, 2715, 2925, 2985, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 5160, 1421, 1457, 2745, 2985, 3045, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 5250, 1457, 1493, 2775, 3045, 3105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5340, 1565, 1625, 2865, 3165, 3265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5490, 1625, 1685, 2925, 3265, 3365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5640, 1685, 1745, 2985, 3365, 3465, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5790, 1745, 1805, 3045, 3465, 3565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5940, 1805, 1865, 3105, 3565, 3665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6090, 1985, 2075, 3265, 3765, 3915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6315, 2075, 2165, 3365, 3915, 4065, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6540, 2165, 2255, 3465, 4065, 4215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6765, 2255, 2345, 3565, 4215, 4365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 6990, 2345, 2435, 3665, 4365, 4515, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 7215, 2615, 2625, 4665, 4680, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 7236, 2625, 2635, 4680, 4695, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 7257, 2655, 2685, 4665, 4710, 4755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 7320, 2685, 2715, 4680, 4755, 4800, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 7383, 2715, 2745, 4695, 4800, 4845, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 7446, 2805, 2865, 4710, 4890, 4980, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 7572, 2865, 2925, 4755, 4980, 5070, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 7698, 2925, 2985, 4800, 5070, 5160, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 7824, 2985, 3045, 4845, 5160, 5250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 7950, 3165, 3265, 4980, 5340, 5490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 8160, 3265, 3365, 5070, 5490, 5640, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 8370, 3365, 3465, 5160, 5640, 5790, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 8580, 3465, 3565, 5250, 5790, 5940, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 8790, 3765, 3915, 5490, 6090, 6315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 9105, 3915, 4065, 5640, 6315, 6540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 9420, 4065, 4215, 5790, 6540, 6765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 9735, 4215, 4365, 5940, 6765, 6990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 10050, 4665, 4680, 7215, 7236, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 10078, 4710, 4755, 7215, 7257, 7320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 10162, 4755, 4800, 7236, 7320, 7383, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 10246, 4890, 4980, 7257, 7446, 7572, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 10414, 4980, 5070, 7320, 7572, 7698, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 10582, 5070, 5160, 7383, 7698, 7824, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 10750, 5340, 5490, 7572, 7950, 8160, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 11030, 5490, 5640, 7698, 8160, 8370, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 11310, 5640, 5790, 7824, 8370, 8580, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 11590, 6090, 6315, 8160, 8790, 9105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 12010, 6315, 6540, 8370, 9105, 9420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 12430, 6540, 6765, 8580, 9420, 9735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 12850, 7257, 7320, 10050, 10078, 10162, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 12958, 7446, 7572, 10078, 10246, 10414, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 13174, 7572, 7698, 10162, 10414, 10582, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 13390, 7950, 8160, 10414, 10750, 11030, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 13750, 8160, 8370, 10582, 11030, 11310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 14110, 8790, 9105, 11030, 11590, 12010, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 14650, 9105, 9420, 11310, 12010, 12430, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsd(pbuffer, 15190, 10246, 10414, 12850, 12958, 13174, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsf(pbuffer, 15460, 10750, 11030, 13174, 13390, 13750, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsg(pbuffer, 15910, 11590, 12010, 13750, 14110, 14650, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1313, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 1565, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 96, pbuffer, 1985, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 186, pbuffer, 2805, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 246, pbuffer, 3165, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 346, pbuffer, 3765, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 496, pbuffer, 4890, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 586, pbuffer, 5340, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 736, pbuffer, 6090, 225, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {4890, 4980});

                pbuffer.scale(2.0 * b_exp, {5340, 5490});

                pbuffer.scale(2.0 * b_exp, {6090, 6315});

                pbuffer.scale(2.0 * b_exp, {7446, 7572});

                pbuffer.scale(2.0 * b_exp, {7950, 8160});

                pbuffer.scale(2.0 * b_exp, {8790, 9105});

                pbuffer.scale(2.0 * b_exp, {10246, 10414});

                pbuffer.scale(2.0 * b_exp, {10750, 11030});

                pbuffer.scale(2.0 * b_exp, {11590, 12010});

                t2cfunc::reduce(cbuffer, 961, pbuffer, 4890, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1051, pbuffer, 5340, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1201, pbuffer, 6090, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1426, pbuffer, 7446, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1552, pbuffer, 7950, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1762, pbuffer, 8790, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2077, pbuffer, 10246, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2245, pbuffer, 10750, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2525, pbuffer, 11590, 420, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1313, 1349});

                pbuffer.scale(2.0 * a_exp, {1565, 1625});

                pbuffer.scale(2.0 * a_exp, {1985, 2075});

                pbuffer.scale(2.0 * a_exp, {2805, 2865});

                pbuffer.scale(2.0 * a_exp, {3165, 3265});

                pbuffer.scale(2.0 * a_exp, {3765, 3915});

                pbuffer.scale(a_exp / b_exp, {4890, 4980});

                pbuffer.scale(a_exp / b_exp, {5340, 5490});

                pbuffer.scale(a_exp / b_exp, {6090, 6315});

                pbuffer.scale(a_exp / b_exp, {7446, 7572});

                pbuffer.scale(a_exp / b_exp, {7950, 8160});

                pbuffer.scale(a_exp / b_exp, {8790, 9105});

                pbuffer.scale(a_exp / b_exp, {10246, 10414});

                pbuffer.scale(a_exp / b_exp, {10750, 11030});

                pbuffer.scale(a_exp / b_exp, {11590, 12010});

                t2cfunc::reduce(cbuffer, 2945, pbuffer, 1313, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2981, pbuffer, 1565, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3041, pbuffer, 1985, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3131, pbuffer, 2805, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3191, pbuffer, 3165, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3291, pbuffer, 3765, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3441, pbuffer, 4890, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3531, pbuffer, 5340, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3681, pbuffer, 6090, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3906, pbuffer, 7446, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4032, pbuffer, 7950, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4242, pbuffer, 8790, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4557, pbuffer, 10246, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4725, pbuffer, 10750, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5005, pbuffer, 11590, 420, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {4890, 4980});

                pbuffer.scale(2.0 * b_exp, {5340, 5490});

                pbuffer.scale(2.0 * b_exp, {6090, 6315});

                pbuffer.scale(2.0 * b_exp, {7446, 7572});

                pbuffer.scale(2.0 * b_exp, {7950, 8160});

                pbuffer.scale(2.0 * b_exp, {8790, 9105});

                pbuffer.scale(2.0 * b_exp, {10246, 10414});

                pbuffer.scale(2.0 * b_exp, {10750, 11030});

                pbuffer.scale(2.0 * b_exp, {11590, 12010});

                pbuffer.scale(4.0 * a_exp * b_exp, {12958, 13174});

                pbuffer.scale(4.0 * a_exp * b_exp, {13390, 13750});

                pbuffer.scale(4.0 * a_exp * b_exp, {14110, 14650});

                pbuffer.scale(4.0 * a_exp * b_exp, {15190, 15460});

                pbuffer.scale(4.0 * a_exp * b_exp, {15460, 15910});

                pbuffer.scale(4.0 * a_exp * b_exp, {15910, 16585});

                t2cfunc::reduce(cbuffer, 5425, pbuffer, 4890, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5515, pbuffer, 5340, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5665, pbuffer, 6090, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5890, pbuffer, 7446, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6016, pbuffer, 7950, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6226, pbuffer, 8790, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6541, pbuffer, 10246, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6709, pbuffer, 10750, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6989, pbuffer, 11590, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7409, pbuffer, 12958, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7625, pbuffer, 13390, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7985, pbuffer, 14110, 540, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8525, pbuffer, 15190, 270, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8795, pbuffer, 15460, 450, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9245, pbuffer, 15910, 675, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 36, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 108, cbuffer, 36, 96, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 288, 0, 108, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 504, cbuffer, 186, 246, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 684, cbuffer, 246, 346, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 984, 504, 684, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2424, cbuffer, 496, 586, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 2694, cbuffer, 586, 736, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 3144, 2424, 2694, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 7572, cbuffer, 961, 1051, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 7842, cbuffer, 1051, 1201, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 8292, 7572, 7842, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8832, cbuffer, 1426, 1552, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9210, cbuffer, 1552, 1762, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 9840, 8832, 9210, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 10596, cbuffer, 2077, 2245, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 11100, cbuffer, 2245, 2525, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 11940, 10596, 11100, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 12948, cbuffer, 2945, 2981, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 13056, cbuffer, 2981, 3041, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 13236, 12948, 13056, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 13452, cbuffer, 3131, 3191, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 13632, cbuffer, 3191, 3291, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 13932, 13452, 13632, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 15372, cbuffer, 3441, 3531, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 15642, cbuffer, 3531, 3681, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 16092, 15372, 15642, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 18252, cbuffer, 3906, 4032, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 18630, cbuffer, 4032, 4242, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 19260, 18252, 18630, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 22284, cbuffer, 4557, 4725, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 22788, cbuffer, 4725, 5005, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 23628, 22284, 22788, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 31548, cbuffer, 5425, 5515, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 31818, cbuffer, 5515, 5665, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 32268, 31548, 31818, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 32808, cbuffer, 5890, 6016, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 33186, cbuffer, 6016, 6226, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 33816, 32808, 33186, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 34572, cbuffer, 6541, 6709, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 35076, cbuffer, 6709, 6989, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 35916, 34572, 35076, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 36924, cbuffer, 7409, 7625, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 37572, cbuffer, 7625, 7985, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 38652, 36924, 37572, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 39948, cbuffer, 8525, 8795, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 40758, cbuffer, 8795, 9245, cfactors, 6, 0, 8);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 42108, 39948, 40758, cfactors, 6, 0, 8);

            t4cfunc::ket_transform<2, 2>(skbuffer, 0, ckbuffer, 288, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 150, ckbuffer, 984, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 3400, ckbuffer, 3144, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 129175, ckbuffer, 8292, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 129550, ckbuffer, 9840, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 130075, ckbuffer, 11940, 0, 6);

            t4cfunc::ket_transform<2, 2>(skbuffer, 130775, ckbuffer, 13236, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 130925, ckbuffer, 13932, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 131925, ckbuffer, 16092, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 133425, ckbuffer, 19260, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 135525, ckbuffer, 23628, 0, 6);

            t4cfunc::ket_transform<2, 2>(skbuffer, 144475, ckbuffer, 32268, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 144850, ckbuffer, 33816, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 145375, ckbuffer, 35916, 0, 6);

            t4cfunc::ket_transform<2, 2>(skbuffer, 146075, ckbuffer, 38652, 0, 7);

            t4cfunc::ket_transform<2, 2>(skbuffer, 146975, ckbuffer, 42108, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 20875, 150, 3400, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 141025, 130925, 131925, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 141775, 131925, 133425, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 142900, 133425, 135525, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 23875, 150, 141025, 141775, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 36250, 3400, 141775, 142900, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 68425, 20875, 23875, 36250, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 400, 0, 129175, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 3775, 150, 129550, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 8275, 3400, 130075, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 21625, 150, 400, 3775, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pgxx(skbuffer, 32875, 3400, 3775, 8275, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dfxx(skbuffer, 63925, 20875, 21625, 32875, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 131175, 130775, 144475, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 132300, 130925, 144850, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 133950, 131925, 145375, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 136225, 133425, 146075, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_skxx(skbuffer, 138325, 135525, 146975, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 1150, 130925, 131175, 132300, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 4900, 131925, 132300, 133950, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 9850, 133425, 133950, 136225, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sixx(skbuffer, 14575, 135525, 136225, 138325, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 26125, 400, 141025, 1150, 4900, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 39625, 3775, 141775, 4900, 9850, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_phxx(skbuffer, 49750, 8275, 142900, 9850, 14575, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 72925, 21625, 23875, 26125, 39625, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dgxx(skbuffer, 86425, 32875, 36250, 39625, 49750, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ffxx(skbuffer, 106675, 63925, 68425, 72925, 86425, r_ab, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 106675, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1225, skbuffer, 109175, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2450, skbuffer, 111675, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 3675, skbuffer, 114175, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 4900, skbuffer, 116675, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 6125, skbuffer, 119175, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 7350, skbuffer, 121675, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 8575, skbuffer, 124175, 2, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 9800, skbuffer, 126675, 2, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 2, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFFDD_hpp */
