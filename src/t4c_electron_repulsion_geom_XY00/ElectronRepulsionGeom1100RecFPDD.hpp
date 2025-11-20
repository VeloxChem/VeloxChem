#ifndef ElectronRepulsionGeom1100RecFPDD_hpp
#define ElectronRepulsionGeom1100RecFPDD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecXXDD.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0100ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecFPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSPXX.hpp"
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
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FP|1/|r-r'||DD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_fpdd(T& distributor,
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

    CSimdArray<double> pbuffer(6998, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(4836, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(21096, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(57225, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(4725, 1);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 11, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 11);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 29, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 32, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 35, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 38, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 41, 0, 1, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 47, 1, 2, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 2, 3, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 3, 4, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 65, 4, 5, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 71, 5, 6, 26, 29, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 77, 6, 7, 29, 32, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 83, 7, 8, 32, 35, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 89, 8, 9, 35, 38, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 11, 14, 41, 47, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 14, 17, 47, 53, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 17, 20, 53, 59, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 20, 23, 59, 65, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 23, 26, 65, 71, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 26, 29, 71, 77, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 29, 32, 77, 83, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 32, 35, 83, 89, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 175, 41, 47, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 190, 47, 53, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 205, 53, 59, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 220, 59, 65, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 235, 65, 71, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 250, 71, 77, 145, 155, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 265, 77, 83, 155, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 280, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 283, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 286, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 289, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 292, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 301, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 310, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 319, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 328, 6, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 337, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 355, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 373, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 391, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 409, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 427, 29, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 445, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 475, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 505, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 535, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 565, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 595, 77, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 625, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 670, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 715, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 760, 135, 220, 235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 805, 145, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 850, 155, 250, 265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 895, 2, 3, 280, 283, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 901, 3, 4, 283, 286, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 907, 4, 5, 286, 289, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 913, 14, 17, 280, 292, 301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 931, 17, 20, 283, 301, 310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 949, 20, 23, 286, 310, 319, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 967, 23, 26, 289, 319, 328, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 985, 41, 47, 292, 337, 355, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1021, 47, 53, 301, 355, 373, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1057, 53, 59, 310, 373, 391, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1093, 59, 65, 319, 391, 409, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1129, 65, 71, 328, 409, 427, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1165, 95, 105, 355, 445, 475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1225, 105, 115, 373, 475, 505, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1285, 115, 125, 391, 505, 535, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1345, 125, 135, 409, 535, 565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1405, 135, 145, 427, 565, 595, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1465, 175, 190, 475, 625, 670, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1555, 190, 205, 505, 670, 715, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1645, 205, 220, 535, 715, 760, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1735, 220, 235, 565, 760, 805, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1825, 235, 250, 595, 805, 850, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1915, 280, 283, 895, 901, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1925, 283, 286, 901, 907, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1935, 292, 301, 895, 913, 931, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1965, 301, 310, 901, 931, 949, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1995, 310, 319, 907, 949, 967, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2025, 337, 355, 913, 985, 1021, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2085, 355, 373, 931, 1021, 1057, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2145, 373, 391, 949, 1057, 1093, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2205, 391, 409, 967, 1093, 1129, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2265, 445, 475, 1021, 1165, 1225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2365, 475, 505, 1057, 1225, 1285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2465, 505, 535, 1093, 1285, 1345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2565, 535, 565, 1129, 1345, 1405, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2665, 625, 670, 1225, 1465, 1555, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2815, 670, 715, 1285, 1555, 1645, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2965, 715, 760, 1345, 1645, 1735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3115, 760, 805, 1405, 1735, 1825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 3265, 895, 901, 1915, 1925, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3280, 913, 931, 1915, 1935, 1965, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3325, 931, 949, 1925, 1965, 1995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3370, 985, 1021, 1935, 2025, 2085, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3460, 1021, 1057, 1965, 2085, 2145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3550, 1057, 1093, 1995, 2145, 2205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3640, 1165, 1225, 2085, 2265, 2365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3790, 1225, 1285, 2145, 2365, 2465, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3940, 1285, 1345, 2205, 2465, 2565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4090, 1465, 1555, 2365, 2665, 2815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4315, 1555, 1645, 2465, 2815, 2965, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4540, 1645, 1735, 2565, 2965, 3115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 4765, 1935, 1965, 3265, 3280, 3325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 4828, 2025, 2085, 3280, 3370, 3460, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 4954, 2085, 2145, 3325, 3460, 3550, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5080, 2265, 2365, 3460, 3640, 3790, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5290, 2365, 2465, 3550, 3790, 3940, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 5500, 2665, 2815, 3790, 4090, 4315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 5815, 2815, 2965, 3940, 4315, 4540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 6130, 3370, 3460, 4765, 4828, 4954, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 6298, 3640, 3790, 4954, 5080, 5290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 6578, 4090, 4315, 5290, 5500, 5815, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 41, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 31, pbuffer, 337, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 49, pbuffer, 445, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 79, pbuffer, 625, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 124, pbuffer, 985, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 1165, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 220, pbuffer, 1465, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {985, 1021});

                pbuffer.scale(2.0 * b_exp, {1165, 1225});

                pbuffer.scale(2.0 * b_exp, {1465, 1555});

                pbuffer.scale(2.0 * b_exp, {2025, 2085});

                pbuffer.scale(2.0 * b_exp, {2265, 2365});

                pbuffer.scale(2.0 * b_exp, {2665, 2815});

                pbuffer.scale(2.0 * b_exp, {3370, 3460});

                pbuffer.scale(2.0 * b_exp, {3640, 3790});

                pbuffer.scale(2.0 * b_exp, {4090, 4315});

                t2cfunc::reduce(cbuffer, 310, pbuffer, 985, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 346, pbuffer, 1165, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 406, pbuffer, 1465, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 496, pbuffer, 2025, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 556, pbuffer, 2265, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 656, pbuffer, 2665, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 806, pbuffer, 3370, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 896, pbuffer, 3640, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1046, pbuffer, 4090, 225, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {41, 47});

                pbuffer.scale(2.0 * a_exp, {95, 105});

                pbuffer.scale(2.0 * a_exp, {175, 190});

                pbuffer.scale(2.0 * a_exp, {337, 355});

                pbuffer.scale(2.0 * a_exp, {445, 475});

                pbuffer.scale(2.0 * a_exp, {625, 670});

                pbuffer.scale(a_exp / b_exp, {985, 1021});

                pbuffer.scale(a_exp / b_exp, {1165, 1225});

                pbuffer.scale(a_exp / b_exp, {1465, 1555});

                pbuffer.scale(a_exp / b_exp, {2025, 2085});

                pbuffer.scale(a_exp / b_exp, {2265, 2365});

                pbuffer.scale(a_exp / b_exp, {2665, 2815});

                pbuffer.scale(a_exp / b_exp, {3370, 3460});

                pbuffer.scale(a_exp / b_exp, {3640, 3790});

                pbuffer.scale(a_exp / b_exp, {4090, 4315});

                t2cfunc::reduce(cbuffer, 1271, pbuffer, 41, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1277, pbuffer, 95, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1287, pbuffer, 175, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1302, pbuffer, 337, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1320, pbuffer, 445, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1350, pbuffer, 625, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1395, pbuffer, 985, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1431, pbuffer, 1165, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1491, pbuffer, 1465, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1581, pbuffer, 2025, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1641, pbuffer, 2265, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1741, pbuffer, 2665, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1891, pbuffer, 3370, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1981, pbuffer, 3640, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2131, pbuffer, 4090, 225, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {985, 1021});

                pbuffer.scale(2.0 * b_exp, {1165, 1225});

                pbuffer.scale(2.0 * b_exp, {1465, 1555});

                pbuffer.scale(2.0 * b_exp, {2025, 2085});

                pbuffer.scale(2.0 * b_exp, {2265, 2365});

                pbuffer.scale(2.0 * b_exp, {2665, 2815});

                pbuffer.scale(2.0 * b_exp, {3370, 3460});

                pbuffer.scale(2.0 * b_exp, {3640, 3790});

                pbuffer.scale(2.0 * b_exp, {4090, 4315});

                pbuffer.scale(4.0 * a_exp * b_exp, {4828, 4954});

                pbuffer.scale(4.0 * a_exp * b_exp, {5080, 5290});

                pbuffer.scale(4.0 * a_exp * b_exp, {5500, 5815});

                pbuffer.scale(4.0 * a_exp * b_exp, {6130, 6298});

                pbuffer.scale(4.0 * a_exp * b_exp, {6298, 6578});

                pbuffer.scale(4.0 * a_exp * b_exp, {6578, 6998});

                t2cfunc::reduce(cbuffer, 2356, pbuffer, 985, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2392, pbuffer, 1165, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2452, pbuffer, 1465, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2542, pbuffer, 2025, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2602, pbuffer, 2265, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2702, pbuffer, 2665, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2852, pbuffer, 3370, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2942, pbuffer, 3640, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3092, pbuffer, 4090, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3317, pbuffer, 4828, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3443, pbuffer, 5080, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3653, pbuffer, 5500, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3968, pbuffer, 6130, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4136, pbuffer, 6298, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4416, pbuffer, 6578, 420, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 6, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 18, cbuffer, 6, 16, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 48, 0, 18, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 84, cbuffer, 31, 49, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 138, cbuffer, 49, 79, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 228, 84, 138, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 660, cbuffer, 124, 160, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 768, cbuffer, 160, 220, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 948, 660, 768, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2892, cbuffer, 310, 346, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 3000, cbuffer, 346, 406, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 3180, 2892, 3000, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3396, cbuffer, 496, 556, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 3576, cbuffer, 556, 656, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 3876, 3396, 3576, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4236, cbuffer, 806, 896, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4506, cbuffer, 896, 1046, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 4956, 4236, 4506, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5496, cbuffer, 1271, 1277, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5514, cbuffer, 1277, 1287, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 5544, 5496, 5514, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5580, cbuffer, 1302, 1320, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5634, cbuffer, 1320, 1350, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 5724, 5580, 5634, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 6156, cbuffer, 1395, 1431, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 6264, cbuffer, 1431, 1491, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 6444, 6156, 6264, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 7308, cbuffer, 1581, 1641, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 7488, cbuffer, 1641, 1741, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 7788, 7308, 7488, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 9228, cbuffer, 1891, 1981, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9498, cbuffer, 1981, 2131, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 9948, 9228, 9498, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 14376, cbuffer, 2356, 2392, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 14484, cbuffer, 2392, 2452, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 14664, 14376, 14484, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 14880, cbuffer, 2542, 2602, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 15060, cbuffer, 2602, 2702, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 15360, 14880, 15060, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 15720, cbuffer, 2852, 2942, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 15990, cbuffer, 2942, 3092, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 16440, 15720, 15990, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 16980, cbuffer, 3317, 3443, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 17358, cbuffer, 3443, 3653, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 17988, 16980, 17358, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 18744, cbuffer, 3968, 4136, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 19248, cbuffer, 4136, 4416, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 20088, 18744, 19248, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<2, 2>(skbuffer, 0, ckbuffer, 48, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 25, ckbuffer, 228, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 1000, ckbuffer, 948, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 48025, ckbuffer, 3180, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 48175, ckbuffer, 3876, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 48425, ckbuffer, 4956, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 48800, ckbuffer, 5544, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 48825, ckbuffer, 5724, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 49125, ckbuffer, 6444, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 49725, ckbuffer, 7788, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 50725, ckbuffer, 9948, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 55225, ckbuffer, 14664, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 55375, ckbuffer, 15360, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 55625, ckbuffer, 16440, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 56000, ckbuffer, 17988, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 56525, ckbuffer, 20088, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 9325, 25, 1000, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 53800, 48825, 49125, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 54025, 49125, 49725, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 54475, 49725, 50725, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 10225, 25, 53800, 54025, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 14275, 1000, 54025, 54475, r_ab, 2, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dpxx(skbuffer, 27775, 9325, 10225, 14275, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 100, 0, 48025, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 1150, 25, 48175, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 2950, 1000, 48425, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 9550, 25, 100, 1150, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 12925, 1000, 1150, 2950, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dpxx(skbuffer, 26425, 9325, 9550, 12925, r_ab, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 48900, 48800, 55225, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 49275, 48825, 55375, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 49975, 49125, 55625, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 51100, 49725, 56000, 2, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 52225, 50725, 56525, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 325, 48825, 48900, 49275, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 1600, 49125, 49275, 49975, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 3700, 49725, 49975, 51100, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 5950, 50725, 51100, 52225, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 10900, 100, 53800, 325, 1600, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 15625, 1150, 54025, 1600, 3700, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 19675, 2950, 54475, 3700, 5950, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 29125, 9550, 10225, 10900, 15625, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 33175, 12925, 14275, 15625, 19675, r_ab, 2, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_fpxx(skbuffer, 41275, 26425, 27775, 29125, 33175, r_ab, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 41275, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 525, skbuffer, 42025, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1050, skbuffer, 42775, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1575, skbuffer, 43525, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2100, skbuffer, 44275, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2625, skbuffer, 45025, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 3150, skbuffer, 45775, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 3675, skbuffer, 46525, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 4200, skbuffer, 47275, 2, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 2, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFPDD_hpp */
