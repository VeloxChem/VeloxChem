#ifndef ElectronRepulsionGeom1100RecFDPP_hpp
#define ElectronRepulsionGeom1100RecFDPP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
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
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSP.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FD|1/|r-r'||PP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_fdpp(T& distributor,
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

    CSimdArray<double> pbuffer(3186, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2070, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(5067, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(34983, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2835, 1);

    // setup Boys fuction data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_data(11, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 10, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 10);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 10, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 25, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 28, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 31, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 34, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 37, 0, 1, 10, 13, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 43, 1, 2, 13, 16, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 49, 2, 3, 16, 19, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 55, 3, 4, 19, 22, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 61, 4, 5, 22, 25, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 67, 5, 6, 25, 28, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 73, 6, 7, 28, 31, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 79, 7, 8, 31, 34, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 85, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 88, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 91, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 94, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 97, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 100, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 103, 1, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 112, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 121, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 130, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 139, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 148, 6, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 157, 7, 28, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 166, 13, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 184, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 202, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 220, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 238, 25, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 256, 28, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 274, 31, 73, 79, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 292, 1, 2, 85, 88, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 298, 2, 3, 88, 91, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 304, 3, 4, 91, 94, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 310, 4, 5, 94, 97, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 316, 5, 6, 97, 100, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 322, 10, 13, 85, 103, 112, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 340, 13, 16, 88, 112, 121, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 358, 16, 19, 91, 121, 130, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 376, 19, 22, 94, 130, 139, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 394, 22, 25, 97, 139, 148, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 412, 25, 28, 100, 148, 157, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 430, 37, 43, 112, 166, 184, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 466, 43, 49, 121, 184, 202, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 502, 49, 55, 130, 202, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 538, 55, 61, 139, 220, 238, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 574, 61, 67, 148, 238, 256, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 610, 67, 73, 157, 256, 274, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 646, 85, 88, 292, 298, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 656, 88, 91, 298, 304, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 666, 91, 94, 304, 310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 676, 94, 97, 310, 316, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 686, 103, 112, 292, 322, 340, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 716, 112, 121, 298, 340, 358, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 746, 121, 130, 304, 358, 376, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 776, 130, 139, 310, 376, 394, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 806, 139, 148, 316, 394, 412, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 836, 166, 184, 340, 430, 466, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 896, 184, 202, 358, 466, 502, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 956, 202, 220, 376, 502, 538, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1016, 220, 238, 394, 538, 574, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1076, 238, 256, 412, 574, 610, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1136, 292, 298, 646, 656, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1151, 298, 304, 656, 666, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1166, 304, 310, 666, 676, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1181, 322, 340, 646, 686, 716, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1226, 340, 358, 656, 716, 746, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1271, 358, 376, 666, 746, 776, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1316, 376, 394, 676, 776, 806, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1361, 430, 466, 716, 836, 896, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1451, 466, 502, 746, 896, 956, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1541, 502, 538, 776, 956, 1016, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1631, 538, 574, 806, 1016, 1076, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 1721, 646, 656, 1136, 1151, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 1742, 656, 666, 1151, 1166, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1763, 686, 716, 1136, 1181, 1226, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1826, 716, 746, 1151, 1226, 1271, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1889, 746, 776, 1166, 1271, 1316, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 1952, 836, 896, 1226, 1361, 1451, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2078, 896, 956, 1271, 1451, 1541, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2204, 956, 1016, 1316, 1541, 1631, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 2330, 1136, 1151, 1721, 1742, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 2358, 1181, 1226, 1721, 1763, 1826, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 2442, 1226, 1271, 1742, 1826, 1889, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 2526, 1361, 1451, 1826, 1952, 2078, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 2694, 1451, 1541, 1889, 2078, 2204, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 2862, 1763, 1826, 2330, 2358, 2442, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 2970, 1952, 2078, 2442, 2526, 2694, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 103, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 166, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 27, pbuffer, 322, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 45, pbuffer, 430, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 81, pbuffer, 686, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 111, pbuffer, 836, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {686, 716});

                pbuffer.scale(2.0 * b_exp, {836, 896});

                pbuffer.scale(2.0 * b_exp, {1181, 1226});

                pbuffer.scale(2.0 * b_exp, {1361, 1451});

                pbuffer.scale(2.0 * b_exp, {1763, 1826});

                pbuffer.scale(2.0 * b_exp, {1952, 2078});

                t2cfunc::reduce(cbuffer, 171, pbuffer, 686, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 201, pbuffer, 836, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 261, pbuffer, 1181, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 306, pbuffer, 1361, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 396, pbuffer, 1763, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 459, pbuffer, 1952, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {103, 112});

                pbuffer.scale(2.0 * a_exp, {166, 184});

                pbuffer.scale(2.0 * a_exp, {322, 340});

                pbuffer.scale(2.0 * a_exp, {430, 466});

                pbuffer.scale(a_exp / b_exp, {686, 716});

                pbuffer.scale(a_exp / b_exp, {836, 896});

                pbuffer.scale(a_exp / b_exp, {1181, 1226});

                pbuffer.scale(a_exp / b_exp, {1361, 1451});

                pbuffer.scale(a_exp / b_exp, {1763, 1826});

                pbuffer.scale(a_exp / b_exp, {1952, 2078});

                t2cfunc::reduce(cbuffer, 585, pbuffer, 103, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 594, pbuffer, 166, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 612, pbuffer, 322, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 630, pbuffer, 430, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 666, pbuffer, 686, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 696, pbuffer, 836, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 756, pbuffer, 1181, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 801, pbuffer, 1361, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 891, pbuffer, 1763, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 954, pbuffer, 1952, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {686, 716});

                pbuffer.scale(2.0 * b_exp, {836, 896});

                pbuffer.scale(2.0 * b_exp, {1181, 1226});

                pbuffer.scale(2.0 * b_exp, {1361, 1451});

                pbuffer.scale(2.0 * b_exp, {1763, 1826});

                pbuffer.scale(2.0 * b_exp, {1952, 2078});

                pbuffer.scale(4.0 * a_exp * b_exp, {2358, 2442});

                pbuffer.scale(4.0 * a_exp * b_exp, {2526, 2694});

                pbuffer.scale(4.0 * a_exp * b_exp, {2862, 2970});

                pbuffer.scale(4.0 * a_exp * b_exp, {2970, 3186});

                t2cfunc::reduce(cbuffer, 1080, pbuffer, 686, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1110, pbuffer, 836, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1170, pbuffer, 1181, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1215, pbuffer, 1361, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1305, pbuffer, 1763, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1368, pbuffer, 1952, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1494, pbuffer, 2358, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1578, pbuffer, 2526, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1746, pbuffer, 2862, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1854, pbuffer, 2970, 216, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 0, cbuffer, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 27, cbuffer, 27, 45, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 243, cbuffer, 81, 111, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1008, cbuffer, 171, 201, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1098, cbuffer, 261, 306, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1233, cbuffer, 396, 459, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1422, cbuffer, 585, 594, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1449, cbuffer, 612, 630, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1665, cbuffer, 666, 696, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2025, cbuffer, 756, 801, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2565, cbuffer, 891, 954, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4077, cbuffer, 1080, 1110, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4167, cbuffer, 1170, 1215, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4302, cbuffer, 1305, 1368, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4491, cbuffer, 1494, 1578, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4743, cbuffer, 1746, 1854, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<1, 1>(skbuffer, 0, ckbuffer, 0, 0, 1);

            t4cfunc::ket_transform<1, 1>(skbuffer, 27, ckbuffer, 27, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 729, ckbuffer, 243, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 30087, ckbuffer, 1008, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 30177, ckbuffer, 1098, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 30312, ckbuffer, 1233, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 30501, ckbuffer, 1422, 0, 1);

            t4cfunc::ket_transform<1, 1>(skbuffer, 30528, ckbuffer, 1449, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 30744, ckbuffer, 1665, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 31104, ckbuffer, 2025, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 31644, ckbuffer, 2565, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 33993, ckbuffer, 4077, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 34083, ckbuffer, 4167, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 34218, ckbuffer, 4302, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 34407, ckbuffer, 4491, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 34659, ckbuffer, 4743, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 5220, 27, 729, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 33156, 30528, 30744, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 33318, 30744, 31104, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 33588, 31104, 31644, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 5868, 27, 33156, 33318, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 8622, 729, 33318, 33588, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ddxx(skbuffer, 16479, 5220, 5868, 8622, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 81, 0, 30087, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 819, 27, 30177, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 1899, 729, 30312, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 5382, 27, 81, 819, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 7812, 729, 819, 1899, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ddxx(skbuffer, 15507, 5220, 5382, 7812, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 30582, 30501, 33993, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 30834, 30528, 34083, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 31239, 30744, 34218, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 31833, 31104, 34407, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 32400, 31644, 34659, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 243, 30528, 30582, 30834, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 1089, 30744, 30834, 31239, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 2304, 31104, 31239, 31833, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 3519, 31644, 31833, 32400, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 6354, 81, 33156, 243, 1089, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 9432, 819, 33318, 1089, 2304, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 11862, 1899, 33588, 2304, 3519, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 17451, 5382, 5868, 6354, 9432, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 20367, 7812, 8622, 9432, 11862, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_fdxx(skbuffer, 25227, 15507, 16479, 17451, 20367, r_ab, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 25227, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 315, skbuffer, 25767, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 630, skbuffer, 26307, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 945, skbuffer, 26847, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1260, skbuffer, 27387, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1575, skbuffer, 27927, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1890, skbuffer, 28467, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2205, skbuffer, 29007, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2520, skbuffer, 29547, 1, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 1, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFDPP_hpp */
