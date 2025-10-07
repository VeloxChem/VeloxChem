#ifndef ElectronRepulsionGeom1100RecFPSF_hpp
#define ElectronRepulsionGeom1100RecFPSF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
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
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FP|1/|r-r'||SF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_fpsf(T& distributor,
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

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(3131, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1560, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(16023, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1323, 1);

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

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // loop over basis function pairs on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            // zero integral buffers

            cbuffer.zero();

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

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 10, 13, 37, 43, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 13, 16, 43, 49, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 16, 19, 49, 55, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 19, 22, 55, 61, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 22, 25, 61, 67, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 25, 28, 67, 73, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 28, 31, 73, 79, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 155, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 158, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 161, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 164, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 173, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 182, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 191, 6, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 200, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 218, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 236, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 254, 25, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 272, 28, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 290, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 320, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 350, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 380, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 410, 67, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 440, 73, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 470, 3, 4, 155, 158, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 476, 4, 5, 158, 161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 482, 16, 19, 155, 164, 173, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 500, 19, 22, 158, 173, 182, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 518, 22, 25, 161, 182, 191, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 536, 43, 49, 164, 200, 218, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 572, 49, 55, 173, 218, 236, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 608, 55, 61, 182, 236, 254, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 644, 61, 67, 191, 254, 272, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 680, 85, 95, 200, 290, 320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 740, 95, 105, 218, 320, 350, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 800, 105, 115, 236, 350, 380, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 860, 115, 125, 254, 380, 410, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 920, 125, 135, 272, 410, 440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 980, 155, 158, 470, 476, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 990, 164, 173, 470, 482, 500, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1020, 173, 182, 476, 500, 518, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1050, 200, 218, 482, 536, 572, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1110, 218, 236, 500, 572, 608, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1170, 236, 254, 518, 608, 644, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1230, 290, 320, 536, 680, 740, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1330, 320, 350, 572, 740, 800, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1430, 350, 380, 608, 800, 860, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1530, 380, 410, 644, 860, 920, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1630, 482, 500, 980, 990, 1020, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1675, 536, 572, 990, 1050, 1110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1765, 572, 608, 1020, 1110, 1170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 1855, 680, 740, 1050, 1230, 1330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2005, 740, 800, 1110, 1330, 1430, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2155, 800, 860, 1170, 1430, 1530, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2305, 1050, 1110, 1630, 1675, 1765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 2431, 1230, 1330, 1675, 1855, 2005, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 2641, 1330, 1430, 1765, 2005, 2155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 2851, 1855, 2005, 2305, 2431, 2641, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 290, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 40, pbuffer, 680, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {680, 740});

                pbuffer.scale(2.0 * b_exp, {1230, 1330});

                pbuffer.scale(2.0 * b_exp, {1855, 2005});

                t2cfunc::reduce(cbuffer, 100, pbuffer, 680, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 1230, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 260, pbuffer, 1855, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {85, 95});

                pbuffer.scale(2.0 * a_exp, {290, 320});

                pbuffer.scale(a_exp / b_exp, {680, 740});

                pbuffer.scale(a_exp / b_exp, {1230, 1330});

                pbuffer.scale(a_exp / b_exp, {1855, 2005});

                t2cfunc::reduce(cbuffer, 410, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 420, pbuffer, 290, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 450, pbuffer, 680, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 510, pbuffer, 1230, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 610, pbuffer, 1855, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {680, 740});

                pbuffer.scale(2.0 * b_exp, {1230, 1330});

                pbuffer.scale(2.0 * b_exp, {1855, 2005});

                pbuffer.scale(4.0 * a_exp * b_exp, {2431, 2641});

                pbuffer.scale(4.0 * a_exp * b_exp, {2851, 3131});

                t2cfunc::reduce(cbuffer, 760, pbuffer, 680, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 820, pbuffer, 1230, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 920, pbuffer, 1855, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1070, pbuffer, 2431, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1280, pbuffer, 2851, 280, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 3>(skbuffer, 0, cbuffer, 0, 0, 0);

            t4cfunc::ket_transform<0, 3>(skbuffer, 7, cbuffer, 10, 0, 1);

            t4cfunc::ket_transform<0, 3>(skbuffer, 280, cbuffer, 40, 0, 2);

            t4cfunc::ket_transform<0, 3>(skbuffer, 13447, cbuffer, 100, 0, 2);

            t4cfunc::ket_transform<0, 3>(skbuffer, 13489, cbuffer, 160, 0, 3);

            t4cfunc::ket_transform<0, 3>(skbuffer, 13559, cbuffer, 260, 0, 4);

            t4cfunc::ket_transform<0, 3>(skbuffer, 13664, cbuffer, 410, 0, 0);

            t4cfunc::ket_transform<0, 3>(skbuffer, 13671, cbuffer, 420, 0, 1);

            t4cfunc::ket_transform<0, 3>(skbuffer, 13755, cbuffer, 450, 0, 2);

            t4cfunc::ket_transform<0, 3>(skbuffer, 13923, cbuffer, 510, 0, 3);

            t4cfunc::ket_transform<0, 3>(skbuffer, 14203, cbuffer, 610, 0, 4);

            t4cfunc::ket_transform<0, 3>(skbuffer, 15463, cbuffer, 760, 0, 2);

            t4cfunc::ket_transform<0, 3>(skbuffer, 15505, cbuffer, 820, 0, 3);

            t4cfunc::ket_transform<0, 3>(skbuffer, 15575, cbuffer, 920, 0, 4);

            t4cfunc::ket_transform<0, 3>(skbuffer, 15680, cbuffer, 1070, 0, 5);

            t4cfunc::ket_transform<0, 3>(skbuffer, 15827, cbuffer, 1280, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 2611, 7, 280, r_ab, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 15064, 13671, 13755, r_ab, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 15127, 13755, 13923, r_ab, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 15253, 13923, 14203, r_ab, 0, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 2863, 7, 15064, 15127, r_ab, 0, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 3997, 280, 15127, 15253, r_ab, 0, 3);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dpxx(skbuffer, 7777, 2611, 2863, 3997, r_ab, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 28, 0, 13447, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 322, 7, 13489, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 826, 280, 13559, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 2674, 7, 28, 322, r_ab, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 3619, 280, 322, 826, r_ab, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dpxx(skbuffer, 7399, 2611, 2674, 3619, r_ab, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 13692, 13664, 15463, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 13797, 13671, 15505, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 13993, 13755, 15575, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 14308, 13923, 15680, 0, 3);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 14623, 14203, 15827, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 91, 13671, 13692, 13797, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 448, 13755, 13797, 13993, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 1036, 13923, 13993, 14308, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 1666, 14203, 14308, 14623, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 3052, 28, 15064, 91, 448, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 4375, 322, 15127, 448, 1036, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 5509, 826, 15253, 1036, 1666, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 8155, 2674, 2863, 3052, 4375, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 9289, 3619, 3997, 4375, 5509, r_ab, 0, 3);

            erirec::comp_bra_geom11_hrr_electron_repulsion_fpxx(skbuffer, 11557, 7399, 7777, 8155, 9289, r_ab, 0, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 11557, 0, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 147, skbuffer, 11767, 0, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 294, skbuffer, 11977, 0, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 441, skbuffer, 12187, 0, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 588, skbuffer, 12397, 0, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 735, skbuffer, 12607, 0, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 882, skbuffer, 12817, 0, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1029, skbuffer, 13027, 0, 3);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1176, skbuffer, 13237, 0, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 0, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFPSF_hpp */
