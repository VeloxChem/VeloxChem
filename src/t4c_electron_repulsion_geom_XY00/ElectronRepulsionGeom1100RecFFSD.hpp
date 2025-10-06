#ifndef ElectronRepulsionGeom1100RecFFSD_hpp
#define ElectronRepulsionGeom1100RecFFSD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
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
#include "ElectronRepulsionPrimRecSLSD.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FF|1/|r-r'||SD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_ffsd(T& distributor,
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

    CSimdArray<double> pbuffer(4180, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1920, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(29620, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2205, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 95, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 98, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 101, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 104, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 107, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 110, 7, 8, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 113, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 122, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 131, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 140, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 149, 6, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 158, 7, 29, 32, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 167, 8, 32, 35, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 176, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 194, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 212, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 230, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 248, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 266, 29, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 284, 32, 77, 83, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 302, 35, 83, 89, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 320, 2, 3, 95, 98, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 326, 3, 4, 98, 101, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 332, 4, 5, 101, 104, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 338, 5, 6, 104, 107, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 344, 6, 7, 107, 110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 350, 14, 17, 95, 113, 122, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 368, 17, 20, 98, 122, 131, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 386, 20, 23, 101, 131, 140, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 404, 23, 26, 104, 140, 149, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 422, 26, 29, 107, 149, 158, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 440, 29, 32, 110, 158, 167, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 458, 41, 47, 113, 176, 194, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 494, 47, 53, 122, 194, 212, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 530, 53, 59, 131, 212, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 566, 59, 65, 140, 230, 248, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 602, 65, 71, 149, 248, 266, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 638, 71, 77, 158, 266, 284, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 674, 77, 83, 167, 284, 302, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 710, 95, 98, 320, 326, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 720, 98, 101, 326, 332, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 730, 101, 104, 332, 338, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 740, 104, 107, 338, 344, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 750, 113, 122, 320, 350, 368, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 780, 122, 131, 326, 368, 386, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 810, 131, 140, 332, 386, 404, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 840, 140, 149, 338, 404, 422, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 870, 149, 158, 344, 422, 440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 900, 176, 194, 350, 458, 494, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 960, 194, 212, 368, 494, 530, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1020, 212, 230, 386, 530, 566, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1080, 230, 248, 404, 566, 602, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1140, 248, 266, 422, 602, 638, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1200, 266, 284, 440, 638, 674, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1260, 320, 326, 710, 720, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1275, 326, 332, 720, 730, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1290, 332, 338, 730, 740, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1305, 350, 368, 710, 750, 780, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1350, 368, 386, 720, 780, 810, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1395, 386, 404, 730, 810, 840, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1440, 404, 422, 740, 840, 870, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1485, 458, 494, 750, 900, 960, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1575, 494, 530, 780, 960, 1020, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1665, 530, 566, 810, 1020, 1080, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1755, 566, 602, 840, 1080, 1140, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1845, 602, 638, 870, 1140, 1200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 1935, 710, 720, 1260, 1275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 1956, 720, 730, 1275, 1290, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1977, 750, 780, 1260, 1305, 1350, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2040, 780, 810, 1275, 1350, 1395, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2103, 810, 840, 1290, 1395, 1440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2166, 900, 960, 1305, 1485, 1575, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2292, 960, 1020, 1350, 1575, 1665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2418, 1020, 1080, 1395, 1665, 1755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2544, 1080, 1140, 1440, 1755, 1845, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 2670, 1260, 1275, 1935, 1956, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 2698, 1305, 1350, 1935, 1977, 2040, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 2782, 1350, 1395, 1956, 2040, 2103, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 2866, 1485, 1575, 1977, 2166, 2292, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 3034, 1575, 1665, 2040, 2292, 2418, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 3202, 1665, 1755, 2103, 2418, 2544, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 3370, 1977, 2040, 2670, 2698, 2782, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 3478, 2166, 2292, 2698, 2866, 3034, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 3694, 2292, 2418, 2782, 3034, 3202, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsd(pbuffer, 3910, 2866, 3034, 3370, 3478, 3694, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 458, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 900, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 96, pbuffer, 1485, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {1485, 1575});

                pbuffer.scale(2.0 * b_exp, {2166, 2292});

                pbuffer.scale(2.0 * b_exp, {2866, 3034});

                t2cfunc::reduce(cbuffer, 186, pbuffer, 1485, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 276, pbuffer, 2166, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 402, pbuffer, 2866, 168, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {458, 494});

                pbuffer.scale(2.0 * a_exp, {900, 960});

                pbuffer.scale(a_exp / b_exp, {1485, 1575});

                pbuffer.scale(a_exp / b_exp, {2166, 2292});

                pbuffer.scale(a_exp / b_exp, {2866, 3034});

                t2cfunc::reduce(cbuffer, 570, pbuffer, 458, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 606, pbuffer, 900, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 666, pbuffer, 1485, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 756, pbuffer, 2166, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 882, pbuffer, 2866, 168, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {1485, 1575});

                pbuffer.scale(2.0 * b_exp, {2166, 2292});

                pbuffer.scale(2.0 * b_exp, {2866, 3034});

                pbuffer.scale(4.0 * a_exp * b_exp, {3478, 3694});

                pbuffer.scale(4.0 * a_exp * b_exp, {3910, 4180});

                t2cfunc::reduce(cbuffer, 1050, pbuffer, 1485, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1140, pbuffer, 2166, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1266, pbuffer, 2866, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1434, pbuffer, 3478, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1650, pbuffer, 3910, 270, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 2>(skbuffer, 0, cbuffer, 0, 0, 2);

            t4cfunc::ket_transform<0, 2>(skbuffer, 30, cbuffer, 36, 0, 3);

            t4cfunc::ket_transform<0, 2>(skbuffer, 680, cbuffer, 96, 0, 4);

            t4cfunc::ket_transform<0, 2>(skbuffer, 25835, cbuffer, 186, 0, 4);

            t4cfunc::ket_transform<0, 2>(skbuffer, 25910, cbuffer, 276, 0, 5);

            t4cfunc::ket_transform<0, 2>(skbuffer, 26015, cbuffer, 402, 0, 6);

            t4cfunc::ket_transform<0, 2>(skbuffer, 26155, cbuffer, 570, 0, 2);

            t4cfunc::ket_transform<0, 2>(skbuffer, 26185, cbuffer, 606, 0, 3);

            t4cfunc::ket_transform<0, 2>(skbuffer, 26385, cbuffer, 666, 0, 4);

            t4cfunc::ket_transform<0, 2>(skbuffer, 26685, cbuffer, 756, 0, 5);

            t4cfunc::ket_transform<0, 2>(skbuffer, 27105, cbuffer, 882, 0, 6);

            t4cfunc::ket_transform<0, 2>(skbuffer, 28895, cbuffer, 1050, 0, 4);

            t4cfunc::ket_transform<0, 2>(skbuffer, 28970, cbuffer, 1140, 0, 5);

            t4cfunc::ket_transform<0, 2>(skbuffer, 29075, cbuffer, 1266, 0, 6);

            t4cfunc::ket_transform<0, 2>(skbuffer, 29215, cbuffer, 1434, 0, 7);

            t4cfunc::ket_transform<0, 2>(skbuffer, 29395, cbuffer, 1650, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 4175, 30, 680, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 28205, 26185, 26385, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 28355, 26385, 26685, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 28580, 26685, 27105, r_ab, 0, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 4775, 30, 28205, 28355, r_ab, 0, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 7250, 680, 28355, 28580, r_ab, 0, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 13685, 4175, 4775, 7250, r_ab, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 80, 0, 25835, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 755, 30, 25910, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 1655, 680, 26015, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 4325, 30, 80, 755, r_ab, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pgxx(skbuffer, 6575, 680, 755, 1655, r_ab, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dfxx(skbuffer, 12785, 4175, 4325, 6575, r_ab, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 26235, 26155, 28895, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 26460, 26185, 28970, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 26790, 26385, 29075, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 27245, 26685, 29215, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_skxx(skbuffer, 27665, 27105, 29395, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 230, 26185, 26235, 26460, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 980, 26385, 26460, 26790, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 1970, 26685, 26790, 27245, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sixx(skbuffer, 2915, 27105, 27245, 27665, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 5225, 80, 28205, 230, 980, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 7925, 755, 28355, 980, 1970, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_phxx(skbuffer, 9950, 1655, 28580, 1970, 2915, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 14585, 4325, 4775, 5225, 7925, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dgxx(skbuffer, 17285, 6575, 7250, 7925, 9950, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ffxx(skbuffer, 21335, 12785, 13685, 14585, 17285, r_ab, 0, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 21335, 0, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 245, skbuffer, 21835, 0, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 490, skbuffer, 22335, 0, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 735, skbuffer, 22835, 0, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 980, skbuffer, 23335, 0, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1225, skbuffer, 23835, 0, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1470, skbuffer, 24335, 0, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1715, skbuffer, 24835, 0, 2);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1960, skbuffer, 25335, 0, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 0, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFFSD_hpp */
