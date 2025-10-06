#ifndef ElectronRepulsionGeom1100RecFFPP_hpp
#define ElectronRepulsionGeom1100RecFFPP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
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
#include "ElectronRepulsionPrimRecSKSS.hpp"
#include "ElectronRepulsionPrimRecSLSD.hpp"
#include "ElectronRepulsionPrimRecSLSP.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FF|1/|r-r'||PP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_ffpp(T& distributor,
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

    CSimdArray<double> pbuffer(4791, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2880, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(7092, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(53316, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3969, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 95, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 98, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 101, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 104, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 107, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 110, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 113, 7, 8, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 116, 1, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 125, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 134, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 143, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 152, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 161, 6, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 170, 7, 29, 32, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 179, 8, 32, 35, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 188, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 206, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 224, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 242, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 260, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 278, 29, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 296, 32, 77, 83, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 314, 35, 83, 89, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 332, 1, 2, 95, 98, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 338, 2, 3, 98, 101, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 344, 3, 4, 101, 104, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 350, 4, 5, 104, 107, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 356, 5, 6, 107, 110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 362, 6, 7, 110, 113, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 368, 11, 14, 95, 116, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 386, 14, 17, 98, 125, 134, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 404, 17, 20, 101, 134, 143, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 422, 20, 23, 104, 143, 152, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 440, 23, 26, 107, 152, 161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 458, 26, 29, 110, 161, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 476, 29, 32, 113, 170, 179, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 494, 41, 47, 125, 188, 206, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 530, 47, 53, 134, 206, 224, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 566, 53, 59, 143, 224, 242, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 602, 59, 65, 152, 242, 260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 638, 65, 71, 161, 260, 278, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 674, 71, 77, 170, 278, 296, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 710, 77, 83, 179, 296, 314, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 746, 95, 98, 332, 338, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 756, 98, 101, 338, 344, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 766, 101, 104, 344, 350, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 776, 104, 107, 350, 356, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 786, 107, 110, 356, 362, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 796, 116, 125, 332, 368, 386, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 826, 125, 134, 338, 386, 404, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 856, 134, 143, 344, 404, 422, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 886, 143, 152, 350, 422, 440, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 916, 152, 161, 356, 440, 458, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 946, 161, 170, 362, 458, 476, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 976, 188, 206, 386, 494, 530, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1036, 206, 224, 404, 530, 566, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1096, 224, 242, 422, 566, 602, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1156, 242, 260, 440, 602, 638, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1216, 260, 278, 458, 638, 674, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1276, 278, 296, 476, 674, 710, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1336, 332, 338, 746, 756, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1351, 338, 344, 756, 766, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1366, 344, 350, 766, 776, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1381, 350, 356, 776, 786, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1396, 368, 386, 746, 796, 826, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1441, 386, 404, 756, 826, 856, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1486, 404, 422, 766, 856, 886, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1531, 422, 440, 776, 886, 916, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1576, 440, 458, 786, 916, 946, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1621, 494, 530, 826, 976, 1036, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1711, 530, 566, 856, 1036, 1096, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1801, 566, 602, 886, 1096, 1156, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1891, 602, 638, 916, 1156, 1216, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1981, 638, 674, 946, 1216, 1276, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 2071, 746, 756, 1336, 1351, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 2092, 756, 766, 1351, 1366, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 2113, 766, 776, 1366, 1381, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2134, 796, 826, 1336, 1396, 1441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2197, 826, 856, 1351, 1441, 1486, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2260, 856, 886, 1366, 1486, 1531, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2323, 886, 916, 1381, 1531, 1576, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2386, 976, 1036, 1441, 1621, 1711, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2512, 1036, 1096, 1486, 1711, 1801, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2638, 1096, 1156, 1531, 1801, 1891, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2764, 1156, 1216, 1576, 1891, 1981, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 2890, 1336, 1351, 2071, 2092, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 2918, 1351, 1366, 2092, 2113, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 2946, 1396, 1441, 2071, 2134, 2197, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 3030, 1441, 1486, 2092, 2197, 2260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 3114, 1486, 1531, 2113, 2260, 2323, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 3198, 1621, 1711, 2197, 2386, 2512, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 3366, 1711, 1801, 2260, 2512, 2638, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 3534, 1801, 1891, 2323, 2638, 2764, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_skss(pbuffer, 3702, 2071, 2092, 2890, 2918, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 3738, 2134, 2197, 2890, 2946, 3030, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 3846, 2197, 2260, 2918, 3030, 3114, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 3954, 2386, 2512, 3030, 3198, 3366, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 4170, 2512, 2638, 3114, 3366, 3534, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsp(pbuffer, 4386, 2946, 3030, 3702, 3738, 3846, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slsd(pbuffer, 4521, 3198, 3366, 3846, 3954, 4170, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 368, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 494, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 54, pbuffer, 796, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 84, pbuffer, 976, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 144, pbuffer, 1396, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 189, pbuffer, 1621, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {1396, 1441});

                pbuffer.scale(2.0 * b_exp, {1621, 1711});

                pbuffer.scale(2.0 * b_exp, {2134, 2197});

                pbuffer.scale(2.0 * b_exp, {2386, 2512});

                pbuffer.scale(2.0 * b_exp, {2946, 3030});

                pbuffer.scale(2.0 * b_exp, {3198, 3366});

                t2cfunc::reduce(cbuffer, 279, pbuffer, 1396, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 324, pbuffer, 1621, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 414, pbuffer, 2134, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 477, pbuffer, 2386, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 603, pbuffer, 2946, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 687, pbuffer, 3198, 168, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {368, 386});

                pbuffer.scale(2.0 * a_exp, {494, 530});

                pbuffer.scale(2.0 * a_exp, {796, 826});

                pbuffer.scale(2.0 * a_exp, {976, 1036});

                pbuffer.scale(a_exp / b_exp, {1396, 1441});

                pbuffer.scale(a_exp / b_exp, {1621, 1711});

                pbuffer.scale(a_exp / b_exp, {2134, 2197});

                pbuffer.scale(a_exp / b_exp, {2386, 2512});

                pbuffer.scale(a_exp / b_exp, {2946, 3030});

                pbuffer.scale(a_exp / b_exp, {3198, 3366});

                t2cfunc::reduce(cbuffer, 855, pbuffer, 368, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 873, pbuffer, 494, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 909, pbuffer, 796, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 939, pbuffer, 976, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 999, pbuffer, 1396, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1044, pbuffer, 1621, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1134, pbuffer, 2134, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1197, pbuffer, 2386, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1323, pbuffer, 2946, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1407, pbuffer, 3198, 168, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {1396, 1441});

                pbuffer.scale(2.0 * b_exp, {1621, 1711});

                pbuffer.scale(2.0 * b_exp, {2134, 2197});

                pbuffer.scale(2.0 * b_exp, {2386, 2512});

                pbuffer.scale(2.0 * b_exp, {2946, 3030});

                pbuffer.scale(2.0 * b_exp, {3198, 3366});

                pbuffer.scale(4.0 * a_exp * b_exp, {3738, 3846});

                pbuffer.scale(4.0 * a_exp * b_exp, {3954, 4170});

                pbuffer.scale(4.0 * a_exp * b_exp, {4386, 4521});

                pbuffer.scale(4.0 * a_exp * b_exp, {4521, 4791});

                t2cfunc::reduce(cbuffer, 1575, pbuffer, 1396, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1620, pbuffer, 1621, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1710, pbuffer, 2134, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1773, pbuffer, 2386, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1899, pbuffer, 2946, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1983, pbuffer, 3198, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2151, pbuffer, 3738, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2259, pbuffer, 3954, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2475, pbuffer, 4386, 135, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2610, pbuffer, 4521, 270, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 0, cbuffer, 0, 18, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 54, cbuffer, 54, 84, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 414, cbuffer, 144, 189, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1521, cbuffer, 279, 324, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1656, cbuffer, 414, 477, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1845, cbuffer, 603, 687, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2097, cbuffer, 855, 873, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2151, cbuffer, 909, 939, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2511, cbuffer, 999, 1044, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3051, cbuffer, 1134, 1197, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3807, cbuffer, 1323, 1407, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 5787, cbuffer, 1575, 1620, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 5922, cbuffer, 1710, 1773, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6111, cbuffer, 1899, 1983, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6363, cbuffer, 2151, 2259, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6687, cbuffer, 2475, 2610, cfactors, 6, 0, 8);

            t4cfunc::ket_transform<1, 1>(skbuffer, 0, ckbuffer, 0, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 54, ckbuffer, 54, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 1224, ckbuffer, 414, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 46503, ckbuffer, 1521, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 46638, ckbuffer, 1656, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 46827, ckbuffer, 1845, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 47079, ckbuffer, 2097, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 47133, ckbuffer, 2151, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 47493, ckbuffer, 2511, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 48033, ckbuffer, 3051, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 48789, ckbuffer, 3807, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 52011, ckbuffer, 5787, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 52146, ckbuffer, 5922, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 52335, ckbuffer, 6111, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 52587, ckbuffer, 6363, 0, 7);

            t4cfunc::ket_transform<1, 1>(skbuffer, 52911, ckbuffer, 6687, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 7515, 54, 1224, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 50769, 47133, 47493, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 51039, 47493, 48033, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 51444, 48033, 48789, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 8595, 54, 50769, 51039, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 13050, 1224, 51039, 51444, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 24633, 7515, 8595, 13050, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 144, 0, 46503, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 1359, 54, 46638, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 2979, 1224, 46827, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 7785, 54, 144, 1359, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pgxx(skbuffer, 11835, 1224, 1359, 2979, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dfxx(skbuffer, 23013, 7515, 7785, 11835, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 47223, 47079, 52011, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 47628, 47133, 52146, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 48222, 47493, 52335, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 49041, 48033, 52587, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_skxx(skbuffer, 49797, 48789, 52911, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 414, 47133, 47223, 47628, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 1764, 47493, 47628, 48222, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 3546, 48033, 48222, 49041, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sixx(skbuffer, 5247, 48789, 49041, 49797, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 9405, 144, 50769, 414, 1764, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 14265, 1359, 51039, 1764, 3546, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_phxx(skbuffer, 17910, 2979, 51444, 3546, 5247, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 26253, 7785, 8595, 9405, 14265, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dgxx(skbuffer, 31113, 11835, 13050, 14265, 17910, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ffxx(skbuffer, 38403, 23013, 24633, 26253, 31113, r_ab, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 38403, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 441, skbuffer, 39303, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 882, skbuffer, 40203, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1323, skbuffer, 41103, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1764, skbuffer, 42003, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2205, skbuffer, 42903, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2646, skbuffer, 43803, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 3087, skbuffer, 44703, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 3528, skbuffer, 45603, 1, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 1, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFFPP_hpp */
