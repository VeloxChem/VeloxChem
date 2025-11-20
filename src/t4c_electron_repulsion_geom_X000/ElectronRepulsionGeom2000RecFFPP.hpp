#ifndef ElectronRepulsionGeom2000RecFFPP_hpp
#define ElectronRepulsionGeom2000RecFFPP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDFXX.hpp"
#include "ElectronRepulsionContrRecDGXX.hpp"
#include "ElectronRepulsionContrRecDHXX.hpp"
#include "ElectronRepulsionContrRecDIXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecPIXX.hpp"
#include "ElectronRepulsionContrRecPKXX.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionGeom1000ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDGXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecFFXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPHXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"
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

/// @brief Computes d^(2)/dA^(2)(FF|1/|r-r'||PP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_ffpp(T& distributor,
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

    CSimdArray<double> cbuffer(2286, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(2286, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(39357, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2646, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 796, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 976, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 1396, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 135, pbuffer, 1621, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {796, 826});

                pbuffer.scale(2.0 * a_exp, {976, 1036});

                pbuffer.scale(2.0 * a_exp, {1396, 1441});

                pbuffer.scale(2.0 * a_exp, {1621, 1711});

                pbuffer.scale(2.0 * a_exp, {2134, 2197});

                pbuffer.scale(2.0 * a_exp, {2386, 2512});

                pbuffer.scale(2.0 * a_exp, {2946, 3030});

                pbuffer.scale(2.0 * a_exp, {3198, 3366});

                t2cfunc::reduce(cbuffer, 225, pbuffer, 796, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 255, pbuffer, 976, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 315, pbuffer, 1396, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 360, pbuffer, 1621, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 450, pbuffer, 2134, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 513, pbuffer, 2386, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 639, pbuffer, 2946, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 723, pbuffer, 3198, 168, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {796, 826});

                pbuffer.scale(2.0 * a_exp, {976, 1036});

                pbuffer.scale(2.0 * a_exp, {1396, 1441});

                pbuffer.scale(2.0 * a_exp, {1621, 1711});

                pbuffer.scale(2.0 * a_exp, {2134, 2197});

                pbuffer.scale(2.0 * a_exp, {2386, 2512});

                pbuffer.scale(2.0 * a_exp, {2946, 3030});

                pbuffer.scale(2.0 * a_exp, {3198, 3366});

                pbuffer.scale(4.0 * a_exp * a_exp, {3738, 3846});

                pbuffer.scale(4.0 * a_exp * a_exp, {3954, 4170});

                pbuffer.scale(4.0 * a_exp * a_exp, {4386, 4521});

                pbuffer.scale(4.0 * a_exp * a_exp, {4521, 4791});

                t2cfunc::reduce(cbuffer, 891, pbuffer, 796, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 921, pbuffer, 976, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 981, pbuffer, 1396, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1026, pbuffer, 1621, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1116, pbuffer, 2134, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1179, pbuffer, 2386, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1305, pbuffer, 2946, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1389, pbuffer, 3198, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1557, pbuffer, 3738, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1665, pbuffer, 3954, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1881, pbuffer, 4386, 135, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2016, pbuffer, 4521, 270, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 0, cbuffer, 0, 30, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 90, cbuffer, 90, 135, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 225, cbuffer, 225, 255, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 315, cbuffer, 315, 360, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 450, cbuffer, 450, 513, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 639, cbuffer, 639, 723, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 891, cbuffer, 891, 921, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 981, cbuffer, 981, 1026, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1116, cbuffer, 1116, 1179, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1305, cbuffer, 1305, 1389, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1557, cbuffer, 1557, 1665, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1881, cbuffer, 1881, 2016, cfactors, 6, 0, 8);

            t4cfunc::ket_transform<1, 1>(skbuffer, 0, ckbuffer, 0, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 630, ckbuffer, 90, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 29088, ckbuffer, 225, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 29178, ckbuffer, 315, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 29313, ckbuffer, 450, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 29502, ckbuffer, 639, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 30996, ckbuffer, 891, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 31086, ckbuffer, 981, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 31221, ckbuffer, 1116, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 31410, ckbuffer, 1305, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 31662, ckbuffer, 1557, 0, 7);

            t4cfunc::ket_transform<1, 1>(skbuffer, 31986, ckbuffer, 1881, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 4221, 0, 630, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 29754, 29088, 29178, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 30024, 29178, 29313, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 30429, 29313, 29502, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 32391, 30996, 31086, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 32661, 31086, 31221, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 33066, 31221, 31410, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pixx(skbuffer, 33633, 31410, 31662, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pkxx(skbuffer, 34389, 31662, 31986, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 35361, 32391, 32661, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_dgxx(skbuffer, 35901, 32661, 33066, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_dhxx(skbuffer, 36711, 33066, 33633, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_dixx(skbuffer, 37845, 33633, 34389, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 4491, 0, 29754, 30024, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 6921, 630, 30024, 30429, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 13968, 4221, 4491, 6921, r_ab, 1, 1);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 90, 29088, 35361, 3, 1, 1);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 765, 29178, 35901, 4, 1, 1);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 1575, 29313, 36711, 5, 1, 1);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 2709, 29502, 37845, 6, 1, 1);

            erirec::comp_bra_geom20_hrr_electron_repulsion_pfxx(skbuffer, 5301, 29754, 90, 765, r_ab, 1, 1);

            erirec::comp_bra_geom20_hrr_electron_repulsion_pgxx(skbuffer, 8136, 30024, 765, 1575, r_ab, 1, 1);

            erirec::comp_bra_geom20_hrr_electron_repulsion_phxx(skbuffer, 10566, 30429, 1575, 2709, r_ab, 1, 1);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dfxx(skbuffer, 15588, 4491, 5301, 8136, r_ab, 1, 1);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dgxx(skbuffer, 18828, 6921, 8136, 10566, r_ab, 1, 1);

            erirec::comp_bra_geom20_hrr_electron_repulsion_ffxx(skbuffer, 23688, 13968, 15588, 18828, r_ab, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 23688, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 441, skbuffer, 24588, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 882, skbuffer, 25488, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1323, skbuffer, 26388, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1764, skbuffer, 27288, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2205, skbuffer, 28188, 1, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 1, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecFFPP_hpp */
