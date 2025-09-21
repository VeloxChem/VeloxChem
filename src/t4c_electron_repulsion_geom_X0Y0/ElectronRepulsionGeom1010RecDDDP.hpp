#ifndef ElectronRepulsionGeom1010RecDDDP_hpp
#define ElectronRepulsionGeom1010RecDDDP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom1010ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DD|1/|r-r'||DP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dddp(T& distributor,
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

    CSimdArray<double> pbuffer(4390, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2924, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(13668, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(19395, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3375, 1);

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

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 10, 13, 37, 43, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 13, 16, 43, 49, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 16, 19, 49, 55, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 19, 22, 55, 61, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 22, 25, 61, 67, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 25, 28, 67, 73, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 28, 31, 73, 79, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 155, 37, 43, 85, 95, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 170, 43, 49, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 185, 49, 55, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 200, 55, 61, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 215, 61, 67, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 230, 67, 73, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 245, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 248, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 251, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 254, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 257, 1, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 266, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 275, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 284, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 293, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 302, 13, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 320, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 338, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 356, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 374, 25, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 392, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 422, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 452, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 482, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 512, 67, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 542, 95, 155, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 587, 105, 170, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 632, 115, 185, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 677, 125, 200, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 722, 135, 215, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 767, 1, 2, 245, 248, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 773, 2, 3, 248, 251, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 779, 3, 4, 251, 254, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 785, 10, 13, 245, 257, 266, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 803, 13, 16, 248, 266, 275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 821, 16, 19, 251, 275, 284, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 839, 19, 22, 254, 284, 293, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 857, 37, 43, 266, 302, 320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 893, 43, 49, 275, 320, 338, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 929, 49, 55, 284, 338, 356, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 965, 55, 61, 293, 356, 374, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1001, 85, 95, 320, 392, 422, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1061, 95, 105, 338, 422, 452, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1121, 105, 115, 356, 452, 482, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1181, 115, 125, 374, 482, 512, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1241, 155, 170, 422, 542, 587, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1331, 170, 185, 452, 587, 632, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1421, 185, 200, 482, 632, 677, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1511, 200, 215, 512, 677, 722, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1601, 245, 248, 767, 773, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1611, 248, 251, 773, 779, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1621, 257, 266, 767, 785, 803, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1651, 266, 275, 773, 803, 821, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1681, 275, 284, 779, 821, 839, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1711, 302, 320, 803, 857, 893, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1771, 320, 338, 821, 893, 929, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1831, 338, 356, 839, 929, 965, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1891, 392, 422, 893, 1001, 1061, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1991, 422, 452, 929, 1061, 1121, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2091, 452, 482, 965, 1121, 1181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2191, 542, 587, 1061, 1241, 1331, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2341, 587, 632, 1121, 1331, 1421, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2491, 632, 677, 1181, 1421, 1511, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2641, 767, 773, 1601, 1611, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2656, 785, 803, 1601, 1621, 1651, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2701, 803, 821, 1611, 1651, 1681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2746, 857, 893, 1651, 1711, 1771, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2836, 893, 929, 1681, 1771, 1831, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2926, 1001, 1061, 1771, 1891, 1991, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3076, 1061, 1121, 1831, 1991, 2091, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3226, 1241, 1331, 1991, 2191, 2341, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3451, 1331, 1421, 2091, 2341, 2491, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3676, 1621, 1651, 2641, 2656, 2701, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3739, 1711, 1771, 2701, 2746, 2836, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3865, 1891, 1991, 2836, 2926, 3076, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 4075, 2191, 2341, 3076, 3226, 3451, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 785, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 857, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 54, pbuffer, 1621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 84, pbuffer, 1711, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {785, 803});

                pbuffer.scale(2.0 * a_exp, {857, 893});

                pbuffer.scale(2.0 * a_exp, {1621, 1651});

                pbuffer.scale(2.0 * a_exp, {1711, 1771});

                pbuffer.scale(2.0 * a_exp, {2656, 2701});

                pbuffer.scale(2.0 * a_exp, {2746, 2836});

                pbuffer.scale(2.0 * a_exp, {3676, 3739});

                pbuffer.scale(2.0 * a_exp, {3739, 3865});

                t2cfunc::reduce(cbuffer, 688, pbuffer, 785, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 706, pbuffer, 857, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 742, pbuffer, 1621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 772, pbuffer, 1711, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 832, pbuffer, 2656, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 877, pbuffer, 2746, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 967, pbuffer, 3676, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1030, pbuffer, 3739, 126, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {785, 803});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {857, 893});

                pbuffer.scale(pfactors, 0, 2.0, {1001, 1061});

                pbuffer.scale(pfactors, 0, 2.0, {1241, 1331});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1621, 1651});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1711, 1771});

                pbuffer.scale(pfactors, 0, 2.0, {1891, 1991});

                pbuffer.scale(pfactors, 0, 2.0, {2191, 2341});

                t2cfunc::reduce(cbuffer, 144, pbuffer, 785, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 162, pbuffer, 857, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 198, pbuffer, 1001, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 258, pbuffer, 1241, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 348, pbuffer, 1621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 378, pbuffer, 1711, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 438, pbuffer, 1891, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 538, pbuffer, 2191, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {785, 803});

                pbuffer.scale(2.0 * a_exp, {857, 893});

                pbuffer.scale(2.0 * a_exp, {1001, 1061});

                pbuffer.scale(2.0 * a_exp, {1241, 1331});

                pbuffer.scale(2.0 * a_exp, {1621, 1651});

                pbuffer.scale(2.0 * a_exp, {1711, 1771});

                pbuffer.scale(2.0 * a_exp, {1891, 1991});

                pbuffer.scale(2.0 * a_exp, {2191, 2341});

                pbuffer.scale(pfactors, 0, 2.0, {2656, 2701});

                pbuffer.scale(pfactors, 0, 2.0, {2746, 2836});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2926, 3076});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3226, 3451});

                pbuffer.scale(pfactors, 0, 2.0, {3676, 3739});

                pbuffer.scale(pfactors, 0, 2.0, {3739, 3865});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3865, 4075});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {4075, 4390});

                t2cfunc::reduce(cbuffer, 1156, pbuffer, 785, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1174, pbuffer, 857, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1210, pbuffer, 1001, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1270, pbuffer, 1241, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1360, pbuffer, 1621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1390, pbuffer, 1711, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1450, pbuffer, 1891, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1550, pbuffer, 2191, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1700, pbuffer, 2656, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1745, pbuffer, 2746, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1835, pbuffer, 2926, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1985, pbuffer, 3226, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2210, pbuffer, 3676, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2273, pbuffer, 3739, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2399, pbuffer, 3865, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2609, pbuffer, 4075, 315, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 342, cbuffer, 0, 18, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1776, cbuffer, 54, 84, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3558, cbuffer, 688, 706, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4992, cbuffer, 742, 772, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 7287, cbuffer, 832, 877, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 10644, cbuffer, 967, 1030, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 144, 162, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 54, cbuffer, 162, 198, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 162, cbuffer, 198, 258, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 396, cbuffer, 0, 0, 54, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 558, cbuffer, 18, 54, 162, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 882, 342, 396, 558, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1206, cbuffer, 348, 378, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1296, cbuffer, 378, 438, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1476, cbuffer, 438, 538, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1866, cbuffer, 54, 1206, 1296, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2136, cbuffer, 84, 1296, 1476, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 2676, 1776, 1866, 2136, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3216, cbuffer, 1156, 1174, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3270, cbuffer, 1174, 1210, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3378, cbuffer, 1210, 1270, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3612, cbuffer, 688, 3216, 3270, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 3774, cbuffer, 706, 3270, 3378, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 4098, 3558, 3612, 3774, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4422, cbuffer, 1360, 1390, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4512, cbuffer, 1390, 1450, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 4692, cbuffer, 1450, 1550, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 5082, cbuffer, 742, 4422, 4512, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 5352, cbuffer, 772, 4512, 4692, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 5892, 4992, 5082, 5352, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 6432, cbuffer, 1700, 1745, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6567, cbuffer, 1745, 1835, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6837, cbuffer, 1835, 1985, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 7422, cbuffer, 832, 6432, 6567, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 7827, cbuffer, 877, 6567, 6837, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 8637, 7287, 7422, 7827, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9447, cbuffer, 2210, 2273, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9636, cbuffer, 2273, 2399, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 10014, cbuffer, 2399, 2609, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 10833, cbuffer, 967, 9447, 9636, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 11400, cbuffer, 1030, 9636, 10014, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 12534, 10644, 10833, 11400, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 0, ckbuffer, 882, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 90, ckbuffer, 990, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 180, ckbuffer, 1098, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1080, ckbuffer, 2676, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1230, ckbuffer, 2856, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1380, ckbuffer, 3036, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 17055, ckbuffer, 4098, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 17145, ckbuffer, 4206, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 17235, ckbuffer, 4314, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 17325, ckbuffer, 5892, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 17475, ckbuffer, 6072, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 17625, ckbuffer, 6252, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 17775, ckbuffer, 8637, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 18000, ckbuffer, 8907, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 18225, ckbuffer, 9177, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 18450, ckbuffer, 12534, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 18765, ckbuffer, 12912, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 19080, ckbuffer, 13290, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 4905, 0, 1080, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 5175, 90, 1230, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 5445, 180, 1380, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 270, 17055, 17325, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 1530, 17325, 17775, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 2880, 17775, 18450, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 5715, 0, 270, 1530, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 8145, 1080, 1530, 2880, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 12195, 4905, 5715, 8145, r_ab, 2, 1);

            t4cfunc::bra_transform<2, 2>(sbuffer, 0, skbuffer, 12195, 2, 1);

            t4cfunc::bra_transform<2, 2>(sbuffer, 375, skbuffer, 12735, 2, 1);

            t4cfunc::bra_transform<2, 2>(sbuffer, 750, skbuffer, 13275, 2, 1);

            t4cfunc::bra_transform<2, 2>(sbuffer, 1125, skbuffer, 13815, 2, 1);

            t4cfunc::bra_transform<2, 2>(sbuffer, 1500, skbuffer, 14355, 2, 1);

            t4cfunc::bra_transform<2, 2>(sbuffer, 1875, skbuffer, 14895, 2, 1);

            t4cfunc::bra_transform<2, 2>(sbuffer, 2250, skbuffer, 15435, 2, 1);

            t4cfunc::bra_transform<2, 2>(sbuffer, 2625, skbuffer, 15975, 2, 1);

            t4cfunc::bra_transform<2, 2>(sbuffer, 3000, skbuffer, 16515, 2, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 2, 2, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDDDP_hpp */
