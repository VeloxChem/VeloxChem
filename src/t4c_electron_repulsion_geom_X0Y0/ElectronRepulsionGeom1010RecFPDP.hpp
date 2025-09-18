#ifndef ElectronRepulsionGeom1010RecFPDP_hpp
#define ElectronRepulsionGeom1010RecFPDP_hpp

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
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FP|1/|r-r'||DP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fpdp(T& distributor,
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

    CSimdArray<double> cbuffer(3182, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(14874, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(28980, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 257, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 302, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 27, pbuffer, 785, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 45, pbuffer, 857, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 81, pbuffer, 1621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 111, pbuffer, 1711, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {257, 266});

                pbuffer.scale(2.0 * a_exp, {302, 320});

                pbuffer.scale(2.0 * a_exp, {785, 803});

                pbuffer.scale(2.0 * a_exp, {857, 893});

                pbuffer.scale(2.0 * a_exp, {1621, 1651});

                pbuffer.scale(2.0 * a_exp, {1711, 1771});

                pbuffer.scale(2.0 * a_exp, {2656, 2701});

                pbuffer.scale(2.0 * a_exp, {2746, 2836});

                pbuffer.scale(2.0 * a_exp, {3676, 3739});

                pbuffer.scale(2.0 * a_exp, {3739, 3865});

                t2cfunc::reduce(cbuffer, 817, pbuffer, 257, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 826, pbuffer, 302, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 844, pbuffer, 785, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 862, pbuffer, 857, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 898, pbuffer, 1621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 928, pbuffer, 1711, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 988, pbuffer, 2656, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1033, pbuffer, 2746, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1123, pbuffer, 3676, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1186, pbuffer, 3739, 126, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {257, 266});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {302, 320});

                pbuffer.scale(pfactors, 0, 2.0, {392, 422});

                pbuffer.scale(pfactors, 0, 2.0, {542, 587});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {785, 803});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {857, 893});

                pbuffer.scale(pfactors, 0, 2.0, {1001, 1061});

                pbuffer.scale(pfactors, 0, 2.0, {1241, 1331});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1621, 1651});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1711, 1771});

                pbuffer.scale(pfactors, 0, 2.0, {1891, 1991});

                pbuffer.scale(pfactors, 0, 2.0, {2191, 2341});

                t2cfunc::reduce(cbuffer, 171, pbuffer, 257, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 180, pbuffer, 302, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 198, pbuffer, 392, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 228, pbuffer, 542, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 273, pbuffer, 785, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 291, pbuffer, 857, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 327, pbuffer, 1001, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 387, pbuffer, 1241, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 477, pbuffer, 1621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 507, pbuffer, 1711, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 567, pbuffer, 1891, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 667, pbuffer, 2191, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {257, 266});

                pbuffer.scale(2.0 * a_exp, {302, 320});

                pbuffer.scale(2.0 * a_exp, {392, 422});

                pbuffer.scale(2.0 * a_exp, {542, 587});

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

                t2cfunc::reduce(cbuffer, 1312, pbuffer, 257, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1321, pbuffer, 302, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1339, pbuffer, 392, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1369, pbuffer, 542, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1414, pbuffer, 785, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1432, pbuffer, 857, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1468, pbuffer, 1001, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1528, pbuffer, 1241, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1618, pbuffer, 1621, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1648, pbuffer, 1711, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1708, pbuffer, 1891, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1808, pbuffer, 2191, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1958, pbuffer, 2656, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2003, pbuffer, 2746, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2093, pbuffer, 2926, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2243, pbuffer, 3226, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2468, pbuffer, 3676, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2531, pbuffer, 3739, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2657, pbuffer, 3865, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2867, pbuffer, 4075, 315, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 171, cbuffer, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 945, cbuffer, 27, 45, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2379, cbuffer, 81, 111, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3990, cbuffer, 817, 826, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4764, cbuffer, 844, 862, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6198, cbuffer, 898, 928, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 8493, cbuffer, 988, 1033, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 11850, cbuffer, 1123, 1186, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 171, 180, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 27, cbuffer, 180, 198, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 81, cbuffer, 198, 228, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 198, cbuffer, 0, 0, 27, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 279, cbuffer, 9, 27, 81, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 441, 171, 198, 279, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 603, cbuffer, 273, 291, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 657, cbuffer, 291, 327, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 765, cbuffer, 327, 387, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 999, cbuffer, 27, 603, 657, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1161, cbuffer, 45, 657, 765, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1485, 945, 999, 1161, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1809, cbuffer, 477, 507, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1899, cbuffer, 507, 567, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 2079, cbuffer, 567, 667, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2469, cbuffer, 81, 1809, 1899, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2739, cbuffer, 111, 1899, 2079, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 3279, 2379, 2469, 2739, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3819, cbuffer, 1312, 1321, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3846, cbuffer, 1321, 1339, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3900, cbuffer, 1339, 1369, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4017, cbuffer, 817, 3819, 3846, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4098, cbuffer, 826, 3846, 3900, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 4260, 3990, 4017, 4098, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4422, cbuffer, 1414, 1432, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4476, cbuffer, 1432, 1468, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 4584, cbuffer, 1468, 1528, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4818, cbuffer, 844, 4422, 4476, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4980, cbuffer, 862, 4476, 4584, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 5304, 4764, 4818, 4980, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5628, cbuffer, 1618, 1648, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5718, cbuffer, 1648, 1708, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 5898, cbuffer, 1708, 1808, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 6288, cbuffer, 898, 5628, 5718, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 6558, cbuffer, 928, 5718, 5898, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 7098, 6198, 6288, 6558, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 7638, cbuffer, 1958, 2003, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 7773, cbuffer, 2003, 2093, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 8043, cbuffer, 2093, 2243, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 8628, cbuffer, 988, 7638, 7773, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 9033, cbuffer, 1033, 7773, 8043, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 9843, 8493, 8628, 9033, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 10653, cbuffer, 2468, 2531, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 10842, cbuffer, 2531, 2657, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 11220, cbuffer, 2657, 2867, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 12039, cbuffer, 1123, 10653, 10842, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 12606, cbuffer, 1186, 10842, 11220, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 13740, 11850, 12039, 12606, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 0, ckbuffer, 441, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 45, ckbuffer, 495, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 90, ckbuffer, 549, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 540, ckbuffer, 1485, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 630, ckbuffer, 1593, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 720, ckbuffer, 1701, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1620, ckbuffer, 3279, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1770, ckbuffer, 3459, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1920, ckbuffer, 3639, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 26505, ckbuffer, 4260, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 26550, ckbuffer, 4314, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 26595, ckbuffer, 4368, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 26640, ckbuffer, 5304, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 26730, ckbuffer, 5412, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 26820, ckbuffer, 5520, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 26910, ckbuffer, 7098, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 27060, ckbuffer, 7278, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 27210, ckbuffer, 7458, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 27360, ckbuffer, 9843, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 27585, ckbuffer, 10113, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 27810, ckbuffer, 10383, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 28035, ckbuffer, 13740, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 28350, ckbuffer, 14118, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 28665, ckbuffer, 14496, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 5445, 0, 540, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 5580, 45, 630, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 5715, 90, 720, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 7065, 540, 1620, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 7335, 630, 1770, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 7605, 720, 1920, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 14355, 5445, 7065, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 14625, 5580, 7335, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 14895, 5715, 7605, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 135, 26505, 26640, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 810, 26640, 26910, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 2070, 26910, 27360, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 3420, 27360, 28035, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 5850, 0, 135, 810, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 7875, 540, 810, 2070, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 10305, 1620, 2070, 3420, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 15165, 5445, 5850, 7875, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 17595, 7065, 7875, 10305, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fpxx(skbuffer, 22455, 14355, 15165, 17595, r_ab, 2, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 22455, 2, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 315, skbuffer, 22905, 2, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 630, skbuffer, 23355, 2, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 945, skbuffer, 23805, 2, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1260, skbuffer, 24255, 2, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1575, skbuffer, 24705, 2, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1890, skbuffer, 25155, 2, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2205, skbuffer, 25605, 2, 1);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2520, skbuffer, 26055, 2, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 2, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFPDP_hpp */
