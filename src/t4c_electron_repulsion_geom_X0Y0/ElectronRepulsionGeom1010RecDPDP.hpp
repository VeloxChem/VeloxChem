#ifndef ElectronRepulsionGeom1010RecDPDP_hpp
#define ElectronRepulsionGeom1010RecDPDP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDP.hpp"
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DP|1/|r-r'||DP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dpdp(T& distributor,
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

    CSimdArray<double> pbuffer(2451, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1849, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(8643, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(10980, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2025, 1);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 9);

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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 9, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 12, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 21, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 24, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 27, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 30, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 33, 0, 1, 9, 12, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 39, 1, 2, 12, 15, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 45, 2, 3, 15, 18, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 51, 3, 4, 18, 21, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 57, 4, 5, 21, 24, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 63, 5, 6, 24, 27, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 69, 6, 7, 27, 30, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 75, 9, 12, 33, 39, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 12, 15, 39, 45, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 15, 18, 45, 51, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 18, 21, 51, 57, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 21, 24, 57, 63, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 24, 27, 63, 69, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 135, 33, 39, 75, 85, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 150, 39, 45, 85, 95, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 165, 45, 51, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 180, 51, 57, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 195, 57, 63, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 210, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 213, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 216, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 219, 1, 9, 12, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 228, 2, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 237, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 246, 4, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 255, 12, 33, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 273, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 291, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 309, 21, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 327, 39, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 357, 45, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 387, 51, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 417, 57, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 447, 85, 135, 150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 492, 95, 150, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 537, 105, 165, 180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 582, 115, 180, 195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 627, 1, 2, 210, 213, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 633, 2, 3, 213, 216, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 639, 9, 12, 210, 219, 228, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 657, 12, 15, 213, 228, 237, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 675, 15, 18, 216, 237, 246, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 693, 33, 39, 228, 255, 273, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 729, 39, 45, 237, 273, 291, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 765, 45, 51, 246, 291, 309, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 801, 75, 85, 273, 327, 357, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 861, 85, 95, 291, 357, 387, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 921, 95, 105, 309, 387, 417, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 981, 135, 150, 357, 447, 492, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1071, 150, 165, 387, 492, 537, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1161, 165, 180, 417, 537, 582, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1251, 210, 213, 627, 633, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1261, 219, 228, 627, 639, 657, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1291, 228, 237, 633, 657, 675, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1321, 255, 273, 657, 693, 729, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1381, 273, 291, 675, 729, 765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1441, 327, 357, 729, 801, 861, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1541, 357, 387, 765, 861, 921, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1641, 447, 492, 861, 981, 1071, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1791, 492, 537, 921, 1071, 1161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1941, 639, 657, 1251, 1261, 1291, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1986, 693, 729, 1291, 1321, 1381, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2076, 801, 861, 1381, 1441, 1541, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 2226, 981, 1071, 1541, 1641, 1791, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 219, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 255, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 27, pbuffer, 639, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 45, pbuffer, 693, 36, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {219, 228});

                pbuffer.scale(2.0 * a_exp, {255, 273});

                pbuffer.scale(2.0 * a_exp, {639, 657});

                pbuffer.scale(2.0 * a_exp, {693, 729});

                pbuffer.scale(2.0 * a_exp, {1261, 1291});

                pbuffer.scale(2.0 * a_exp, {1321, 1381});

                pbuffer.scale(2.0 * a_exp, {1941, 1986});

                pbuffer.scale(2.0 * a_exp, {1986, 2076});

                t2cfunc::reduce(cbuffer, 387, pbuffer, 219, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 396, pbuffer, 255, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 414, pbuffer, 639, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 432, pbuffer, 693, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 468, pbuffer, 1261, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 498, pbuffer, 1321, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 558, pbuffer, 1941, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 603, pbuffer, 1986, 90, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {219, 228});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {255, 273});

                pbuffer.scale(pfactors, 0, 2.0, {327, 357});

                pbuffer.scale(pfactors, 0, 2.0, {447, 492});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {639, 657});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {693, 729});

                pbuffer.scale(pfactors, 0, 2.0, {801, 861});

                pbuffer.scale(pfactors, 0, 2.0, {981, 1071});

                t2cfunc::reduce(cbuffer, 81, pbuffer, 219, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 255, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 108, pbuffer, 327, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 138, pbuffer, 447, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 183, pbuffer, 639, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 201, pbuffer, 693, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 237, pbuffer, 801, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 297, pbuffer, 981, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {219, 228});

                pbuffer.scale(2.0 * a_exp, {255, 273});

                pbuffer.scale(2.0 * a_exp, {327, 357});

                pbuffer.scale(2.0 * a_exp, {447, 492});

                pbuffer.scale(2.0 * a_exp, {639, 657});

                pbuffer.scale(2.0 * a_exp, {693, 729});

                pbuffer.scale(2.0 * a_exp, {801, 861});

                pbuffer.scale(2.0 * a_exp, {981, 1071});

                pbuffer.scale(pfactors, 0, 2.0, {1261, 1291});

                pbuffer.scale(pfactors, 0, 2.0, {1321, 1381});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1441, 1541});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1641, 1791});

                pbuffer.scale(pfactors, 0, 2.0, {1941, 1986});

                pbuffer.scale(pfactors, 0, 2.0, {1986, 2076});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2076, 2226});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2226, 2451});

                t2cfunc::reduce(cbuffer, 693, pbuffer, 219, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 702, pbuffer, 255, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 720, pbuffer, 327, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 750, pbuffer, 447, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 795, pbuffer, 639, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 813, pbuffer, 693, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 849, pbuffer, 801, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 909, pbuffer, 981, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 999, pbuffer, 1261, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1029, pbuffer, 1321, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1089, pbuffer, 1441, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1189, pbuffer, 1641, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1339, pbuffer, 1941, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1384, pbuffer, 1986, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1474, pbuffer, 2076, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1624, pbuffer, 2226, 225, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 171, cbuffer, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 945, cbuffer, 27, 45, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1980, cbuffer, 387, 396, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2754, cbuffer, 414, 432, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4188, cbuffer, 468, 498, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6483, cbuffer, 558, 603, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 81, 90, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 27, cbuffer, 90, 108, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 81, cbuffer, 108, 138, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 198, cbuffer, 0, 0, 27, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 279, cbuffer, 9, 27, 81, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 441, 171, 198, 279, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 603, cbuffer, 183, 201, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 657, cbuffer, 201, 237, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 765, cbuffer, 237, 297, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 999, cbuffer, 27, 603, 657, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1161, cbuffer, 45, 657, 765, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1485, 945, 999, 1161, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1809, cbuffer, 693, 702, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1836, cbuffer, 702, 720, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1890, cbuffer, 720, 750, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2007, cbuffer, 387, 1809, 1836, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2088, cbuffer, 396, 1836, 1890, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 2250, 1980, 2007, 2088, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2412, cbuffer, 795, 813, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2466, cbuffer, 813, 849, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 2574, cbuffer, 849, 909, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2808, cbuffer, 414, 2412, 2466, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2970, cbuffer, 432, 2466, 2574, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 3294, 2754, 2808, 2970, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3618, cbuffer, 999, 1029, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3708, cbuffer, 1029, 1089, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3888, cbuffer, 1089, 1189, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4278, cbuffer, 468, 3618, 3708, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4548, cbuffer, 498, 3708, 3888, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 5088, 4188, 4278, 4548, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5628, cbuffer, 1339, 1384, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5763, cbuffer, 1384, 1474, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6033, cbuffer, 1474, 1624, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 6618, cbuffer, 558, 5628, 5763, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 7023, cbuffer, 603, 5763, 6033, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 7833, 6483, 6618, 7023, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 0, ckbuffer, 441, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 45, ckbuffer, 495, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 90, ckbuffer, 549, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 540, ckbuffer, 1485, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 630, ckbuffer, 1593, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 720, ckbuffer, 1701, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 2970, ckbuffer, 0, 1, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 3105, ckbuffer, 162, 1, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 3240, ckbuffer, 324, 1, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 9450, ckbuffer, 2250, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 9495, ckbuffer, 2304, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 9540, ckbuffer, 2358, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 9585, ckbuffer, 3294, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 9675, ckbuffer, 3402, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 9765, ckbuffer, 3510, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 9855, ckbuffer, 5088, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 10005, ckbuffer, 5268, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 10155, ckbuffer, 5448, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 10305, ckbuffer, 7833, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 10530, ckbuffer, 8103, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 10755, ckbuffer, 8373, 0, 4);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 135, 9450, 9585, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 810, 9585, 9855, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 1620, 9855, 10305, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 3375, 0, 135, 810, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 4590, 540, 810, 1620, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 7020, 2970, 3375, 4590, r_ab, 2, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 7020, 2, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 225, skbuffer, 7290, 2, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 450, skbuffer, 7560, 2, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 675, skbuffer, 7830, 2, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 900, skbuffer, 8100, 2, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1125, skbuffer, 8370, 2, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1350, skbuffer, 8640, 2, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1575, skbuffer, 8910, 2, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1800, skbuffer, 9180, 2, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 2, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDPDP_hpp */
