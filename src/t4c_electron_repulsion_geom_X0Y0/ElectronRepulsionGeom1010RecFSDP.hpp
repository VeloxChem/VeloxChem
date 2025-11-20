#ifndef ElectronRepulsionGeom1010RecFSDP_hpp
#define ElectronRepulsionGeom1010RecFSDP_hpp

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
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSSXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FS|1/|r-r'||DP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fsdp(T& distributor,
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

    CSimdArray<double> cbuffer(1935, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(9045, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(14175, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(945, 1);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 9, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 9);
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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 9, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 219, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 255, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 639, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 54, pbuffer, 693, 36, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {9, 12});

                pbuffer.scale(2.0 * a_exp, {33, 39});

                pbuffer.scale(2.0 * a_exp, {219, 228});

                pbuffer.scale(2.0 * a_exp, {255, 273});

                pbuffer.scale(2.0 * a_exp, {639, 657});

                pbuffer.scale(2.0 * a_exp, {693, 729});

                pbuffer.scale(2.0 * a_exp, {1261, 1291});

                pbuffer.scale(2.0 * a_exp, {1321, 1381});

                pbuffer.scale(2.0 * a_exp, {1941, 1986});

                pbuffer.scale(2.0 * a_exp, {1986, 2076});

                t2cfunc::reduce(cbuffer, 430, pbuffer, 9, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 433, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 439, pbuffer, 219, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 448, pbuffer, 255, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 466, pbuffer, 639, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 484, pbuffer, 693, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 520, pbuffer, 1261, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 550, pbuffer, 1321, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 610, pbuffer, 1941, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 655, pbuffer, 1986, 90, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {9, 12});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {33, 39});

                pbuffer.scale(pfactors, 0, 2.0, {75, 85});

                pbuffer.scale(pfactors, 0, 2.0, {135, 150});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {219, 228});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {255, 273});

                pbuffer.scale(pfactors, 0, 2.0, {327, 357});

                pbuffer.scale(pfactors, 0, 2.0, {447, 492});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {639, 657});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {693, 729});

                pbuffer.scale(pfactors, 0, 2.0, {801, 861});

                pbuffer.scale(pfactors, 0, 2.0, {981, 1071});

                t2cfunc::reduce(cbuffer, 90, pbuffer, 9, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 93, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 99, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 109, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 124, pbuffer, 219, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 133, pbuffer, 255, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 151, pbuffer, 327, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 181, pbuffer, 447, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 226, pbuffer, 639, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 244, pbuffer, 693, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 280, pbuffer, 801, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 340, pbuffer, 981, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {9, 12});

                pbuffer.scale(2.0 * a_exp, {33, 39});

                pbuffer.scale(2.0 * a_exp, {75, 85});

                pbuffer.scale(2.0 * a_exp, {135, 150});

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

                t2cfunc::reduce(cbuffer, 745, pbuffer, 9, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 748, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 754, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 764, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 779, pbuffer, 219, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 788, pbuffer, 255, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 806, pbuffer, 327, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 836, pbuffer, 447, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 881, pbuffer, 639, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 899, pbuffer, 693, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 935, pbuffer, 801, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 995, pbuffer, 981, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1085, pbuffer, 1261, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1115, pbuffer, 1321, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1175, pbuffer, 1441, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1275, pbuffer, 1641, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1425, pbuffer, 1941, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1470, pbuffer, 1986, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1560, pbuffer, 2076, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1710, pbuffer, 2226, 225, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 57, cbuffer, 0, 3, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 372, cbuffer, 9, 18, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1146, cbuffer, 36, 54, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2067, cbuffer, 430, 433, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2382, cbuffer, 439, 448, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3156, cbuffer, 466, 484, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4590, cbuffer, 520, 550, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6885, cbuffer, 610, 655, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 90, 93, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9, cbuffer, 93, 99, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 27, cbuffer, 99, 109, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 66, cbuffer, 0, 0, 9, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 93, cbuffer, 3, 9, 27, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 147, 57, 66, 93, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 201, cbuffer, 124, 133, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 228, cbuffer, 133, 151, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 282, cbuffer, 151, 181, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 399, cbuffer, 9, 201, 228, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 480, cbuffer, 18, 228, 282, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 642, 372, 399, 480, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 804, cbuffer, 226, 244, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 858, cbuffer, 244, 280, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 966, cbuffer, 280, 340, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1200, cbuffer, 36, 804, 858, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1362, cbuffer, 54, 858, 966, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1686, 1146, 1200, 1362, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2010, cbuffer, 745, 748, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2019, cbuffer, 748, 754, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 2037, cbuffer, 754, 764, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2076, cbuffer, 430, 2010, 2019, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2103, cbuffer, 433, 2019, 2037, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 2157, 2067, 2076, 2103, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2211, cbuffer, 779, 788, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2238, cbuffer, 788, 806, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 2292, cbuffer, 806, 836, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2409, cbuffer, 439, 2211, 2238, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2490, cbuffer, 448, 2238, 2292, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 2652, 2382, 2409, 2490, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2814, cbuffer, 881, 899, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2868, cbuffer, 899, 935, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 2976, cbuffer, 935, 995, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3210, cbuffer, 466, 2814, 2868, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 3372, cbuffer, 484, 2868, 2976, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 3696, 3156, 3210, 3372, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4020, cbuffer, 1085, 1115, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4110, cbuffer, 1115, 1175, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 4290, cbuffer, 1175, 1275, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4680, cbuffer, 520, 4020, 4110, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4950, cbuffer, 550, 4110, 4290, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 5490, 4590, 4680, 4950, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 6030, cbuffer, 1425, 1470, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6165, cbuffer, 1470, 1560, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6435, cbuffer, 1560, 1710, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 7020, cbuffer, 610, 6030, 6165, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 7425, cbuffer, 655, 6165, 6435, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 8235, 6885, 7020, 7425, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 0, ckbuffer, 147, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 15, ckbuffer, 165, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 30, ckbuffer, 183, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 180, ckbuffer, 642, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 225, ckbuffer, 696, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 270, ckbuffer, 750, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 720, ckbuffer, 1686, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 810, ckbuffer, 1794, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 900, ckbuffer, 1902, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 12600, ckbuffer, 2157, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 12615, ckbuffer, 2175, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 12630, ckbuffer, 2193, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 12645, ckbuffer, 2652, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 12690, ckbuffer, 2706, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 12735, ckbuffer, 2760, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 12780, ckbuffer, 3696, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 12870, ckbuffer, 3804, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 12960, ckbuffer, 3912, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 13050, ckbuffer, 5490, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 13200, ckbuffer, 5670, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 13350, ckbuffer, 5850, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 13500, ckbuffer, 8235, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 13725, ckbuffer, 8505, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 13950, ckbuffer, 8775, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 3150, 0, 180, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 3195, 15, 225, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 3240, 30, 270, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 3690, 180, 720, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 3825, 225, 810, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 3960, 270, 900, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 7740, 3150, 3690, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 7830, 3195, 3825, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 7920, 3240, 3960, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 45, 12600, 12645, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 315, 12645, 12780, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 990, 12780, 13050, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 1800, 13050, 13500, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 3285, 0, 45, 315, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 4095, 180, 315, 990, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 5310, 720, 990, 1800, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 8010, 3150, 3285, 4095, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 8820, 3690, 4095, 5310, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fsxx(skbuffer, 11250, 7740, 8010, 8820, r_ab, 2, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 11250, 2, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 105, skbuffer, 11400, 2, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 210, skbuffer, 11550, 2, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 315, skbuffer, 11700, 2, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 420, skbuffer, 11850, 2, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 525, skbuffer, 12000, 2, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 630, skbuffer, 12150, 2, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 735, skbuffer, 12300, 2, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 840, skbuffer, 12450, 2, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 2, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFSDP_hpp */
