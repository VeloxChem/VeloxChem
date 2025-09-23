#ifndef ElectronRepulsionGeom1010RecDPFS_hpp
#define ElectronRepulsionGeom1010RecDPFS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXDS.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXFS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
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
#include "ElectronRepulsionPrimRecSGSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DP|1/|r-r'||FS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dpfs(T& distributor,
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

    CSimdArray<double> pbuffer(2485, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1935, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(11610, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(5124, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 210, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 213, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 216, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 219, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 222, 1, 9, 12, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 231, 2, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 240, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 249, 4, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 258, 12, 33, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 276, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 294, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 312, 21, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 330, 39, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 360, 45, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 390, 51, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 420, 57, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 450, 85, 135, 150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 495, 95, 150, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 540, 105, 165, 180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 585, 115, 180, 195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 630, 0, 1, 210, 213, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 636, 1, 2, 213, 216, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 642, 2, 3, 216, 219, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 648, 9, 12, 213, 222, 231, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 666, 12, 15, 216, 231, 240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 684, 15, 18, 219, 240, 249, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 702, 33, 39, 231, 258, 276, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 738, 39, 45, 240, 276, 294, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 774, 45, 51, 249, 294, 312, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 810, 75, 85, 276, 330, 360, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 870, 85, 95, 294, 360, 390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 930, 95, 105, 312, 390, 420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 990, 135, 150, 360, 450, 495, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1080, 150, 165, 390, 495, 540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1170, 165, 180, 420, 540, 585, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1260, 210, 213, 630, 636, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1270, 213, 216, 636, 642, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1280, 222, 231, 636, 648, 666, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1310, 231, 240, 642, 666, 684, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1340, 258, 276, 666, 702, 738, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1400, 276, 294, 684, 738, 774, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1460, 330, 360, 738, 810, 870, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1560, 360, 390, 774, 870, 930, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1660, 450, 495, 870, 990, 1080, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1810, 495, 540, 930, 1080, 1170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1960, 630, 636, 1260, 1270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1975, 648, 666, 1270, 1280, 1310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2020, 702, 738, 1310, 1340, 1400, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2110, 810, 870, 1400, 1460, 1560, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 2260, 990, 1080, 1560, 1660, 1810, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 210, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 222, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12, pbuffer, 258, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 630, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 648, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 54, pbuffer, 702, 36, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {210, 213});

                pbuffer.scale(2.0 * a_exp, {222, 231});

                pbuffer.scale(2.0 * a_exp, {258, 276});

                pbuffer.scale(2.0 * a_exp, {630, 636});

                pbuffer.scale(2.0 * a_exp, {648, 666});

                pbuffer.scale(2.0 * a_exp, {702, 738});

                pbuffer.scale(2.0 * a_exp, {1260, 1270});

                pbuffer.scale(2.0 * a_exp, {1280, 1310});

                pbuffer.scale(2.0 * a_exp, {1340, 1400});

                pbuffer.scale(2.0 * a_exp, {1960, 1975});

                pbuffer.scale(2.0 * a_exp, {1975, 2020});

                pbuffer.scale(2.0 * a_exp, {2020, 2110});

                t2cfunc::reduce(cbuffer, 405, pbuffer, 210, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 408, pbuffer, 222, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 417, pbuffer, 258, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 435, pbuffer, 630, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 441, pbuffer, 648, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 459, pbuffer, 702, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 495, pbuffer, 1260, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 505, pbuffer, 1280, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 535, pbuffer, 1340, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 595, pbuffer, 1960, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 610, pbuffer, 1975, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 655, pbuffer, 2020, 90, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {210, 213});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {222, 231});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {258, 276});

                pbuffer.scale(pfactors, 0, 2.0, {330, 360});

                pbuffer.scale(pfactors, 0, 2.0, {450, 495});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {630, 636});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {648, 666});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {702, 738});

                pbuffer.scale(pfactors, 0, 2.0, {810, 870});

                pbuffer.scale(pfactors, 0, 2.0, {990, 1080});

                t2cfunc::reduce(cbuffer, 90, pbuffer, 210, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 93, pbuffer, 222, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 102, pbuffer, 258, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 120, pbuffer, 330, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 150, pbuffer, 450, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 195, pbuffer, 630, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 201, pbuffer, 648, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 219, pbuffer, 702, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 255, pbuffer, 810, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 315, pbuffer, 990, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {210, 213});

                pbuffer.scale(2.0 * a_exp, {222, 231});

                pbuffer.scale(2.0 * a_exp, {258, 276});

                pbuffer.scale(2.0 * a_exp, {330, 360});

                pbuffer.scale(2.0 * a_exp, {450, 495});

                pbuffer.scale(2.0 * a_exp, {630, 636});

                pbuffer.scale(2.0 * a_exp, {648, 666});

                pbuffer.scale(2.0 * a_exp, {702, 738});

                pbuffer.scale(2.0 * a_exp, {810, 870});

                pbuffer.scale(2.0 * a_exp, {990, 1080});

                pbuffer.scale(pfactors, 0, 2.0, {1260, 1270});

                pbuffer.scale(pfactors, 0, 2.0, {1280, 1310});

                pbuffer.scale(pfactors, 0, 2.0, {1340, 1400});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1460, 1560});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1660, 1810});

                pbuffer.scale(pfactors, 0, 2.0, {1960, 1975});

                pbuffer.scale(pfactors, 0, 2.0, {1975, 2020});

                pbuffer.scale(pfactors, 0, 2.0, {2020, 2110});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2110, 2260});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2260, 2485});

                t2cfunc::reduce(cbuffer, 745, pbuffer, 210, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 748, pbuffer, 222, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 757, pbuffer, 258, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 775, pbuffer, 330, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 805, pbuffer, 450, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 850, pbuffer, 630, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 856, pbuffer, 648, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 874, pbuffer, 702, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 910, pbuffer, 810, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 970, pbuffer, 990, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1060, pbuffer, 1260, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1070, pbuffer, 1280, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1100, pbuffer, 1340, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1160, pbuffer, 1460, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1260, pbuffer, 1660, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1410, pbuffer, 1960, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1425, pbuffer, 1975, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1470, pbuffer, 2020, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1560, pbuffer, 2110, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1710, pbuffer, 2260, 225, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 180, cbuffer, 0, 3, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 216, cbuffer, 3, 12, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 486, 180, 216, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1170, cbuffer, 30, 36, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1242, cbuffer, 36, 54, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 1782, 1170, 1242, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2610, cbuffer, 405, 408, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2646, cbuffer, 408, 417, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 2916, 2610, 2646, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 3600, cbuffer, 435, 441, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3672, cbuffer, 441, 459, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 4212, 3600, 3672, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 5460, cbuffer, 495, 505, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 5580, cbuffer, 505, 535, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 6480, 5460, 5580, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 8460, cbuffer, 595, 610, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 8640, cbuffer, 610, 655, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 9990, 8460, 8640, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 90, 93, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9, cbuffer, 93, 102, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 36, cbuffer, 102, 120, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 90, cbuffer, 120, 150, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 189, cbuffer, 0, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 243, cbuffer, 3, 9, 36, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 324, cbuffer, 12, 36, 90, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 504, 180, 189, 243, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 558, 216, 243, 324, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 720, 486, 504, 558, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 810, cbuffer, 195, 201, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 828, cbuffer, 201, 219, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 882, cbuffer, 219, 255, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 990, cbuffer, 255, 315, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1188, cbuffer, 30, 810, 828, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1296, cbuffer, 36, 828, 882, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1458, cbuffer, 54, 882, 990, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1818, 1170, 1188, 1296, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1926, 1242, 1296, 1458, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 2250, 1782, 1818, 1926, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2430, cbuffer, 745, 748, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2439, cbuffer, 748, 757, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2466, cbuffer, 757, 775, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 2520, cbuffer, 775, 805, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2619, cbuffer, 405, 2430, 2439, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2673, cbuffer, 408, 2439, 2466, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2754, cbuffer, 417, 2466, 2520, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 2934, 2610, 2619, 2673, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 2988, 2646, 2673, 2754, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 3150, 2916, 2934, 2988, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 3240, cbuffer, 850, 856, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3258, cbuffer, 856, 874, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3312, cbuffer, 874, 910, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3420, cbuffer, 910, 970, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 3618, cbuffer, 435, 3240, 3258, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3726, cbuffer, 441, 3258, 3312, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 3888, cbuffer, 459, 3312, 3420, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 4248, 3600, 3618, 3726, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 4356, 3672, 3726, 3888, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 4680, 4212, 4248, 4356, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 4860, cbuffer, 1060, 1070, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4890, cbuffer, 1070, 1100, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4980, cbuffer, 1100, 1160, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 5160, cbuffer, 1160, 1260, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 5490, cbuffer, 495, 4860, 4890, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 5670, cbuffer, 505, 4890, 4980, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 5940, cbuffer, 535, 4980, 5160, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 6540, 5460, 5490, 5670, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 6720, 5580, 5670, 5940, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 7260, 6480, 6540, 6720, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 7560, cbuffer, 1410, 1425, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 7605, cbuffer, 1425, 1470, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 7740, cbuffer, 1470, 1560, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 8010, cbuffer, 1560, 1710, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 8505, cbuffer, 595, 7560, 7605, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 8775, cbuffer, 610, 7605, 7740, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 9180, cbuffer, 655, 7740, 8010, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 10080, 8460, 8505, 8775, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 10350, 8640, 8775, 9180, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 11160, 9990, 10080, 10350, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 0, ckbuffer, 720, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21, ckbuffer, 750, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 42, ckbuffer, 780, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 252, ckbuffer, 2250, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 294, ckbuffer, 2310, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 336, ckbuffer, 2370, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 4410, ckbuffer, 3150, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 4431, ckbuffer, 3180, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 4452, ckbuffer, 3210, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 4473, ckbuffer, 4680, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 4515, ckbuffer, 4740, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 4557, ckbuffer, 4800, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 4599, ckbuffer, 7260, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 4669, ckbuffer, 7360, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 4739, ckbuffer, 7460, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 4809, ckbuffer, 11160, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 4914, ckbuffer, 11310, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 5019, ckbuffer, 11460, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1386, 0, 252, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1449, 21, 294, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1512, 42, 336, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 63, 4410, 4473, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 378, 4473, 4599, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 756, 4599, 4809, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 1575, 0, 63, 378, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 2142, 252, 378, 756, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 3276, 1386, 1575, 2142, r_ab, 3, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 3276, 3, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 105, skbuffer, 3402, 3, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 210, skbuffer, 3528, 3, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 315, skbuffer, 3654, 3, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 420, skbuffer, 3780, 3, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 525, skbuffer, 3906, 3, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 630, skbuffer, 4032, 3, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 735, skbuffer, 4158, 3, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 840, skbuffer, 4284, 3, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 3, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDPFS_hpp */
