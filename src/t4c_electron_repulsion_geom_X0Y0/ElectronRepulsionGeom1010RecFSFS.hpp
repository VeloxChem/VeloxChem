#ifndef ElectronRepulsionGeom1010RecFSFS_hpp
#define ElectronRepulsionGeom1010RecFSFS_hpp

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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FS|1/|r-r'||FS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fsfs(T& distributor,
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

    CSimdArray<double> cbuffer(2025, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(12150, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(6615, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(441, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1, pbuffer, 9, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 210, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13, pbuffer, 222, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 22, pbuffer, 258, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 40, pbuffer, 630, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 46, pbuffer, 648, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 64, pbuffer, 702, 36, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {9, 12});

                pbuffer.scale(2.0 * a_exp, {33, 39});

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

                t2cfunc::reduce(cbuffer, 450, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 451, pbuffer, 9, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 454, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 460, pbuffer, 210, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 463, pbuffer, 222, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 472, pbuffer, 258, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 490, pbuffer, 630, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 496, pbuffer, 648, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 514, pbuffer, 702, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 550, pbuffer, 1260, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 560, pbuffer, 1280, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 590, pbuffer, 1340, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 650, pbuffer, 1960, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 665, pbuffer, 1975, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 710, pbuffer, 2020, 90, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {0, 1});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {9, 12});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {33, 39});

                pbuffer.scale(pfactors, 0, 2.0, {75, 85});

                pbuffer.scale(pfactors, 0, 2.0, {135, 150});

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

                t2cfunc::reduce(cbuffer, 100, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 101, pbuffer, 9, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 104, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 110, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 120, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 135, pbuffer, 210, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 138, pbuffer, 222, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 147, pbuffer, 258, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 165, pbuffer, 330, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 195, pbuffer, 450, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 240, pbuffer, 630, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 246, pbuffer, 648, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 264, pbuffer, 702, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 300, pbuffer, 810, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 360, pbuffer, 990, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {9, 12});

                pbuffer.scale(2.0 * a_exp, {33, 39});

                pbuffer.scale(2.0 * a_exp, {75, 85});

                pbuffer.scale(2.0 * a_exp, {135, 150});

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

                t2cfunc::reduce(cbuffer, 800, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 801, pbuffer, 9, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 804, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 810, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 820, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 835, pbuffer, 210, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 838, pbuffer, 222, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 847, pbuffer, 258, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 865, pbuffer, 330, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 895, pbuffer, 450, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 940, pbuffer, 630, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 946, pbuffer, 648, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 964, pbuffer, 702, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1000, pbuffer, 810, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1060, pbuffer, 990, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1150, pbuffer, 1260, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1160, pbuffer, 1280, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1190, pbuffer, 1340, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1250, pbuffer, 1460, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1350, pbuffer, 1660, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1500, pbuffer, 1960, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1515, pbuffer, 1975, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1560, pbuffer, 2020, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1650, pbuffer, 2110, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1800, pbuffer, 2260, 225, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 60, cbuffer, 0, 1, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 72, cbuffer, 1, 4, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 162, 60, 72, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 450, cbuffer, 10, 13, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 486, cbuffer, 13, 22, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 756, 450, 486, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1440, cbuffer, 40, 46, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1512, cbuffer, 46, 64, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 2052, 1440, 1512, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2760, cbuffer, 450, 451, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2772, cbuffer, 451, 454, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 2862, 2760, 2772, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 3150, cbuffer, 460, 463, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3186, cbuffer, 463, 472, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 3456, 3150, 3186, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 4140, cbuffer, 490, 496, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4212, cbuffer, 496, 514, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 4752, 4140, 4212, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 6000, cbuffer, 550, 560, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6120, cbuffer, 560, 590, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 7020, 6000, 6120, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 9000, cbuffer, 650, 665, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 9180, cbuffer, 665, 710, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 10530, 9000, 9180, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 100, 101, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3, cbuffer, 101, 104, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 12, cbuffer, 104, 110, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 30, cbuffer, 110, 120, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 63, cbuffer, 0, 0, 3, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 81, cbuffer, 1, 3, 12, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 108, cbuffer, 4, 12, 30, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 168, 60, 63, 81, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 186, 72, 81, 108, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 240, 162, 168, 186, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 270, cbuffer, 135, 138, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 279, cbuffer, 138, 147, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 306, cbuffer, 147, 165, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 360, cbuffer, 165, 195, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 459, cbuffer, 10, 270, 279, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 513, cbuffer, 13, 279, 306, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 594, cbuffer, 22, 306, 360, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 774, 450, 459, 513, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 828, 486, 513, 594, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 990, 756, 774, 828, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1080, cbuffer, 240, 246, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1098, cbuffer, 246, 264, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1152, cbuffer, 264, 300, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1260, cbuffer, 300, 360, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1458, cbuffer, 40, 1080, 1098, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1566, cbuffer, 46, 1098, 1152, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1728, cbuffer, 64, 1152, 1260, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 2088, 1440, 1458, 1566, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 2196, 1512, 1566, 1728, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 2520, 2052, 2088, 2196, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2700, cbuffer, 800, 801, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2703, cbuffer, 801, 804, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2712, cbuffer, 804, 810, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 2730, cbuffer, 810, 820, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2763, cbuffer, 450, 2700, 2703, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2781, cbuffer, 451, 2703, 2712, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2808, cbuffer, 454, 2712, 2730, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 2868, 2760, 2763, 2781, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 2886, 2772, 2781, 2808, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 2940, 2862, 2868, 2886, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2970, cbuffer, 835, 838, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2979, cbuffer, 838, 847, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3006, cbuffer, 847, 865, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3060, cbuffer, 865, 895, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 3159, cbuffer, 460, 2970, 2979, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3213, cbuffer, 463, 2979, 3006, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 3294, cbuffer, 472, 3006, 3060, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 3474, 3150, 3159, 3213, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 3528, 3186, 3213, 3294, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 3690, 3456, 3474, 3528, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 3780, cbuffer, 940, 946, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3798, cbuffer, 946, 964, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3852, cbuffer, 964, 1000, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3960, cbuffer, 1000, 1060, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 4158, cbuffer, 490, 3780, 3798, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4266, cbuffer, 496, 3798, 3852, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4428, cbuffer, 514, 3852, 3960, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 4788, 4140, 4158, 4266, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 4896, 4212, 4266, 4428, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 5220, 4752, 4788, 4896, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 5400, cbuffer, 1150, 1160, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5430, cbuffer, 1160, 1190, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5520, cbuffer, 1190, 1250, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 5700, cbuffer, 1250, 1350, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 6030, cbuffer, 550, 5400, 5430, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 6210, cbuffer, 560, 5430, 5520, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 6480, cbuffer, 590, 5520, 5700, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 7080, 6000, 6030, 6210, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 7260, 6120, 6210, 6480, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 7800, 7020, 7080, 7260, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 8100, cbuffer, 1500, 1515, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 8145, cbuffer, 1515, 1560, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 8280, cbuffer, 1560, 1650, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 8550, cbuffer, 1650, 1800, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 9045, cbuffer, 650, 8100, 8145, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 9315, cbuffer, 665, 8145, 8280, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 9720, cbuffer, 710, 8280, 8550, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 10620, 9000, 9045, 9315, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 10890, 9180, 9315, 9720, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 11700, 10530, 10620, 10890, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 0, ckbuffer, 240, 0, 0);

            t4cfunc::ket_transform<3, 0>(skbuffer, 7, ckbuffer, 250, 0, 0);

            t4cfunc::ket_transform<3, 0>(skbuffer, 14, ckbuffer, 260, 0, 0);

            t4cfunc::ket_transform<3, 0>(skbuffer, 84, ckbuffer, 990, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 105, ckbuffer, 1020, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 126, ckbuffer, 1050, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 336, ckbuffer, 2520, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 378, ckbuffer, 2580, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 420, ckbuffer, 2640, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 5880, ckbuffer, 2940, 0, 0);

            t4cfunc::ket_transform<3, 0>(skbuffer, 5887, ckbuffer, 2950, 0, 0);

            t4cfunc::ket_transform<3, 0>(skbuffer, 5894, ckbuffer, 2960, 0, 0);

            t4cfunc::ket_transform<3, 0>(skbuffer, 5901, ckbuffer, 3690, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 5922, ckbuffer, 3720, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 5943, ckbuffer, 3750, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 5964, ckbuffer, 5220, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 6006, ckbuffer, 5280, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 6048, ckbuffer, 5340, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 6090, ckbuffer, 7800, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 6160, ckbuffer, 7900, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 6230, ckbuffer, 8000, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 6300, ckbuffer, 11700, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 6405, ckbuffer, 11850, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 6510, ckbuffer, 12000, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 1470, 0, 84, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 1491, 7, 105, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 1512, 14, 126, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1722, 84, 336, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1785, 105, 378, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1848, 126, 420, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 3612, 1470, 1722, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 3654, 1491, 1785, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 3696, 1512, 1848, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 21, 5880, 5901, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 147, 5901, 5964, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 462, 5964, 6090, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 840, 6090, 6300, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 1533, 0, 21, 147, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 1911, 84, 147, 462, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 2478, 336, 462, 840, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 3738, 1470, 1533, 1911, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 4116, 1722, 1911, 2478, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fsxx(skbuffer, 5250, 3612, 3738, 4116, r_ab, 3, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 5250, 3, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 49, skbuffer, 5320, 3, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 98, skbuffer, 5390, 3, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 147, skbuffer, 5460, 3, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 196, skbuffer, 5530, 3, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 245, skbuffer, 5600, 3, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 294, skbuffer, 5670, 3, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 343, skbuffer, 5740, 3, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 392, skbuffer, 5810, 3, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 3, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFSFS_hpp */
