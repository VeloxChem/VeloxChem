#ifndef ElectronRepulsionGeom1010RecDDFS_hpp
#define ElectronRepulsionGeom1010RecDDFS_hpp

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
#include "ElectronRepulsionPrimRecSHSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DD|1/|r-r'||FS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_ddfs(T& distributor,
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

    CSimdArray<double> pbuffer(4445, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(3060, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(18360, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(9051, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1575, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 245, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 248, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 251, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 254, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 257, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 260, 1, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 269, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 278, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 287, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 296, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 305, 13, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 323, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 341, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 359, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 377, 25, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 395, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 425, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 455, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 485, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 515, 67, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 545, 95, 155, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 590, 105, 170, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 635, 115, 185, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 680, 125, 200, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 725, 135, 215, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 770, 0, 1, 245, 248, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 776, 1, 2, 248, 251, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 782, 2, 3, 251, 254, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 788, 3, 4, 254, 257, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 794, 10, 13, 248, 260, 269, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 812, 13, 16, 251, 269, 278, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 830, 16, 19, 254, 278, 287, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 848, 19, 22, 257, 287, 296, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 866, 37, 43, 269, 305, 323, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 902, 43, 49, 278, 323, 341, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 938, 49, 55, 287, 341, 359, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 974, 55, 61, 296, 359, 377, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1010, 85, 95, 323, 395, 425, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1070, 95, 105, 341, 425, 455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1130, 105, 115, 359, 455, 485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1190, 115, 125, 377, 485, 515, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1250, 155, 170, 425, 545, 590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1340, 170, 185, 455, 590, 635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1430, 185, 200, 485, 635, 680, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1520, 200, 215, 515, 680, 725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1610, 245, 248, 770, 776, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1620, 248, 251, 776, 782, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1630, 251, 254, 782, 788, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1640, 260, 269, 776, 794, 812, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1670, 269, 278, 782, 812, 830, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1700, 278, 287, 788, 830, 848, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1730, 305, 323, 812, 866, 902, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1790, 323, 341, 830, 902, 938, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1850, 341, 359, 848, 938, 974, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1910, 395, 425, 902, 1010, 1070, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2010, 425, 455, 938, 1070, 1130, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2110, 455, 485, 974, 1130, 1190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2210, 545, 590, 1070, 1250, 1340, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2360, 590, 635, 1130, 1340, 1430, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2510, 635, 680, 1190, 1430, 1520, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2660, 770, 776, 1610, 1620, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2675, 776, 782, 1620, 1630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2690, 794, 812, 1620, 1640, 1670, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2735, 812, 830, 1630, 1670, 1700, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2780, 866, 902, 1670, 1730, 1790, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2870, 902, 938, 1700, 1790, 1850, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2960, 1010, 1070, 1790, 1910, 2010, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3110, 1070, 1130, 1850, 2010, 2110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3260, 1250, 1340, 2010, 2210, 2360, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3485, 1340, 1430, 2110, 2360, 2510, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 3710, 1610, 1620, 2660, 2675, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3731, 1640, 1670, 2675, 2690, 2735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3794, 1730, 1790, 2735, 2780, 2870, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3920, 1910, 2010, 2870, 2960, 3110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 4130, 2210, 2360, 3110, 3260, 3485, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 770, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 794, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 24, pbuffer, 866, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 1610, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 70, pbuffer, 1640, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 1730, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {770, 776});

                pbuffer.scale(2.0 * a_exp, {794, 812});

                pbuffer.scale(2.0 * a_exp, {866, 902});

                pbuffer.scale(2.0 * a_exp, {1610, 1620});

                pbuffer.scale(2.0 * a_exp, {1640, 1670});

                pbuffer.scale(2.0 * a_exp, {1730, 1790});

                pbuffer.scale(2.0 * a_exp, {2660, 2675});

                pbuffer.scale(2.0 * a_exp, {2690, 2735});

                pbuffer.scale(2.0 * a_exp, {2780, 2870});

                pbuffer.scale(2.0 * a_exp, {3710, 3731});

                pbuffer.scale(2.0 * a_exp, {3731, 3794});

                pbuffer.scale(2.0 * a_exp, {3794, 3920});

                t2cfunc::reduce(cbuffer, 720, pbuffer, 770, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 726, pbuffer, 794, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 744, pbuffer, 866, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 780, pbuffer, 1610, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 790, pbuffer, 1640, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 820, pbuffer, 1730, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 880, pbuffer, 2660, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 895, pbuffer, 2690, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 940, pbuffer, 2780, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1030, pbuffer, 3710, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1051, pbuffer, 3731, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1114, pbuffer, 3794, 126, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {770, 776});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {794, 812});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {866, 902});

                pbuffer.scale(pfactors, 0, 2.0, {1010, 1070});

                pbuffer.scale(pfactors, 0, 2.0, {1250, 1340});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1610, 1620});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1640, 1670});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1730, 1790});

                pbuffer.scale(pfactors, 0, 2.0, {1910, 2010});

                pbuffer.scale(pfactors, 0, 2.0, {2210, 2360});

                t2cfunc::reduce(cbuffer, 160, pbuffer, 770, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 166, pbuffer, 794, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 184, pbuffer, 866, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 220, pbuffer, 1010, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 280, pbuffer, 1250, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 370, pbuffer, 1610, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 380, pbuffer, 1640, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 410, pbuffer, 1730, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 470, pbuffer, 1910, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 570, pbuffer, 2210, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {770, 776});

                pbuffer.scale(2.0 * a_exp, {794, 812});

                pbuffer.scale(2.0 * a_exp, {866, 902});

                pbuffer.scale(2.0 * a_exp, {1010, 1070});

                pbuffer.scale(2.0 * a_exp, {1250, 1340});

                pbuffer.scale(2.0 * a_exp, {1610, 1620});

                pbuffer.scale(2.0 * a_exp, {1640, 1670});

                pbuffer.scale(2.0 * a_exp, {1730, 1790});

                pbuffer.scale(2.0 * a_exp, {1910, 2010});

                pbuffer.scale(2.0 * a_exp, {2210, 2360});

                pbuffer.scale(pfactors, 0, 2.0, {2660, 2675});

                pbuffer.scale(pfactors, 0, 2.0, {2690, 2735});

                pbuffer.scale(pfactors, 0, 2.0, {2780, 2870});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2960, 3110});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3260, 3485});

                pbuffer.scale(pfactors, 0, 2.0, {3710, 3731});

                pbuffer.scale(pfactors, 0, 2.0, {3731, 3794});

                pbuffer.scale(pfactors, 0, 2.0, {3794, 3920});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3920, 4130});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {4130, 4445});

                t2cfunc::reduce(cbuffer, 1240, pbuffer, 770, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1246, pbuffer, 794, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1264, pbuffer, 866, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1300, pbuffer, 1010, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1360, pbuffer, 1250, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1450, pbuffer, 1610, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1460, pbuffer, 1640, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1490, pbuffer, 1730, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1550, pbuffer, 1910, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1650, pbuffer, 2210, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1800, pbuffer, 2660, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1815, pbuffer, 2690, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1860, pbuffer, 2780, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1950, pbuffer, 2960, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2100, pbuffer, 3260, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2325, pbuffer, 3710, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2346, pbuffer, 3731, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2409, pbuffer, 3794, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2535, pbuffer, 3920, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2745, pbuffer, 4130, 315, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 360, cbuffer, 0, 6, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 432, cbuffer, 6, 24, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 972, 360, 432, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2220, cbuffer, 60, 70, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2340, cbuffer, 70, 100, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 3240, 2220, 2340, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 4680, cbuffer, 720, 726, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4752, cbuffer, 726, 744, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 5292, 4680, 4752, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 6540, cbuffer, 780, 790, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6660, cbuffer, 790, 820, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 7560, 6540, 6660, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 9540, cbuffer, 880, 895, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 9720, cbuffer, 895, 940, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 11070, 9540, 9720, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 13950, cbuffer, 1030, 1051, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 14202, cbuffer, 1051, 1114, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 16092, 13950, 14202, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 160, 166, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 18, cbuffer, 166, 184, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 72, cbuffer, 184, 220, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 180, cbuffer, 220, 280, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 378, cbuffer, 0, 0, 18, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 486, cbuffer, 6, 18, 72, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 648, cbuffer, 24, 72, 180, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1008, 360, 378, 486, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1116, 432, 486, 648, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 1440, 972, 1008, 1116, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1620, cbuffer, 370, 380, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1650, cbuffer, 380, 410, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1740, cbuffer, 410, 470, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1920, cbuffer, 470, 570, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2250, cbuffer, 60, 1620, 1650, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2430, cbuffer, 70, 1650, 1740, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2700, cbuffer, 100, 1740, 1920, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 3300, 2220, 2250, 2430, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 3480, 2340, 2430, 2700, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 4020, 3240, 3300, 3480, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 4320, cbuffer, 1240, 1246, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4338, cbuffer, 1246, 1264, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4392, cbuffer, 1264, 1300, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 4500, cbuffer, 1300, 1360, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 4698, cbuffer, 720, 4320, 4338, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4806, cbuffer, 726, 4338, 4392, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4968, cbuffer, 744, 4392, 4500, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 5328, 4680, 4698, 4806, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 5436, 4752, 4806, 4968, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 5760, 5292, 5328, 5436, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 5940, cbuffer, 1450, 1460, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5970, cbuffer, 1460, 1490, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6060, cbuffer, 1490, 1550, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6240, cbuffer, 1550, 1650, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 6570, cbuffer, 780, 5940, 5970, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 6750, cbuffer, 790, 5970, 6060, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 7020, cbuffer, 820, 6060, 6240, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 7620, 6540, 6570, 6750, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 7800, 6660, 6750, 7020, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 8340, 7560, 7620, 7800, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 8640, cbuffer, 1800, 1815, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 8685, cbuffer, 1815, 1860, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 8820, cbuffer, 1860, 1950, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 9090, cbuffer, 1950, 2100, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 9585, cbuffer, 880, 8640, 8685, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 9855, cbuffer, 895, 8685, 8820, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 10260, cbuffer, 940, 8820, 9090, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 11160, 9540, 9585, 9855, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 11430, 9720, 9855, 10260, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 12240, 11070, 11160, 11430, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 12690, cbuffer, 2325, 2346, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 12753, cbuffer, 2346, 2409, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 12942, cbuffer, 2409, 2535, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 13320, cbuffer, 2535, 2745, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 14013, cbuffer, 1030, 12690, 12753, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 14391, cbuffer, 1051, 12753, 12942, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 14958, cbuffer, 1114, 12942, 13320, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 16218, 13950, 14013, 14391, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 16596, 14202, 14391, 14958, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 17730, 16092, 16218, 16596, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 0, ckbuffer, 1440, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 42, ckbuffer, 1500, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 84, ckbuffer, 1560, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 504, ckbuffer, 4020, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 574, ckbuffer, 4120, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 644, ckbuffer, 4220, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 7959, ckbuffer, 5760, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 8001, ckbuffer, 5820, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 8043, ckbuffer, 5880, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 8085, ckbuffer, 8340, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 8155, ckbuffer, 8440, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 8225, ckbuffer, 8540, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 8295, ckbuffer, 12240, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 8400, ckbuffer, 12390, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 8505, ckbuffer, 12540, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 8610, ckbuffer, 17730, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 8757, ckbuffer, 17940, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 8904, ckbuffer, 18150, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2289, 0, 504, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2415, 42, 574, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2541, 84, 644, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 126, 7959, 8085, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 714, 8085, 8295, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 1344, 8295, 8610, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 2667, 0, 126, 714, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 3801, 504, 714, 1344, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 5691, 2289, 2667, 3801, r_ab, 3, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 0, skbuffer, 5691, 3, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 175, skbuffer, 5943, 3, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 350, skbuffer, 6195, 3, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 525, skbuffer, 6447, 3, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 700, skbuffer, 6699, 3, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 875, skbuffer, 6951, 3, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 1050, skbuffer, 7203, 3, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 1225, skbuffer, 7455, 3, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 1400, skbuffer, 7707, 3, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 2, 3, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDDFS_hpp */
