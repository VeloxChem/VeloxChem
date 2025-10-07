#ifndef ElectronRepulsionGeom1010RecFPFS_hpp
#define ElectronRepulsionGeom1010RecFPFS_hpp

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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FP|1/|r-r'||FS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fpfs(T& distributor,
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

    CSimdArray<double> cbuffer(3330, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(19980, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(13524, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1323, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 245, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 260, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12, pbuffer, 305, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 770, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 36, pbuffer, 794, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 54, pbuffer, 866, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 1610, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 1640, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 130, pbuffer, 1730, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {245, 248});

                pbuffer.scale(2.0 * a_exp, {260, 269});

                pbuffer.scale(2.0 * a_exp, {305, 323});

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

                t2cfunc::reduce(cbuffer, 855, pbuffer, 245, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 858, pbuffer, 260, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 867, pbuffer, 305, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 885, pbuffer, 770, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 891, pbuffer, 794, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 909, pbuffer, 866, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 945, pbuffer, 1610, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 955, pbuffer, 1640, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 985, pbuffer, 1730, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1045, pbuffer, 2660, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1060, pbuffer, 2690, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1105, pbuffer, 2780, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1195, pbuffer, 3710, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1216, pbuffer, 3731, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1279, pbuffer, 3794, 126, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {245, 248});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {260, 269});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {305, 323});

                pbuffer.scale(pfactors, 0, 2.0, {395, 425});

                pbuffer.scale(pfactors, 0, 2.0, {545, 590});

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

                t2cfunc::reduce(cbuffer, 190, pbuffer, 245, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 193, pbuffer, 260, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 202, pbuffer, 305, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 220, pbuffer, 395, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 250, pbuffer, 545, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 295, pbuffer, 770, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 301, pbuffer, 794, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 319, pbuffer, 866, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 355, pbuffer, 1010, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 415, pbuffer, 1250, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 505, pbuffer, 1610, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 515, pbuffer, 1640, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 545, pbuffer, 1730, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 605, pbuffer, 1910, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 705, pbuffer, 2210, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {245, 248});

                pbuffer.scale(2.0 * a_exp, {260, 269});

                pbuffer.scale(2.0 * a_exp, {305, 323});

                pbuffer.scale(2.0 * a_exp, {395, 425});

                pbuffer.scale(2.0 * a_exp, {545, 590});

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

                t2cfunc::reduce(cbuffer, 1405, pbuffer, 245, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1408, pbuffer, 260, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1417, pbuffer, 305, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1435, pbuffer, 395, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1465, pbuffer, 545, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1510, pbuffer, 770, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1516, pbuffer, 794, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1534, pbuffer, 866, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1570, pbuffer, 1010, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1630, pbuffer, 1250, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1720, pbuffer, 1610, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1730, pbuffer, 1640, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1760, pbuffer, 1730, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1820, pbuffer, 1910, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1920, pbuffer, 2210, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2070, pbuffer, 2660, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2085, pbuffer, 2690, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2130, pbuffer, 2780, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2220, pbuffer, 2960, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2370, pbuffer, 3260, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2595, pbuffer, 3710, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2616, pbuffer, 3731, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2679, pbuffer, 3794, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2805, pbuffer, 3920, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3015, pbuffer, 4130, 315, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 180, cbuffer, 0, 3, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 216, cbuffer, 3, 12, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 486, 180, 216, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1170, cbuffer, 30, 36, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1242, cbuffer, 36, 54, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 1782, 1170, 1242, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 3030, cbuffer, 90, 100, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3150, cbuffer, 100, 130, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 4050, 3030, 3150, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 5310, cbuffer, 855, 858, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 5346, cbuffer, 858, 867, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 5616, 5310, 5346, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 6300, cbuffer, 885, 891, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6372, cbuffer, 891, 909, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 6912, 6300, 6372, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 8160, cbuffer, 945, 955, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 8280, cbuffer, 955, 985, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 9180, 8160, 8280, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 11160, cbuffer, 1045, 1060, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 11340, cbuffer, 1060, 1105, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 12690, 11160, 11340, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 15570, cbuffer, 1195, 1216, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 15822, cbuffer, 1216, 1279, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 17712, 15570, 15822, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 190, 193, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9, cbuffer, 193, 202, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 36, cbuffer, 202, 220, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 90, cbuffer, 220, 250, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 189, cbuffer, 0, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 243, cbuffer, 3, 9, 36, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 324, cbuffer, 12, 36, 90, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 504, 180, 189, 243, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 558, 216, 243, 324, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 720, 486, 504, 558, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 810, cbuffer, 295, 301, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 828, cbuffer, 301, 319, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 882, cbuffer, 319, 355, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 990, cbuffer, 355, 415, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1188, cbuffer, 30, 810, 828, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1296, cbuffer, 36, 828, 882, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1458, cbuffer, 54, 882, 990, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1818, 1170, 1188, 1296, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1926, 1242, 1296, 1458, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 2250, 1782, 1818, 1926, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2430, cbuffer, 505, 515, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2460, cbuffer, 515, 545, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2550, cbuffer, 545, 605, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 2730, cbuffer, 605, 705, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 3060, cbuffer, 90, 2430, 2460, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3240, cbuffer, 100, 2460, 2550, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 3510, cbuffer, 130, 2550, 2730, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 4110, 3030, 3060, 3240, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 4290, 3150, 3240, 3510, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 4830, 4050, 4110, 4290, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 5130, cbuffer, 1405, 1408, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5139, cbuffer, 1408, 1417, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5166, cbuffer, 1417, 1435, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 5220, cbuffer, 1435, 1465, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 5319, cbuffer, 855, 5130, 5139, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 5373, cbuffer, 858, 5139, 5166, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 5454, cbuffer, 867, 5166, 5220, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 5634, 5310, 5319, 5373, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 5688, 5346, 5373, 5454, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 5850, 5616, 5634, 5688, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 5940, cbuffer, 1510, 1516, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5958, cbuffer, 1516, 1534, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6012, cbuffer, 1534, 1570, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6120, cbuffer, 1570, 1630, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 6318, cbuffer, 885, 5940, 5958, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 6426, cbuffer, 891, 5958, 6012, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 6588, cbuffer, 909, 6012, 6120, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 6948, 6300, 6318, 6426, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 7056, 6372, 6426, 6588, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 7380, 6912, 6948, 7056, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 7560, cbuffer, 1720, 1730, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 7590, cbuffer, 1730, 1760, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 7680, cbuffer, 1760, 1820, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 7860, cbuffer, 1820, 1920, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 8190, cbuffer, 945, 7560, 7590, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 8370, cbuffer, 955, 7590, 7680, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 8640, cbuffer, 985, 7680, 7860, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 9240, 8160, 8190, 8370, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 9420, 8280, 8370, 8640, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 9960, 9180, 9240, 9420, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 10260, cbuffer, 2070, 2085, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 10305, cbuffer, 2085, 2130, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 10440, cbuffer, 2130, 2220, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 10710, cbuffer, 2220, 2370, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 11205, cbuffer, 1045, 10260, 10305, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 11475, cbuffer, 1060, 10305, 10440, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 11880, cbuffer, 1105, 10440, 10710, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 12780, 11160, 11205, 11475, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 13050, 11340, 11475, 11880, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 13860, 12690, 12780, 13050, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 14310, cbuffer, 2595, 2616, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 14373, cbuffer, 2616, 2679, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 14562, cbuffer, 2679, 2805, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 14940, cbuffer, 2805, 3015, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 15633, cbuffer, 1195, 14310, 14373, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 16011, cbuffer, 1216, 14373, 14562, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 16578, cbuffer, 1279, 14562, 14940, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 17838, 15570, 15633, 16011, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 18216, 15822, 16011, 16578, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 19350, 17712, 17838, 18216, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 0, ckbuffer, 720, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21, ckbuffer, 750, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 42, ckbuffer, 780, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 252, ckbuffer, 2250, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 294, ckbuffer, 2310, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 336, ckbuffer, 2370, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 756, ckbuffer, 4830, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 826, ckbuffer, 4930, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 896, ckbuffer, 5030, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12369, ckbuffer, 5850, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12390, ckbuffer, 5880, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12411, ckbuffer, 5910, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12432, ckbuffer, 7380, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12474, ckbuffer, 7440, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12516, ckbuffer, 7500, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12558, ckbuffer, 9960, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12628, ckbuffer, 10060, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12698, ckbuffer, 10160, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12768, ckbuffer, 13860, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12873, ckbuffer, 14010, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12978, ckbuffer, 14160, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 13083, ckbuffer, 19350, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 13230, ckbuffer, 19560, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 13377, ckbuffer, 19770, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 2541, 0, 252, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 2604, 21, 294, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 2667, 42, 336, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 3297, 252, 756, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 3423, 294, 826, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 3549, 336, 896, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 6699, 2541, 3297, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 6825, 2604, 3423, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 6951, 2667, 3549, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 63, 12369, 12432, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 378, 12432, 12558, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 966, 12558, 12768, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 1596, 12768, 13083, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 2730, 0, 63, 378, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 3675, 252, 378, 966, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 4809, 756, 966, 1596, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 7077, 2541, 2730, 3675, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 8211, 3297, 3675, 4809, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fpxx(skbuffer, 10479, 6699, 7077, 8211, r_ab, 3, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 10479, 3, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 147, skbuffer, 10689, 3, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 294, skbuffer, 10899, 3, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 441, skbuffer, 11109, 3, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 588, skbuffer, 11319, 3, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 735, skbuffer, 11529, 3, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 882, skbuffer, 11739, 3, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1029, skbuffer, 11949, 3, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1176, skbuffer, 12159, 3, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 3, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFPFS_hpp */
