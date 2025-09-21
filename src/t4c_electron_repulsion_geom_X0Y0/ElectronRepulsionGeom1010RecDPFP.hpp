#ifndef ElectronRepulsionGeom1010RecDPFP_hpp
#define ElectronRepulsionGeom1010RecDPFP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXDP.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXFP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSH.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DP|1/|r-r'||FP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dpfp(T& distributor,
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

    CSimdArray<double> pbuffer(3956, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(3182, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(24510, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(15372, 1);

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

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 245, 85, 95, 155, 170, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 266, 95, 105, 170, 185, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 287, 105, 115, 185, 200, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 308, 115, 125, 200, 215, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 329, 125, 135, 215, 230, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 350, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 353, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 356, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 359, 1, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 368, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 377, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 386, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 395, 13, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 413, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 431, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 449, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 467, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 497, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 527, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 557, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 587, 95, 155, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 632, 105, 170, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 677, 115, 185, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 722, 125, 200, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 767, 170, 245, 266, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 830, 185, 266, 287, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 893, 200, 287, 308, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 956, 215, 308, 329, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1019, 1, 2, 350, 353, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1025, 2, 3, 353, 356, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1031, 10, 13, 350, 359, 368, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1049, 13, 16, 353, 368, 377, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1067, 16, 19, 356, 377, 386, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1085, 37, 43, 368, 395, 413, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1121, 43, 49, 377, 413, 431, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1157, 49, 55, 386, 431, 449, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1193, 85, 95, 413, 467, 497, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1253, 95, 105, 431, 497, 527, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1313, 105, 115, 449, 527, 557, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1373, 155, 170, 497, 587, 632, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1463, 170, 185, 527, 632, 677, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1553, 185, 200, 557, 677, 722, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1643, 245, 266, 632, 767, 830, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1769, 266, 287, 677, 830, 893, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1895, 287, 308, 722, 893, 956, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2021, 350, 353, 1019, 1025, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2031, 359, 368, 1019, 1031, 1049, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2061, 368, 377, 1025, 1049, 1067, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2091, 395, 413, 1049, 1085, 1121, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2151, 413, 431, 1067, 1121, 1157, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2211, 467, 497, 1121, 1193, 1253, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2311, 497, 527, 1157, 1253, 1313, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2411, 587, 632, 1253, 1373, 1463, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2561, 632, 677, 1313, 1463, 1553, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 2711, 767, 830, 1463, 1643, 1769, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 2921, 830, 893, 1553, 1769, 1895, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3131, 1031, 1049, 2021, 2031, 2061, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3176, 1085, 1121, 2061, 2091, 2151, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3266, 1193, 1253, 2151, 2211, 2311, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3416, 1373, 1463, 2311, 2411, 2561, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 3641, 1643, 1769, 2561, 2711, 2921, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 359, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 395, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 27, pbuffer, 467, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 57, pbuffer, 1031, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 75, pbuffer, 1085, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 111, pbuffer, 1193, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {359, 368});

                pbuffer.scale(2.0 * a_exp, {395, 413});

                pbuffer.scale(2.0 * a_exp, {467, 497});

                pbuffer.scale(2.0 * a_exp, {1031, 1049});

                pbuffer.scale(2.0 * a_exp, {1085, 1121});

                pbuffer.scale(2.0 * a_exp, {1193, 1253});

                pbuffer.scale(2.0 * a_exp, {2031, 2061});

                pbuffer.scale(2.0 * a_exp, {2091, 2151});

                pbuffer.scale(2.0 * a_exp, {2211, 2311});

                pbuffer.scale(2.0 * a_exp, {3131, 3176});

                pbuffer.scale(2.0 * a_exp, {3176, 3266});

                pbuffer.scale(2.0 * a_exp, {3266, 3416});

                t2cfunc::reduce(cbuffer, 666, pbuffer, 359, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 675, pbuffer, 395, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 693, pbuffer, 467, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 723, pbuffer, 1031, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 741, pbuffer, 1085, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 777, pbuffer, 1193, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 837, pbuffer, 2031, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 867, pbuffer, 2091, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 927, pbuffer, 2211, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1027, pbuffer, 3131, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1072, pbuffer, 3176, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1162, pbuffer, 3266, 150, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {359, 368});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {395, 413});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {467, 497});

                pbuffer.scale(pfactors, 0, 2.0, {587, 632});

                pbuffer.scale(pfactors, 0, 2.0, {767, 830});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1031, 1049});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1085, 1121});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1193, 1253});

                pbuffer.scale(pfactors, 0, 2.0, {1373, 1463});

                pbuffer.scale(pfactors, 0, 2.0, {1643, 1769});

                t2cfunc::reduce(cbuffer, 171, pbuffer, 359, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 180, pbuffer, 395, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 198, pbuffer, 467, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 228, pbuffer, 587, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 273, pbuffer, 767, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 336, pbuffer, 1031, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 354, pbuffer, 1085, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 390, pbuffer, 1193, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 450, pbuffer, 1373, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 540, pbuffer, 1643, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {359, 368});

                pbuffer.scale(2.0 * a_exp, {395, 413});

                pbuffer.scale(2.0 * a_exp, {467, 497});

                pbuffer.scale(2.0 * a_exp, {587, 632});

                pbuffer.scale(2.0 * a_exp, {767, 830});

                pbuffer.scale(2.0 * a_exp, {1031, 1049});

                pbuffer.scale(2.0 * a_exp, {1085, 1121});

                pbuffer.scale(2.0 * a_exp, {1193, 1253});

                pbuffer.scale(2.0 * a_exp, {1373, 1463});

                pbuffer.scale(2.0 * a_exp, {1643, 1769});

                pbuffer.scale(pfactors, 0, 2.0, {2031, 2061});

                pbuffer.scale(pfactors, 0, 2.0, {2091, 2151});

                pbuffer.scale(pfactors, 0, 2.0, {2211, 2311});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2411, 2561});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2711, 2921});

                pbuffer.scale(pfactors, 0, 2.0, {3131, 3176});

                pbuffer.scale(pfactors, 0, 2.0, {3176, 3266});

                pbuffer.scale(pfactors, 0, 2.0, {3266, 3416});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3416, 3641});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3641, 3956});

                t2cfunc::reduce(cbuffer, 1312, pbuffer, 359, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1321, pbuffer, 395, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1339, pbuffer, 467, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1369, pbuffer, 587, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1414, pbuffer, 767, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1477, pbuffer, 1031, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1495, pbuffer, 1085, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1531, pbuffer, 1193, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1591, pbuffer, 1373, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1681, pbuffer, 1643, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1807, pbuffer, 2031, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1837, pbuffer, 2091, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1897, pbuffer, 2211, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1997, pbuffer, 2411, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2147, pbuffer, 2711, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2357, pbuffer, 3131, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2402, pbuffer, 3176, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2492, pbuffer, 3266, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2642, pbuffer, 3416, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2867, pbuffer, 3641, 315, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 306, cbuffer, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 414, cbuffer, 9, 27, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 900, 306, 414, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2322, cbuffer, 57, 75, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2538, cbuffer, 75, 111, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 3510, 2322, 2538, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 5436, cbuffer, 666, 675, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5544, cbuffer, 675, 693, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 6030, 5436, 5544, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 7452, cbuffer, 723, 741, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 7668, cbuffer, 741, 777, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 8640, 7452, 7668, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 11280, cbuffer, 837, 867, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 11640, cbuffer, 867, 927, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 13260, 11280, 11640, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 17490, cbuffer, 1027, 1072, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 18030, cbuffer, 1072, 1162, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 20460, 17490, 18030, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 171, 180, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 27, cbuffer, 180, 198, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 81, cbuffer, 198, 228, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 171, cbuffer, 228, 273, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 333, cbuffer, 0, 0, 27, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 468, cbuffer, 9, 27, 81, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 630, cbuffer, 27, 81, 171, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 954, 306, 333, 468, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 1116, 414, 468, 630, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 1440, 900, 954, 1116, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1710, cbuffer, 336, 354, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1764, cbuffer, 354, 390, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1872, cbuffer, 390, 450, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 2052, cbuffer, 450, 540, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2376, cbuffer, 57, 1710, 1764, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2646, cbuffer, 75, 1764, 1872, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2970, cbuffer, 111, 1872, 2052, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 3618, 2322, 2376, 2646, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 3942, 2538, 2646, 2970, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 4590, 3510, 3618, 3942, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5130, cbuffer, 1312, 1321, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5157, cbuffer, 1321, 1339, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 5211, cbuffer, 1339, 1369, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 5301, cbuffer, 1369, 1414, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 5463, cbuffer, 666, 5130, 5157, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 5598, cbuffer, 675, 5157, 5211, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 5760, cbuffer, 693, 5211, 5301, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 6084, 5436, 5463, 5598, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 6246, 5544, 5598, 5760, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 6570, 6030, 6084, 6246, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 6840, cbuffer, 1477, 1495, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6894, cbuffer, 1495, 1531, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 7002, cbuffer, 1531, 1591, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 7182, cbuffer, 1591, 1681, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 7506, cbuffer, 723, 6840, 6894, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 7776, cbuffer, 741, 6894, 7002, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 8100, cbuffer, 777, 7002, 7182, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 8748, 7452, 7506, 7776, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 9072, 7668, 7776, 8100, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 9720, 8640, 8748, 9072, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 10260, cbuffer, 1807, 1837, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 10350, cbuffer, 1837, 1897, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 10530, cbuffer, 1897, 1997, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 10830, cbuffer, 1997, 2147, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 11370, cbuffer, 837, 10260, 10350, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 11820, cbuffer, 867, 10350, 10530, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 12360, cbuffer, 927, 10530, 10830, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 13440, 11280, 11370, 11820, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 13980, 11640, 11820, 12360, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 15060, 13260, 13440, 13980, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 15960, cbuffer, 2357, 2402, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 16095, cbuffer, 2402, 2492, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 16365, cbuffer, 2492, 2642, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 16815, cbuffer, 2642, 2867, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 17625, cbuffer, 1027, 15960, 16095, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 18300, cbuffer, 1072, 16095, 16365, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 19110, cbuffer, 1162, 16365, 16815, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 20730, 17490, 17625, 18300, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 21540, 18030, 18300, 19110, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 23160, 20460, 20730, 21540, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 0, ckbuffer, 1440, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 63, ckbuffer, 1530, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 126, ckbuffer, 1620, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 756, ckbuffer, 4590, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 882, ckbuffer, 4770, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 1008, ckbuffer, 4950, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 13230, ckbuffer, 6570, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 13293, ckbuffer, 6660, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 13356, ckbuffer, 6750, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 13419, ckbuffer, 9720, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 13545, ckbuffer, 9900, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 13671, ckbuffer, 10080, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 13797, ckbuffer, 15060, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 14007, ckbuffer, 15360, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 14217, ckbuffer, 15660, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 14427, ckbuffer, 23160, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 14742, ckbuffer, 23610, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 15057, ckbuffer, 24060, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 4158, 0, 756, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 4347, 63, 882, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 4536, 126, 1008, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 189, 13230, 13419, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 1134, 13419, 13797, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 2268, 13797, 14427, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 4725, 0, 189, 1134, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 6426, 756, 1134, 2268, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 9828, 4158, 4725, 6426, r_ab, 3, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 9828, 3, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 315, skbuffer, 10206, 3, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 630, skbuffer, 10584, 3, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 945, skbuffer, 10962, 3, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1260, skbuffer, 11340, 3, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1575, skbuffer, 11718, 3, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1890, skbuffer, 12096, 3, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2205, skbuffer, 12474, 3, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 2520, skbuffer, 12852, 3, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 3, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDPFP_hpp */
