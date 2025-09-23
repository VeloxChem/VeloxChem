#ifndef ElectronRepulsionGeom1010RecFSFP_hpp
#define ElectronRepulsionGeom1010RecFSFP_hpp

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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FS|1/|r-r'||FP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fsfp(T& distributor,
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

    CSimdArray<double> cbuffer(3330, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(25650, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(19845, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 10, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 19, pbuffer, 359, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 28, pbuffer, 395, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 46, pbuffer, 467, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 76, pbuffer, 1031, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 94, pbuffer, 1085, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 130, pbuffer, 1193, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {10, 13});

                pbuffer.scale(2.0 * a_exp, {37, 43});

                pbuffer.scale(2.0 * a_exp, {85, 95});

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

                t2cfunc::reduce(cbuffer, 740, pbuffer, 10, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 743, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 749, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 759, pbuffer, 359, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 768, pbuffer, 395, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 786, pbuffer, 467, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 816, pbuffer, 1031, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 834, pbuffer, 1085, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 870, pbuffer, 1193, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 930, pbuffer, 2031, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 960, pbuffer, 2091, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1020, pbuffer, 2211, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1120, pbuffer, 3131, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1165, pbuffer, 3176, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1255, pbuffer, 3266, 150, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {10, 13});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {37, 43});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {85, 95});

                pbuffer.scale(pfactors, 0, 2.0, {155, 170});

                pbuffer.scale(pfactors, 0, 2.0, {245, 266});

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

                t2cfunc::reduce(cbuffer, 190, pbuffer, 10, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 193, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 199, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 209, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 224, pbuffer, 245, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 245, pbuffer, 359, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 254, pbuffer, 395, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 272, pbuffer, 467, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 302, pbuffer, 587, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 347, pbuffer, 767, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 410, pbuffer, 1031, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 428, pbuffer, 1085, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 464, pbuffer, 1193, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 524, pbuffer, 1373, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 614, pbuffer, 1643, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {10, 13});

                pbuffer.scale(2.0 * a_exp, {37, 43});

                pbuffer.scale(2.0 * a_exp, {85, 95});

                pbuffer.scale(2.0 * a_exp, {155, 170});

                pbuffer.scale(2.0 * a_exp, {245, 266});

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

                t2cfunc::reduce(cbuffer, 1405, pbuffer, 10, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1408, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1414, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1424, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1439, pbuffer, 245, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1460, pbuffer, 359, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1469, pbuffer, 395, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1487, pbuffer, 467, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1517, pbuffer, 587, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1562, pbuffer, 767, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1625, pbuffer, 1031, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1643, pbuffer, 1085, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1679, pbuffer, 1193, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1739, pbuffer, 1373, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1829, pbuffer, 1643, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1955, pbuffer, 2031, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1985, pbuffer, 2091, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2045, pbuffer, 2211, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2145, pbuffer, 2411, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2295, pbuffer, 2711, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2505, pbuffer, 3131, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2550, pbuffer, 3176, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2640, pbuffer, 3266, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2790, pbuffer, 3416, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3015, pbuffer, 3641, 315, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 102, cbuffer, 0, 3, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 138, cbuffer, 3, 9, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 300, 102, 138, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 876, cbuffer, 19, 28, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 984, cbuffer, 28, 46, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 1470, 876, 984, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2892, cbuffer, 76, 94, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3108, cbuffer, 94, 130, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 4080, 2892, 3108, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 5802, cbuffer, 740, 743, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5838, cbuffer, 743, 749, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 6000, 5802, 5838, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6576, cbuffer, 759, 768, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 6684, cbuffer, 768, 786, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 7170, 6576, 6684, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 8592, cbuffer, 816, 834, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8808, cbuffer, 834, 870, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 9780, 8592, 8808, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 12420, cbuffer, 930, 960, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 12780, cbuffer, 960, 1020, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 14400, 12420, 12780, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 18630, cbuffer, 1120, 1165, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 19170, cbuffer, 1165, 1255, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 21600, 18630, 19170, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 190, 193, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9, cbuffer, 193, 199, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 27, cbuffer, 199, 209, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 57, cbuffer, 209, 224, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 111, cbuffer, 0, 0, 9, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 156, cbuffer, 3, 9, 27, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 210, cbuffer, 9, 27, 57, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 318, 102, 111, 156, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 372, 138, 156, 210, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 480, 300, 318, 372, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 570, cbuffer, 245, 254, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 597, cbuffer, 254, 272, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 651, cbuffer, 272, 302, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 741, cbuffer, 302, 347, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 903, cbuffer, 19, 570, 597, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1038, cbuffer, 28, 597, 651, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 1200, cbuffer, 46, 651, 741, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1524, 876, 903, 1038, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 1686, 984, 1038, 1200, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 2010, 1470, 1524, 1686, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2280, cbuffer, 410, 428, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2334, cbuffer, 428, 464, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 2442, cbuffer, 464, 524, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 2622, cbuffer, 524, 614, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2946, cbuffer, 76, 2280, 2334, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 3216, cbuffer, 94, 2334, 2442, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 3540, cbuffer, 130, 2442, 2622, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 4188, 2892, 2946, 3216, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 4512, 3108, 3216, 3540, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 5160, 4080, 4188, 4512, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5700, cbuffer, 1405, 1408, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5709, cbuffer, 1408, 1414, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 5727, cbuffer, 1414, 1424, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 5757, cbuffer, 1424, 1439, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 5811, cbuffer, 740, 5700, 5709, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 5856, cbuffer, 743, 5709, 5727, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 5910, cbuffer, 749, 5727, 5757, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 6018, 5802, 5811, 5856, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 6072, 5838, 5856, 5910, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 6180, 6000, 6018, 6072, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 6270, cbuffer, 1460, 1469, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6297, cbuffer, 1469, 1487, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6351, cbuffer, 1487, 1517, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 6441, cbuffer, 1517, 1562, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 6603, cbuffer, 759, 6270, 6297, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 6738, cbuffer, 768, 6297, 6351, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 6900, cbuffer, 786, 6351, 6441, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 7224, 6576, 6603, 6738, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 7386, 6684, 6738, 6900, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 7710, 7170, 7224, 7386, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 7980, cbuffer, 1625, 1643, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 8034, cbuffer, 1643, 1679, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 8142, cbuffer, 1679, 1739, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 8322, cbuffer, 1739, 1829, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 8646, cbuffer, 816, 7980, 8034, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 8916, cbuffer, 834, 8034, 8142, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 9240, cbuffer, 870, 8142, 8322, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 9888, 8592, 8646, 8916, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 10212, 8808, 8916, 9240, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 10860, 9780, 9888, 10212, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 11400, cbuffer, 1955, 1985, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 11490, cbuffer, 1985, 2045, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 11670, cbuffer, 2045, 2145, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 11970, cbuffer, 2145, 2295, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 12510, cbuffer, 930, 11400, 11490, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 12960, cbuffer, 960, 11490, 11670, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 13500, cbuffer, 1020, 11670, 11970, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 14580, 12420, 12510, 12960, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 15120, 12780, 12960, 13500, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 16200, 14400, 14580, 15120, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 17100, cbuffer, 2505, 2550, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 17235, cbuffer, 2550, 2640, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 17505, cbuffer, 2640, 2790, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 17955, cbuffer, 2790, 3015, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 18765, cbuffer, 1120, 17100, 17235, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 19440, cbuffer, 1165, 17235, 17505, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 20250, cbuffer, 1255, 17505, 17955, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 21870, 18630, 18765, 19440, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 22680, 19170, 19440, 20250, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 24300, 21600, 21870, 22680, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 0, ckbuffer, 480, 0, 0);

            t4cfunc::ket_transform<3, 1>(skbuffer, 21, ckbuffer, 510, 0, 0);

            t4cfunc::ket_transform<3, 1>(skbuffer, 42, ckbuffer, 540, 0, 0);

            t4cfunc::ket_transform<3, 1>(skbuffer, 252, ckbuffer, 2010, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 315, ckbuffer, 2100, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 378, ckbuffer, 2190, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 1008, ckbuffer, 5160, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 1134, ckbuffer, 5340, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 1260, ckbuffer, 5520, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 17640, ckbuffer, 6180, 0, 0);

            t4cfunc::ket_transform<3, 1>(skbuffer, 17661, ckbuffer, 6210, 0, 0);

            t4cfunc::ket_transform<3, 1>(skbuffer, 17682, ckbuffer, 6240, 0, 0);

            t4cfunc::ket_transform<3, 1>(skbuffer, 17703, ckbuffer, 7710, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 17766, ckbuffer, 7800, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 17829, ckbuffer, 7890, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 17892, ckbuffer, 10860, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 18018, ckbuffer, 11040, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 18144, ckbuffer, 11220, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 18270, ckbuffer, 16200, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 18480, ckbuffer, 16500, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 18690, ckbuffer, 16800, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 18900, ckbuffer, 24300, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 19215, ckbuffer, 24750, 0, 4);

            t4cfunc::ket_transform<3, 1>(skbuffer, 19530, ckbuffer, 25200, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 4410, 0, 252, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 4473, 21, 315, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 4536, 42, 378, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 5166, 252, 1008, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 5355, 315, 1134, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 5544, 378, 1260, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 10836, 4410, 5166, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 10962, 4473, 5355, r_ab, 3, 1);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 11088, 4536, 5544, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 63, 17640, 17703, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 441, 17703, 17892, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 1386, 17892, 18270, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 2520, 18270, 18900, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 4599, 0, 63, 441, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 5733, 252, 441, 1386, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 7434, 1008, 1386, 2520, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 11214, 4410, 4599, 5733, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 12348, 5166, 5733, 7434, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fsxx(skbuffer, 15750, 10836, 11214, 12348, r_ab, 3, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 15750, 3, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 147, skbuffer, 15960, 3, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 294, skbuffer, 16170, 3, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 441, skbuffer, 16380, 3, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 588, skbuffer, 16590, 3, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 735, skbuffer, 16800, 3, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 882, skbuffer, 17010, 3, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1029, skbuffer, 17220, 3, 1);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1176, skbuffer, 17430, 3, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 3, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFSFP_hpp */
