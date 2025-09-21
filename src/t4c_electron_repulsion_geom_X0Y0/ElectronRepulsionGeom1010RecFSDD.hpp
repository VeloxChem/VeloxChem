#ifndef ElectronRepulsionGeom1010RecFSDD_hpp
#define ElectronRepulsionGeom1010RecFSDD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"
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
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FS|1/|r-r'||DD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fsdd(T& distributor,
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

    CSimdArray<double> pbuffer(3835, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(3060, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(16335, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(23625, 1);

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

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 245, 85, 95, 155, 170, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 266, 95, 105, 170, 185, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 287, 105, 115, 185, 200, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 308, 115, 125, 200, 215, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 329, 125, 135, 215, 230, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 350, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 353, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 356, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 365, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 374, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 383, 13, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 401, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 419, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 437, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 455, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 485, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 515, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 545, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 575, 95, 155, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 620, 105, 170, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 665, 115, 185, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 710, 125, 200, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 755, 170, 245, 266, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 818, 185, 266, 287, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 881, 200, 287, 308, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 944, 215, 308, 329, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1007, 2, 3, 350, 353, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1013, 13, 16, 350, 356, 365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1031, 16, 19, 353, 365, 374, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1049, 37, 43, 356, 383, 401, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1085, 43, 49, 365, 401, 419, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1121, 49, 55, 374, 419, 437, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1157, 85, 95, 401, 455, 485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1217, 95, 105, 419, 485, 515, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1277, 105, 115, 437, 515, 545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1337, 155, 170, 485, 575, 620, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1427, 170, 185, 515, 620, 665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1517, 185, 200, 545, 665, 710, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1607, 245, 266, 620, 755, 818, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1733, 266, 287, 665, 818, 881, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1859, 287, 308, 710, 881, 944, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1985, 356, 365, 1007, 1013, 1031, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2015, 383, 401, 1013, 1049, 1085, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2075, 401, 419, 1031, 1085, 1121, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2135, 455, 485, 1085, 1157, 1217, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2235, 485, 515, 1121, 1217, 1277, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2335, 575, 620, 1217, 1337, 1427, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2485, 620, 665, 1277, 1427, 1517, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 2635, 755, 818, 1427, 1607, 1733, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 2845, 818, 881, 1517, 1733, 1859, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3055, 1049, 1085, 1985, 2015, 2075, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3145, 1157, 1217, 2075, 2135, 2235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3295, 1337, 1427, 2235, 2335, 2485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 3520, 1607, 1733, 2485, 2635, 2845, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 383, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 34, pbuffer, 455, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 64, pbuffer, 1049, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 1157, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {37, 43});

                pbuffer.scale(2.0 * a_exp, {85, 95});

                pbuffer.scale(2.0 * a_exp, {383, 401});

                pbuffer.scale(2.0 * a_exp, {455, 485});

                pbuffer.scale(2.0 * a_exp, {1049, 1085});

                pbuffer.scale(2.0 * a_exp, {1157, 1217});

                pbuffer.scale(2.0 * a_exp, {2015, 2075});

                pbuffer.scale(2.0 * a_exp, {2135, 2235});

                pbuffer.scale(2.0 * a_exp, {3055, 3145});

                pbuffer.scale(2.0 * a_exp, {3145, 3295});

                t2cfunc::reduce(cbuffer, 680, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 686, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 696, pbuffer, 383, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 714, pbuffer, 455, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 744, pbuffer, 1049, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 780, pbuffer, 1157, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 840, pbuffer, 2015, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 900, pbuffer, 2135, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1000, pbuffer, 3055, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1090, pbuffer, 3145, 150, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {37, 43});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {85, 95});

                pbuffer.scale(pfactors, 0, 2.0, {155, 170});

                pbuffer.scale(pfactors, 0, 2.0, {245, 266});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {383, 401});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {455, 485});

                pbuffer.scale(pfactors, 0, 2.0, {575, 620});

                pbuffer.scale(pfactors, 0, 2.0, {755, 818});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1049, 1085});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1157, 1217});

                pbuffer.scale(pfactors, 0, 2.0, {1337, 1427});

                pbuffer.scale(pfactors, 0, 2.0, {1607, 1733});

                t2cfunc::reduce(cbuffer, 160, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 166, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 176, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 191, pbuffer, 245, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 212, pbuffer, 383, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 230, pbuffer, 455, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 260, pbuffer, 575, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 305, pbuffer, 755, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 368, pbuffer, 1049, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 404, pbuffer, 1157, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 464, pbuffer, 1337, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 554, pbuffer, 1607, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {37, 43});

                pbuffer.scale(2.0 * a_exp, {85, 95});

                pbuffer.scale(2.0 * a_exp, {155, 170});

                pbuffer.scale(2.0 * a_exp, {245, 266});

                pbuffer.scale(2.0 * a_exp, {383, 401});

                pbuffer.scale(2.0 * a_exp, {455, 485});

                pbuffer.scale(2.0 * a_exp, {575, 620});

                pbuffer.scale(2.0 * a_exp, {755, 818});

                pbuffer.scale(2.0 * a_exp, {1049, 1085});

                pbuffer.scale(2.0 * a_exp, {1157, 1217});

                pbuffer.scale(2.0 * a_exp, {1337, 1427});

                pbuffer.scale(2.0 * a_exp, {1607, 1733});

                pbuffer.scale(pfactors, 0, 2.0, {2015, 2075});

                pbuffer.scale(pfactors, 0, 2.0, {2135, 2235});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2335, 2485});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2635, 2845});

                pbuffer.scale(pfactors, 0, 2.0, {3055, 3145});

                pbuffer.scale(pfactors, 0, 2.0, {3145, 3295});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3295, 3520});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3520, 3835});

                t2cfunc::reduce(cbuffer, 1240, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1246, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1256, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1271, pbuffer, 245, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1292, pbuffer, 383, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1310, pbuffer, 455, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1340, pbuffer, 575, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1385, pbuffer, 755, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1448, pbuffer, 1049, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1484, pbuffer, 1157, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1544, pbuffer, 1337, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1634, pbuffer, 1607, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1760, pbuffer, 2015, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1820, pbuffer, 2135, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1920, pbuffer, 2335, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2070, pbuffer, 2635, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2280, pbuffer, 3055, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2370, pbuffer, 3145, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2520, pbuffer, 3295, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2745, pbuffer, 3520, 315, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 93, cbuffer, 0, 6, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 642, cbuffer, 16, 34, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2010, cbuffer, 64, 100, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 3723, cbuffer, 680, 686, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4272, cbuffer, 696, 714, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5640, cbuffer, 744, 780, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8190, cbuffer, 840, 900, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 12285, cbuffer, 1000, 1090, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 0, cbuffer, 160, 166, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 18, cbuffer, 166, 176, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 48, cbuffer, 176, 191, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 111, cbuffer, 0, 0, 18, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 165, cbuffer, 6, 18, 48, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 255, 93, 111, 165, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 363, cbuffer, 212, 230, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 417, cbuffer, 230, 260, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 507, cbuffer, 260, 305, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 696, cbuffer, 16, 363, 417, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 858, cbuffer, 34, 417, 507, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 1128, 642, 696, 858, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1452, cbuffer, 368, 404, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1560, cbuffer, 404, 464, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 1740, cbuffer, 464, 554, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2118, cbuffer, 64, 1452, 1560, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2442, cbuffer, 100, 1560, 1740, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 2982, 2010, 2118, 2442, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3630, cbuffer, 1240, 1246, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3648, cbuffer, 1246, 1256, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 3678, cbuffer, 1256, 1271, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 3741, cbuffer, 680, 3630, 3648, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 3795, cbuffer, 686, 3648, 3678, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 3885, 3723, 3741, 3795, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3993, cbuffer, 1292, 1310, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 4047, cbuffer, 1310, 1340, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 4137, cbuffer, 1340, 1385, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4326, cbuffer, 696, 3993, 4047, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 4488, cbuffer, 714, 4047, 4137, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 4758, 4272, 4326, 4488, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5082, cbuffer, 1448, 1484, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 5190, cbuffer, 1484, 1544, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 5370, cbuffer, 1544, 1634, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 5748, cbuffer, 744, 5082, 5190, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 6072, cbuffer, 780, 5190, 5370, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 6612, 5640, 5748, 6072, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 7260, cbuffer, 1760, 1820, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 7440, cbuffer, 1820, 1920, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 7740, cbuffer, 1920, 2070, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 8370, cbuffer, 840, 7260, 7440, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 8910, cbuffer, 900, 7440, 7740, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 9810, 8190, 8370, 8910, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 10890, cbuffer, 2280, 2370, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 11160, cbuffer, 2370, 2520, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 11610, cbuffer, 2520, 2745, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 12555, cbuffer, 1000, 10890, 11160, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 13365, cbuffer, 1090, 11160, 11610, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 14715, 12285, 12555, 13365, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 0, ckbuffer, 255, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 25, ckbuffer, 291, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 50, ckbuffer, 327, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 300, ckbuffer, 1128, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 375, ckbuffer, 1236, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 450, ckbuffer, 1344, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 1200, ckbuffer, 2982, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 1350, ckbuffer, 3198, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 1500, ckbuffer, 3414, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21000, ckbuffer, 3885, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21025, ckbuffer, 3921, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21050, ckbuffer, 3957, 0, 0);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21075, ckbuffer, 4758, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21150, ckbuffer, 4866, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21225, ckbuffer, 4974, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21300, ckbuffer, 6612, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21450, ckbuffer, 6828, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21600, ckbuffer, 7044, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 21750, ckbuffer, 9810, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 22000, ckbuffer, 10170, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 22250, ckbuffer, 10530, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 22500, ckbuffer, 14715, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 22875, ckbuffer, 15255, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 23250, ckbuffer, 15795, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 5250, 0, 300, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 5325, 25, 375, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 5400, 50, 450, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 6150, 300, 1200, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 6375, 375, 1350, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 6600, 450, 1500, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 12900, 5250, 6150, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 13050, 5325, 6375, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 13200, 5400, 6600, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 75, 21000, 21075, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 525, 21075, 21300, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 1650, 21300, 21750, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 3000, 21750, 22500, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 5475, 0, 75, 525, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 6825, 300, 525, 1650, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 8850, 1200, 1650, 3000, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 13350, 5250, 5475, 6825, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 14700, 6150, 6825, 8850, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fsxx(skbuffer, 18750, 12900, 13350, 14700, r_ab, 2, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 18750, 2, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 175, skbuffer, 19000, 2, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 350, skbuffer, 19250, 2, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 525, skbuffer, 19500, 2, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 700, skbuffer, 19750, 2, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 875, skbuffer, 20000, 2, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1050, skbuffer, 20250, 2, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1225, skbuffer, 20500, 2, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 1400, skbuffer, 20750, 2, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 2, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFSDD_hpp */
