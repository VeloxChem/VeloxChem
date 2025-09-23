#ifndef ElectronRepulsionGeom1010RecDFFS_hpp
#define ElectronRepulsionGeom1010RecDFFS_hpp

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
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
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
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DF|1/|r-r'||FS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dffs(T& distributor,
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

    CSimdArray<double> pbuffer(7385, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(4455, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(26730, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(14112, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2205, 1);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 11, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 11);
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

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 10, pfactors, 16, bf_data, 10);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 29, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 32, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 35, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 38, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 41, 0, 1, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 47, 1, 2, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 2, 3, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 3, 4, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 65, 4, 5, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 71, 5, 6, 26, 29, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 77, 6, 7, 29, 32, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 83, 7, 8, 32, 35, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 89, 8, 9, 35, 38, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 11, 14, 41, 47, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 14, 17, 47, 53, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 17, 20, 53, 59, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 20, 23, 59, 65, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 23, 26, 65, 71, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 26, 29, 71, 77, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 29, 32, 77, 83, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 32, 35, 83, 89, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 175, 41, 47, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 190, 47, 53, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 205, 53, 59, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 220, 59, 65, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 235, 65, 71, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 250, 71, 77, 145, 155, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 265, 77, 83, 155, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 280, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 283, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 286, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 289, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 292, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 295, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 298, 1, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 307, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 316, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 325, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 334, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 343, 6, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 352, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 370, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 388, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 406, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 424, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 442, 29, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 460, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 490, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 520, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 550, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 580, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 610, 77, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 640, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 685, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 730, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 775, 135, 220, 235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 820, 145, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 865, 155, 250, 265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 910, 0, 1, 280, 283, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 916, 1, 2, 283, 286, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 922, 2, 3, 286, 289, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 928, 3, 4, 289, 292, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 934, 4, 5, 292, 295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 940, 11, 14, 283, 298, 307, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 958, 14, 17, 286, 307, 316, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 976, 17, 20, 289, 316, 325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 994, 20, 23, 292, 325, 334, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1012, 23, 26, 295, 334, 343, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1030, 41, 47, 307, 352, 370, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1066, 47, 53, 316, 370, 388, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1102, 53, 59, 325, 388, 406, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1138, 59, 65, 334, 406, 424, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1174, 65, 71, 343, 424, 442, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1210, 95, 105, 370, 460, 490, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1270, 105, 115, 388, 490, 520, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1330, 115, 125, 406, 520, 550, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1390, 125, 135, 424, 550, 580, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1450, 135, 145, 442, 580, 610, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1510, 175, 190, 490, 640, 685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1600, 190, 205, 520, 685, 730, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1690, 205, 220, 550, 730, 775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1780, 220, 235, 580, 775, 820, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1870, 235, 250, 610, 820, 865, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1960, 280, 283, 910, 916, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1970, 283, 286, 916, 922, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1980, 286, 289, 922, 928, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1990, 289, 292, 928, 934, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2000, 298, 307, 916, 940, 958, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2030, 307, 316, 922, 958, 976, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2060, 316, 325, 928, 976, 994, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2090, 325, 334, 934, 994, 1012, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2120, 352, 370, 958, 1030, 1066, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2180, 370, 388, 976, 1066, 1102, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2240, 388, 406, 994, 1102, 1138, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2300, 406, 424, 1012, 1138, 1174, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2360, 460, 490, 1066, 1210, 1270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2460, 490, 520, 1102, 1270, 1330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2560, 520, 550, 1138, 1330, 1390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2660, 550, 580, 1174, 1390, 1450, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2760, 640, 685, 1270, 1510, 1600, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2910, 685, 730, 1330, 1600, 1690, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3060, 730, 775, 1390, 1690, 1780, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3210, 775, 820, 1450, 1780, 1870, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 3360, 910, 916, 1960, 1970, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 3375, 916, 922, 1970, 1980, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 3390, 922, 928, 1980, 1990, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3405, 940, 958, 1970, 2000, 2030, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3450, 958, 976, 1980, 2030, 2060, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3495, 976, 994, 1990, 2060, 2090, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3540, 1030, 1066, 2030, 2120, 2180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3630, 1066, 1102, 2060, 2180, 2240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3720, 1102, 1138, 2090, 2240, 2300, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3810, 1210, 1270, 2180, 2360, 2460, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3960, 1270, 1330, 2240, 2460, 2560, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4110, 1330, 1390, 2300, 2560, 2660, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4260, 1510, 1600, 2460, 2760, 2910, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4485, 1600, 1690, 2560, 2910, 3060, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4710, 1690, 1780, 2660, 3060, 3210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 4935, 1960, 1970, 3360, 3375, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 4956, 1970, 1980, 3375, 3390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 4977, 2000, 2030, 3375, 3405, 3450, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 5040, 2030, 2060, 3390, 3450, 3495, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 5103, 2120, 2180, 3450, 3540, 3630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 5229, 2180, 2240, 3495, 3630, 3720, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5355, 2360, 2460, 3630, 3810, 3960, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5565, 2460, 2560, 3720, 3960, 4110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 5775, 2760, 2910, 3960, 4260, 4485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 6090, 2910, 3060, 4110, 4485, 4710, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 6405, 3360, 3375, 4935, 4956, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 6433, 3405, 3450, 4956, 4977, 5040, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 6517, 3540, 3630, 5040, 5103, 5229, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 6685, 3810, 3960, 5229, 5355, 5565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 6965, 4260, 4485, 5565, 5775, 6090, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1960, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 2000, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 40, pbuffer, 2120, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 3360, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 115, pbuffer, 3405, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 3540, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1960, 1970});

                pbuffer.scale(2.0 * a_exp, {2000, 2030});

                pbuffer.scale(2.0 * a_exp, {2120, 2180});

                pbuffer.scale(2.0 * a_exp, {3360, 3375});

                pbuffer.scale(2.0 * a_exp, {3405, 3450});

                pbuffer.scale(2.0 * a_exp, {3540, 3630});

                pbuffer.scale(2.0 * a_exp, {4935, 4956});

                pbuffer.scale(2.0 * a_exp, {4977, 5040});

                pbuffer.scale(2.0 * a_exp, {5103, 5229});

                pbuffer.scale(2.0 * a_exp, {6405, 6433});

                pbuffer.scale(2.0 * a_exp, {6433, 6517});

                pbuffer.scale(2.0 * a_exp, {6517, 6685});

                t2cfunc::reduce(cbuffer, 1125, pbuffer, 1960, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1135, pbuffer, 2000, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1165, pbuffer, 2120, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1225, pbuffer, 3360, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1240, pbuffer, 3405, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1285, pbuffer, 3540, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1375, pbuffer, 4935, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1396, pbuffer, 4977, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1459, pbuffer, 5103, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1585, pbuffer, 6405, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1613, pbuffer, 6433, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1697, pbuffer, 6517, 168, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1960, 1970});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2000, 2030});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2120, 2180});

                pbuffer.scale(pfactors, 0, 2.0, {2360, 2460});

                pbuffer.scale(pfactors, 0, 2.0, {2760, 2910});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3360, 3375});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3405, 3450});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3540, 3630});

                pbuffer.scale(pfactors, 0, 2.0, {3810, 3960});

                pbuffer.scale(pfactors, 0, 2.0, {4260, 4485});

                t2cfunc::reduce(cbuffer, 250, pbuffer, 1960, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 260, pbuffer, 2000, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 290, pbuffer, 2120, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 350, pbuffer, 2360, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 450, pbuffer, 2760, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 600, pbuffer, 3360, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 615, pbuffer, 3405, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 660, pbuffer, 3540, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 750, pbuffer, 3810, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 900, pbuffer, 4260, 225, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1960, 1970});

                pbuffer.scale(2.0 * a_exp, {2000, 2030});

                pbuffer.scale(2.0 * a_exp, {2120, 2180});

                pbuffer.scale(2.0 * a_exp, {2360, 2460});

                pbuffer.scale(2.0 * a_exp, {2760, 2910});

                pbuffer.scale(2.0 * a_exp, {3360, 3375});

                pbuffer.scale(2.0 * a_exp, {3405, 3450});

                pbuffer.scale(2.0 * a_exp, {3540, 3630});

                pbuffer.scale(2.0 * a_exp, {3810, 3960});

                pbuffer.scale(2.0 * a_exp, {4260, 4485});

                pbuffer.scale(pfactors, 0, 2.0, {4935, 4956});

                pbuffer.scale(pfactors, 0, 2.0, {4977, 5040});

                pbuffer.scale(pfactors, 0, 2.0, {5103, 5229});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5355, 5565});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5775, 6090});

                pbuffer.scale(pfactors, 0, 2.0, {6405, 6433});

                pbuffer.scale(pfactors, 0, 2.0, {6433, 6517});

                pbuffer.scale(pfactors, 0, 2.0, {6517, 6685});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {6685, 6965});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {6965, 7385});

                t2cfunc::reduce(cbuffer, 1865, pbuffer, 1960, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1875, pbuffer, 2000, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1905, pbuffer, 2120, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1965, pbuffer, 2360, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2065, pbuffer, 2760, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2215, pbuffer, 3360, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2230, pbuffer, 3405, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2275, pbuffer, 3540, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2365, pbuffer, 3810, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2515, pbuffer, 4260, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2740, pbuffer, 4935, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2761, pbuffer, 4977, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2824, pbuffer, 5103, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2950, pbuffer, 5355, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3160, pbuffer, 5775, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3475, pbuffer, 6405, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3503, pbuffer, 6433, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3587, pbuffer, 6517, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3755, pbuffer, 6685, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4035, pbuffer, 6965, 420, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 600, cbuffer, 0, 10, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 720, cbuffer, 10, 40, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 1620, 600, 720, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 3600, cbuffer, 100, 115, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3780, cbuffer, 115, 160, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 5130, 3600, 3780, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 7350, cbuffer, 1125, 1135, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 7470, cbuffer, 1135, 1165, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 8370, 7350, 7470, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 10350, cbuffer, 1225, 1240, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 10530, cbuffer, 1240, 1285, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 11880, 10350, 10530, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 14760, cbuffer, 1375, 1396, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 15012, cbuffer, 1396, 1459, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 16902, 14760, 15012, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 20850, cbuffer, 1585, 1613, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 21186, cbuffer, 1613, 1697, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 23706, 20850, 21186, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 250, 260, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 30, cbuffer, 260, 290, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 120, cbuffer, 290, 350, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 300, cbuffer, 350, 450, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 630, cbuffer, 0, 0, 30, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 810, cbuffer, 10, 30, 120, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1080, cbuffer, 40, 120, 300, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1680, 600, 630, 810, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1860, 720, 810, 1080, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 2400, 1620, 1680, 1860, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2700, cbuffer, 600, 615, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2745, cbuffer, 615, 660, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2880, cbuffer, 660, 750, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3150, cbuffer, 750, 900, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 3645, cbuffer, 100, 2700, 2745, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3915, cbuffer, 115, 2745, 2880, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4320, cbuffer, 160, 2880, 3150, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 5220, 3600, 3645, 3915, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 5490, 3780, 3915, 4320, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 6300, 5130, 5220, 5490, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 6750, cbuffer, 1865, 1875, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 6780, cbuffer, 1875, 1905, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6870, cbuffer, 1905, 1965, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 7050, cbuffer, 1965, 2065, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 7380, cbuffer, 1125, 6750, 6780, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 7560, cbuffer, 1135, 6780, 6870, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 7830, cbuffer, 1165, 6870, 7050, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 8430, 7350, 7380, 7560, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 8610, 7470, 7560, 7830, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 9150, 8370, 8430, 8610, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 9450, cbuffer, 2215, 2230, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9495, cbuffer, 2230, 2275, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9630, cbuffer, 2275, 2365, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 9900, cbuffer, 2365, 2515, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 10395, cbuffer, 1225, 9450, 9495, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 10665, cbuffer, 1240, 9495, 9630, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 11070, cbuffer, 1285, 9630, 9900, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 11970, 10350, 10395, 10665, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 12240, 10530, 10665, 11070, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 13050, 11880, 11970, 12240, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 13500, cbuffer, 2740, 2761, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 13563, cbuffer, 2761, 2824, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 13752, cbuffer, 2824, 2950, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 14130, cbuffer, 2950, 3160, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 14823, cbuffer, 1375, 13500, 13563, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 15201, cbuffer, 1396, 13563, 13752, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 15768, cbuffer, 1459, 13752, 14130, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 17028, 14760, 14823, 15201, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 17406, 15012, 15201, 15768, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 18540, 16902, 17028, 17406, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 19170, cbuffer, 3475, 3503, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 19254, cbuffer, 3503, 3587, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 19506, cbuffer, 3587, 3755, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 20010, cbuffer, 3755, 4035, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 20934, cbuffer, 1585, 19170, 19254, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 21438, cbuffer, 1613, 19254, 19506, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 22194, cbuffer, 1697, 19506, 20010, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 23874, 20850, 20934, 21438, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 24378, 21186, 21438, 22194, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 25890, 23706, 23874, 24378, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<3, 0>(skbuffer, 0, ckbuffer, 2400, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 70, ckbuffer, 2500, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 140, ckbuffer, 2600, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 840, ckbuffer, 6300, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 945, ckbuffer, 6450, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1050, ckbuffer, 6600, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12558, ckbuffer, 9150, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12628, ckbuffer, 9250, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12698, ckbuffer, 9350, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12768, ckbuffer, 13050, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12873, ckbuffer, 13200, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 12978, ckbuffer, 13350, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 13083, ckbuffer, 18540, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 13230, ckbuffer, 18750, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 13377, ckbuffer, 18960, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 13524, ckbuffer, 25890, 0, 6);

            t4cfunc::ket_transform<3, 0>(skbuffer, 13720, ckbuffer, 26170, 0, 6);

            t4cfunc::ket_transform<3, 0>(skbuffer, 13916, ckbuffer, 26450, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 3423, 0, 840, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 3633, 70, 945, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 3843, 140, 1050, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 210, 12558, 12768, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 1155, 12768, 13083, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 2100, 13083, 13524, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 4053, 0, 210, 1155, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 5943, 840, 1155, 2100, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 8778, 3423, 4053, 5943, r_ab, 3, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 0, skbuffer, 8778, 3, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 245, skbuffer, 9198, 3, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 490, skbuffer, 9618, 3, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 735, skbuffer, 10038, 3, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 980, skbuffer, 10458, 3, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 1225, skbuffer, 10878, 3, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 1470, skbuffer, 11298, 3, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 1715, skbuffer, 11718, 3, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 1960, skbuffer, 12138, 3, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 3, 3, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDFFS_hpp */
