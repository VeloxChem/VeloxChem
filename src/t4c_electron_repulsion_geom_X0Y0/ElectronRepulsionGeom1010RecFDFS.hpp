#ifndef ElectronRepulsionGeom1010RecFDFS_hpp
#define ElectronRepulsionGeom1010RecFDFS_hpp

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
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FD|1/|r-r'||FS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fdfs(T& distributor,
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

    CSimdArray<double> cbuffer(4995, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(29970, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(23058, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 910, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 940, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 24, pbuffer, 1030, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 1960, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 70, pbuffer, 2000, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 2120, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 3360, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 175, pbuffer, 3405, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 220, pbuffer, 3540, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {910, 916});

                pbuffer.scale(2.0 * a_exp, {940, 958});

                pbuffer.scale(2.0 * a_exp, {1030, 1066});

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

                t2cfunc::reduce(cbuffer, 1395, pbuffer, 910, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1401, pbuffer, 940, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1419, pbuffer, 1030, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1455, pbuffer, 1960, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1465, pbuffer, 2000, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1495, pbuffer, 2120, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1555, pbuffer, 3360, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1570, pbuffer, 3405, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1615, pbuffer, 3540, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1705, pbuffer, 4935, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1726, pbuffer, 4977, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1789, pbuffer, 5103, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1915, pbuffer, 6405, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1943, pbuffer, 6433, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2027, pbuffer, 6517, 168, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {910, 916});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {940, 958});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1030, 1066});

                pbuffer.scale(pfactors, 0, 2.0, {1210, 1270});

                pbuffer.scale(pfactors, 0, 2.0, {1510, 1600});

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

                t2cfunc::reduce(cbuffer, 310, pbuffer, 910, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 316, pbuffer, 940, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 334, pbuffer, 1030, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 370, pbuffer, 1210, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 430, pbuffer, 1510, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 520, pbuffer, 1960, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 530, pbuffer, 2000, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 560, pbuffer, 2120, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 620, pbuffer, 2360, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 720, pbuffer, 2760, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 870, pbuffer, 3360, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 885, pbuffer, 3405, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 930, pbuffer, 3540, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1020, pbuffer, 3810, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1170, pbuffer, 4260, 225, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {910, 916});

                pbuffer.scale(2.0 * a_exp, {940, 958});

                pbuffer.scale(2.0 * a_exp, {1030, 1066});

                pbuffer.scale(2.0 * a_exp, {1210, 1270});

                pbuffer.scale(2.0 * a_exp, {1510, 1600});

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

                t2cfunc::reduce(cbuffer, 2195, pbuffer, 910, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2201, pbuffer, 940, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2219, pbuffer, 1030, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2255, pbuffer, 1210, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2315, pbuffer, 1510, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2405, pbuffer, 1960, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2415, pbuffer, 2000, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2445, pbuffer, 2120, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2505, pbuffer, 2360, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2605, pbuffer, 2760, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2755, pbuffer, 3360, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2770, pbuffer, 3405, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2815, pbuffer, 3540, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2905, pbuffer, 3810, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3055, pbuffer, 4260, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3280, pbuffer, 4935, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3301, pbuffer, 4977, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3364, pbuffer, 5103, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3490, pbuffer, 5355, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3700, pbuffer, 5775, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4015, pbuffer, 6405, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4043, pbuffer, 6433, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4127, pbuffer, 6517, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4295, pbuffer, 6685, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4575, pbuffer, 6965, 420, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 360, cbuffer, 0, 6, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 432, cbuffer, 6, 24, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 972, 360, 432, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2220, cbuffer, 60, 70, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2340, cbuffer, 70, 100, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 3240, 2220, 2340, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 5220, cbuffer, 160, 175, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 5400, cbuffer, 175, 220, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 6750, 5220, 5400, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 8730, cbuffer, 1395, 1401, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 8802, cbuffer, 1401, 1419, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 9342, 8730, 8802, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 10590, cbuffer, 1455, 1465, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 10710, cbuffer, 1465, 1495, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 11610, 10590, 10710, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 13590, cbuffer, 1555, 1570, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 13770, cbuffer, 1570, 1615, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 15120, 13590, 13770, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 18000, cbuffer, 1705, 1726, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 18252, cbuffer, 1726, 1789, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 20142, 18000, 18252, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 24090, cbuffer, 1915, 1943, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 24426, cbuffer, 1943, 2027, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 26946, 24090, 24426, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 310, 316, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 18, cbuffer, 316, 334, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 72, cbuffer, 334, 370, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 180, cbuffer, 370, 430, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 378, cbuffer, 0, 0, 18, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 486, cbuffer, 6, 18, 72, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 648, cbuffer, 24, 72, 180, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1008, 360, 378, 486, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1116, 432, 486, 648, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 1440, 972, 1008, 1116, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1620, cbuffer, 520, 530, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1650, cbuffer, 530, 560, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1740, cbuffer, 560, 620, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1920, cbuffer, 620, 720, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2250, cbuffer, 60, 1620, 1650, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2430, cbuffer, 70, 1650, 1740, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2700, cbuffer, 100, 1740, 1920, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 3300, 2220, 2250, 2430, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 3480, 2340, 2430, 2700, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 4020, 3240, 3300, 3480, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 4320, cbuffer, 870, 885, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4365, cbuffer, 885, 930, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4500, cbuffer, 930, 1020, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 4770, cbuffer, 1020, 1170, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 5265, cbuffer, 160, 4320, 4365, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 5535, cbuffer, 175, 4365, 4500, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 5940, cbuffer, 220, 4500, 4770, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 6840, 5220, 5265, 5535, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 7110, 5400, 5535, 5940, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 7920, 6750, 6840, 7110, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 8370, cbuffer, 2195, 2201, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 8388, cbuffer, 2201, 2219, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 8442, cbuffer, 2219, 2255, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 8550, cbuffer, 2255, 2315, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 8748, cbuffer, 1395, 8370, 8388, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 8856, cbuffer, 1401, 8388, 8442, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 9018, cbuffer, 1419, 8442, 8550, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 9378, 8730, 8748, 8856, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 9486, 8802, 8856, 9018, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 9810, 9342, 9378, 9486, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 9990, cbuffer, 2405, 2415, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 10020, cbuffer, 2415, 2445, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 10110, cbuffer, 2445, 2505, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 10290, cbuffer, 2505, 2605, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 10620, cbuffer, 1455, 9990, 10020, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 10800, cbuffer, 1465, 10020, 10110, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 11070, cbuffer, 1495, 10110, 10290, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 11670, 10590, 10620, 10800, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 11850, 10710, 10800, 11070, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 12390, 11610, 11670, 11850, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 12690, cbuffer, 2755, 2770, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 12735, cbuffer, 2770, 2815, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 12870, cbuffer, 2815, 2905, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 13140, cbuffer, 2905, 3055, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 13635, cbuffer, 1555, 12690, 12735, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 13905, cbuffer, 1570, 12735, 12870, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 14310, cbuffer, 1615, 12870, 13140, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 15210, 13590, 13635, 13905, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 15480, 13770, 13905, 14310, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 16290, 15120, 15210, 15480, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 16740, cbuffer, 3280, 3301, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 16803, cbuffer, 3301, 3364, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 16992, cbuffer, 3364, 3490, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 17370, cbuffer, 3490, 3700, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 18063, cbuffer, 1705, 16740, 16803, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 18441, cbuffer, 1726, 16803, 16992, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 19008, cbuffer, 1789, 16992, 17370, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 20268, 18000, 18063, 18441, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 20646, 18252, 18441, 19008, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 21780, 20142, 20268, 20646, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 22410, cbuffer, 4015, 4043, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 22494, cbuffer, 4043, 4127, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 22746, cbuffer, 4127, 4295, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 23250, cbuffer, 4295, 4575, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 24174, cbuffer, 1915, 22410, 22494, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 24678, cbuffer, 1943, 22494, 22746, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 25434, cbuffer, 2027, 22746, 23250, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 27114, 24090, 24174, 24678, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 27618, 24426, 24678, 25434, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 29130, 26946, 27114, 27618, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<3, 0>(skbuffer, 0, ckbuffer, 1440, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 42, ckbuffer, 1500, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 84, ckbuffer, 1560, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 504, ckbuffer, 4020, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 574, ckbuffer, 4120, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 644, ckbuffer, 4220, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1344, ckbuffer, 7920, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1449, ckbuffer, 8070, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1554, ckbuffer, 8220, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21378, ckbuffer, 9810, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21420, ckbuffer, 9870, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21462, ckbuffer, 9930, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21504, ckbuffer, 12390, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21574, ckbuffer, 12490, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21644, ckbuffer, 12590, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21714, ckbuffer, 16290, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21819, ckbuffer, 16440, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21924, ckbuffer, 16590, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 22029, ckbuffer, 21780, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 22176, ckbuffer, 21990, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 22323, ckbuffer, 22200, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 22470, ckbuffer, 29130, 0, 6);

            t4cfunc::ket_transform<3, 0>(skbuffer, 22666, ckbuffer, 29410, 0, 6);

            t4cfunc::ket_transform<3, 0>(skbuffer, 22862, ckbuffer, 29690, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 3927, 0, 504, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 4053, 42, 574, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 4179, 84, 644, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 5439, 504, 1344, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 5649, 574, 1449, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 5859, 644, 1554, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 10794, 3927, 5439, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 11046, 4053, 5649, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 11298, 4179, 5859, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 126, 21378, 21504, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 714, 21504, 21714, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 1659, 21714, 22029, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 2604, 22029, 22470, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 4305, 0, 126, 714, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 6069, 504, 714, 1659, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 7959, 1344, 1659, 2604, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 11550, 3927, 4305, 6069, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 13818, 5439, 6069, 7959, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fdxx(skbuffer, 17598, 10794, 11550, 13818, r_ab, 3, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 17598, 3, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 245, skbuffer, 18018, 3, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 490, skbuffer, 18438, 3, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 735, skbuffer, 18858, 3, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 980, skbuffer, 19278, 3, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1225, skbuffer, 19698, 3, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1470, skbuffer, 20118, 3, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1715, skbuffer, 20538, 3, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1960, skbuffer, 20958, 3, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 3, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFDFS_hpp */
