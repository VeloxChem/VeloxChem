#ifndef ElectronRepulsionGeom1010RecFDDP_hpp
#define ElectronRepulsionGeom1010RecFDDP_hpp

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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FD|1/|r-r'||DP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fddp(T& distributor,
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

    CSimdArray<double> pbuffer(7302, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(4773, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(22311, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(49410, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(4725, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 280, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 283, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 286, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 289, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 292, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 295, 1, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 304, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 313, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 322, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 331, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 340, 6, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 349, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 367, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 385, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 403, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 421, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 439, 29, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 457, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 487, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 517, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 547, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 577, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 607, 77, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 637, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 682, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 727, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 772, 135, 220, 235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 817, 145, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 862, 155, 250, 265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 907, 1, 2, 280, 283, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 913, 2, 3, 283, 286, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 919, 3, 4, 286, 289, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 925, 4, 5, 289, 292, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 931, 11, 14, 280, 295, 304, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 949, 14, 17, 283, 304, 313, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 967, 17, 20, 286, 313, 322, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 985, 20, 23, 289, 322, 331, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1003, 23, 26, 292, 331, 340, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1021, 41, 47, 304, 349, 367, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1057, 47, 53, 313, 367, 385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1093, 53, 59, 322, 385, 403, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1129, 59, 65, 331, 403, 421, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1165, 65, 71, 340, 421, 439, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1201, 95, 105, 367, 457, 487, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1261, 105, 115, 385, 487, 517, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1321, 115, 125, 403, 517, 547, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1381, 125, 135, 421, 547, 577, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1441, 135, 145, 439, 577, 607, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1501, 175, 190, 487, 637, 682, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1591, 190, 205, 517, 682, 727, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1681, 205, 220, 547, 727, 772, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1771, 220, 235, 577, 772, 817, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1861, 235, 250, 607, 817, 862, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1951, 280, 283, 907, 913, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1961, 283, 286, 913, 919, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1971, 286, 289, 919, 925, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1981, 295, 304, 907, 931, 949, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2011, 304, 313, 913, 949, 967, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2041, 313, 322, 919, 967, 985, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2071, 322, 331, 925, 985, 1003, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2101, 349, 367, 949, 1021, 1057, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2161, 367, 385, 967, 1057, 1093, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2221, 385, 403, 985, 1093, 1129, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2281, 403, 421, 1003, 1129, 1165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2341, 457, 487, 1057, 1201, 1261, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2441, 487, 517, 1093, 1261, 1321, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2541, 517, 547, 1129, 1321, 1381, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2641, 547, 577, 1165, 1381, 1441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2741, 637, 682, 1261, 1501, 1591, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2891, 682, 727, 1321, 1591, 1681, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3041, 727, 772, 1381, 1681, 1771, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3191, 772, 817, 1441, 1771, 1861, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 3341, 907, 913, 1951, 1961, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 3356, 913, 919, 1961, 1971, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3371, 931, 949, 1951, 1981, 2011, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3416, 949, 967, 1961, 2011, 2041, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 3461, 967, 985, 1971, 2041, 2071, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3506, 1021, 1057, 2011, 2101, 2161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3596, 1057, 1093, 2041, 2161, 2221, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 3686, 1093, 1129, 2071, 2221, 2281, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3776, 1201, 1261, 2161, 2341, 2441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3926, 1261, 1321, 2221, 2441, 2541, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4076, 1321, 1381, 2281, 2541, 2641, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4226, 1501, 1591, 2441, 2741, 2891, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4451, 1591, 1681, 2541, 2891, 3041, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4676, 1681, 1771, 2641, 3041, 3191, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 4901, 1951, 1961, 3341, 3356, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 4922, 1981, 2011, 3341, 3371, 3416, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 4985, 2011, 2041, 3356, 3416, 3461, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 5048, 2101, 2161, 3416, 3506, 3596, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 5174, 2161, 2221, 3461, 3596, 3686, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5300, 2341, 2441, 3596, 3776, 3926, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5510, 2441, 2541, 3686, 3926, 4076, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 5720, 2741, 2891, 3926, 4226, 4451, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 6035, 2891, 3041, 4076, 4451, 4676, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 6350, 3371, 3416, 4901, 4922, 4985, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 6434, 3506, 3596, 4985, 5048, 5174, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 6602, 3776, 3926, 5174, 5300, 5510, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 6882, 4226, 4451, 5510, 5720, 6035, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 931, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 1021, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 54, pbuffer, 1981, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 84, pbuffer, 2101, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 144, pbuffer, 3371, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 189, pbuffer, 3506, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {931, 949});

                pbuffer.scale(2.0 * a_exp, {1021, 1057});

                pbuffer.scale(2.0 * a_exp, {1981, 2011});

                pbuffer.scale(2.0 * a_exp, {2101, 2161});

                pbuffer.scale(2.0 * a_exp, {3371, 3416});

                pbuffer.scale(2.0 * a_exp, {3506, 3596});

                pbuffer.scale(2.0 * a_exp, {4922, 4985});

                pbuffer.scale(2.0 * a_exp, {5048, 5174});

                pbuffer.scale(2.0 * a_exp, {6350, 6434});

                pbuffer.scale(2.0 * a_exp, {6434, 6602});

                t2cfunc::reduce(cbuffer, 1333, pbuffer, 931, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1351, pbuffer, 1021, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1387, pbuffer, 1981, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1417, pbuffer, 2101, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1477, pbuffer, 3371, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1522, pbuffer, 3506, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1612, pbuffer, 4922, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1675, pbuffer, 5048, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1801, pbuffer, 6350, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1885, pbuffer, 6434, 168, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {931, 949});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1021, 1057});

                pbuffer.scale(pfactors, 0, 2.0, {1201, 1261});

                pbuffer.scale(pfactors, 0, 2.0, {1501, 1591});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1981, 2011});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2101, 2161});

                pbuffer.scale(pfactors, 0, 2.0, {2341, 2441});

                pbuffer.scale(pfactors, 0, 2.0, {2741, 2891});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3371, 3416});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3506, 3596});

                pbuffer.scale(pfactors, 0, 2.0, {3776, 3926});

                pbuffer.scale(pfactors, 0, 2.0, {4226, 4451});

                t2cfunc::reduce(cbuffer, 279, pbuffer, 931, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 297, pbuffer, 1021, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 333, pbuffer, 1201, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 393, pbuffer, 1501, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 483, pbuffer, 1981, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 513, pbuffer, 2101, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 573, pbuffer, 2341, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 673, pbuffer, 2741, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 823, pbuffer, 3371, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 868, pbuffer, 3506, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 958, pbuffer, 3776, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1108, pbuffer, 4226, 225, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {931, 949});

                pbuffer.scale(2.0 * a_exp, {1021, 1057});

                pbuffer.scale(2.0 * a_exp, {1201, 1261});

                pbuffer.scale(2.0 * a_exp, {1501, 1591});

                pbuffer.scale(2.0 * a_exp, {1981, 2011});

                pbuffer.scale(2.0 * a_exp, {2101, 2161});

                pbuffer.scale(2.0 * a_exp, {2341, 2441});

                pbuffer.scale(2.0 * a_exp, {2741, 2891});

                pbuffer.scale(2.0 * a_exp, {3371, 3416});

                pbuffer.scale(2.0 * a_exp, {3506, 3596});

                pbuffer.scale(2.0 * a_exp, {3776, 3926});

                pbuffer.scale(2.0 * a_exp, {4226, 4451});

                pbuffer.scale(pfactors, 0, 2.0, {4922, 4985});

                pbuffer.scale(pfactors, 0, 2.0, {5048, 5174});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5300, 5510});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5720, 6035});

                pbuffer.scale(pfactors, 0, 2.0, {6350, 6434});

                pbuffer.scale(pfactors, 0, 2.0, {6434, 6602});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {6602, 6882});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {6882, 7302});

                t2cfunc::reduce(cbuffer, 2053, pbuffer, 931, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2071, pbuffer, 1021, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2107, pbuffer, 1201, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2167, pbuffer, 1501, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2257, pbuffer, 1981, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2287, pbuffer, 2101, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2347, pbuffer, 2341, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2447, pbuffer, 2741, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2597, pbuffer, 3371, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2642, pbuffer, 3506, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2732, pbuffer, 3776, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2882, pbuffer, 4226, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3107, pbuffer, 4922, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3170, pbuffer, 5048, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3296, pbuffer, 5300, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3506, pbuffer, 5720, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3821, pbuffer, 6350, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3905, pbuffer, 6434, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4073, pbuffer, 6602, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4353, pbuffer, 6882, 420, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 342, cbuffer, 0, 18, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1776, cbuffer, 54, 84, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4071, cbuffer, 144, 189, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6573, cbuffer, 1333, 1351, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 8007, cbuffer, 1387, 1417, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 10302, cbuffer, 1477, 1522, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 13659, cbuffer, 1612, 1675, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 18279, cbuffer, 1801, 1885, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 279, 297, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 54, cbuffer, 297, 333, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 162, cbuffer, 333, 393, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 396, cbuffer, 0, 0, 54, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 558, cbuffer, 18, 54, 162, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 882, 342, 396, 558, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1206, cbuffer, 483, 513, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1296, cbuffer, 513, 573, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1476, cbuffer, 573, 673, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1866, cbuffer, 54, 1206, 1296, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2136, cbuffer, 84, 1296, 1476, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 2676, 1776, 1866, 2136, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3216, cbuffer, 823, 868, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3351, cbuffer, 868, 958, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3621, cbuffer, 958, 1108, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4206, cbuffer, 144, 3216, 3351, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4611, cbuffer, 189, 3351, 3621, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 5421, 4071, 4206, 4611, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 6231, cbuffer, 2053, 2071, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6285, cbuffer, 2071, 2107, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6393, cbuffer, 2107, 2167, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 6627, cbuffer, 1333, 6231, 6285, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 6789, cbuffer, 1351, 6285, 6393, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 7113, 6573, 6627, 6789, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 7437, cbuffer, 2257, 2287, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 7527, cbuffer, 2287, 2347, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 7707, cbuffer, 2347, 2447, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 8097, cbuffer, 1387, 7437, 7527, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 8367, cbuffer, 1417, 7527, 7707, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 8907, 8007, 8097, 8367, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9447, cbuffer, 2597, 2642, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9582, cbuffer, 2642, 2732, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 9852, cbuffer, 2732, 2882, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 10437, cbuffer, 1477, 9447, 9582, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 10842, cbuffer, 1522, 9582, 9852, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 11652, 10302, 10437, 10842, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 12462, cbuffer, 3107, 3170, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 12651, cbuffer, 3170, 3296, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 13029, cbuffer, 3296, 3506, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 13848, cbuffer, 1612, 12462, 12651, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 14415, cbuffer, 1675, 12651, 13029, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 15549, 13659, 13848, 14415, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 16683, cbuffer, 3821, 3905, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 16935, cbuffer, 3905, 4073, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 17439, cbuffer, 4073, 4353, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 18531, cbuffer, 1801, 16683, 16935, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 19287, cbuffer, 1885, 16935, 17439, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 20799, 18279, 18531, 19287, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<2, 1>(skbuffer, 0, ckbuffer, 882, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 90, ckbuffer, 990, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 180, ckbuffer, 1098, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1080, ckbuffer, 2676, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1230, ckbuffer, 2856, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1380, ckbuffer, 3036, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 2880, ckbuffer, 5421, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 3105, ckbuffer, 5691, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 3330, ckbuffer, 5961, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 45810, ckbuffer, 7113, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 45900, ckbuffer, 7221, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 45990, ckbuffer, 7329, 0, 2);

            t4cfunc::ket_transform<2, 1>(skbuffer, 46080, ckbuffer, 8907, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 46230, ckbuffer, 9087, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 46380, ckbuffer, 9267, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 46530, ckbuffer, 11652, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 46755, ckbuffer, 11922, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 46980, ckbuffer, 12192, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 47205, ckbuffer, 15549, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 47520, ckbuffer, 15927, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 47835, ckbuffer, 16305, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 48150, ckbuffer, 20799, 0, 6);

            t4cfunc::ket_transform<2, 1>(skbuffer, 48570, ckbuffer, 21303, 0, 6);

            t4cfunc::ket_transform<2, 1>(skbuffer, 48990, ckbuffer, 21807, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 8415, 0, 1080, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 8685, 90, 1230, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 8955, 180, 1380, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 11655, 1080, 2880, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 12105, 1230, 3105, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 12555, 1380, 3330, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 23130, 8415, 11655, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 23670, 8685, 12105, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 24210, 8955, 12555, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 270, 45810, 46080, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 1530, 46080, 46530, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 3555, 46530, 47205, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 5580, 47205, 48150, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 9225, 0, 270, 1530, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 13005, 1080, 1530, 3555, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 17055, 2880, 3555, 5580, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 24750, 8415, 9225, 13005, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 29610, 11655, 13005, 17055, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fdxx(skbuffer, 37710, 23130, 24750, 29610, r_ab, 2, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 37710, 2, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 525, skbuffer, 38610, 2, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1050, skbuffer, 39510, 2, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1575, skbuffer, 40410, 2, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2100, skbuffer, 41310, 2, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2625, skbuffer, 42210, 2, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 3150, skbuffer, 43110, 2, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 3675, skbuffer, 44010, 2, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 4200, skbuffer, 44910, 2, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 2, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFDDP_hpp */
