#ifndef ElectronRepulsionGeom1010RecFFDP_hpp
#define ElectronRepulsionGeom1010RecFFDP_hpp

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
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPHXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSIXX.hpp"
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
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSF.hpp"
#include "ElectronRepulsionPrimRecSKSG.hpp"
#include "ElectronRepulsionPrimRecSKSP.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FF|1/|r-r'||DP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_ffdp(T& distributor,
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

    CSimdArray<double> pbuffer(11466, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(6708, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(31356, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(75465, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(6615, 1);

    // setup Boys fuction data

    const CBoysFunc<11> bf_table;

    CSimdArray<double> bf_data(13, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 12, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 12, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 12, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 12);
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

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 11, pfactors, 16, bf_data, 11);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 12, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 21, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 24, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 27, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 30, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 33, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 36, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 39, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 42, 10, 11, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 45, 0, 1, 12, 15, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 51, 1, 2, 15, 18, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 57, 2, 3, 18, 21, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 63, 3, 4, 21, 24, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 69, 4, 5, 24, 27, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 75, 5, 6, 27, 30, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 81, 6, 7, 30, 33, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 87, 7, 8, 33, 36, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 93, 8, 9, 36, 39, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 99, 9, 10, 39, 42, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 12, 15, 45, 51, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 15, 18, 51, 57, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 18, 21, 57, 63, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 21, 24, 63, 69, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 24, 27, 69, 75, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 27, 30, 75, 81, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 30, 33, 81, 87, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 175, 33, 36, 87, 93, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 185, 36, 39, 93, 99, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 195, 45, 51, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 210, 51, 57, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 225, 57, 63, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 240, 63, 69, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 255, 69, 75, 145, 155, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 270, 75, 81, 155, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 285, 81, 87, 165, 175, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 300, 87, 93, 175, 185, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 315, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 318, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 321, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 324, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 327, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 330, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 333, 1, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 342, 2, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 351, 3, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 360, 4, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 369, 5, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 378, 6, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 387, 7, 30, 33, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 396, 15, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 414, 18, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 432, 21, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 450, 24, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 468, 27, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 486, 30, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 504, 33, 81, 87, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 522, 51, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 552, 57, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 582, 63, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 612, 69, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 642, 75, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 672, 81, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 702, 87, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 732, 115, 195, 210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 777, 125, 210, 225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 822, 135, 225, 240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 867, 145, 240, 255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 912, 155, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 957, 165, 270, 285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1002, 175, 285, 300, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1047, 1, 2, 315, 318, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1053, 2, 3, 318, 321, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1059, 3, 4, 321, 324, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1065, 4, 5, 324, 327, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1071, 5, 6, 327, 330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1077, 12, 15, 315, 333, 342, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1095, 15, 18, 318, 342, 351, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1113, 18, 21, 321, 351, 360, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1131, 21, 24, 324, 360, 369, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1149, 24, 27, 327, 369, 378, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1167, 27, 30, 330, 378, 387, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1185, 45, 51, 342, 396, 414, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1221, 51, 57, 351, 414, 432, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1257, 57, 63, 360, 432, 450, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1293, 63, 69, 369, 450, 468, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1329, 69, 75, 378, 468, 486, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1365, 75, 81, 387, 486, 504, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1401, 105, 115, 414, 522, 552, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1461, 115, 125, 432, 552, 582, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1521, 125, 135, 450, 582, 612, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1581, 135, 145, 468, 612, 642, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1641, 145, 155, 486, 642, 672, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1701, 155, 165, 504, 672, 702, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1761, 195, 210, 552, 732, 777, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1851, 210, 225, 582, 777, 822, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1941, 225, 240, 612, 822, 867, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2031, 240, 255, 642, 867, 912, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2121, 255, 270, 672, 912, 957, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2211, 270, 285, 702, 957, 1002, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2301, 315, 318, 1047, 1053, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2311, 318, 321, 1053, 1059, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2321, 321, 324, 1059, 1065, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2331, 324, 327, 1065, 1071, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2341, 333, 342, 1047, 1077, 1095, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2371, 342, 351, 1053, 1095, 1113, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2401, 351, 360, 1059, 1113, 1131, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2431, 360, 369, 1065, 1131, 1149, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2461, 369, 378, 1071, 1149, 1167, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2491, 396, 414, 1095, 1185, 1221, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2551, 414, 432, 1113, 1221, 1257, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2611, 432, 450, 1131, 1257, 1293, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2671, 450, 468, 1149, 1293, 1329, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2731, 468, 486, 1167, 1329, 1365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2791, 522, 552, 1221, 1401, 1461, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2891, 552, 582, 1257, 1461, 1521, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2991, 582, 612, 1293, 1521, 1581, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3091, 612, 642, 1329, 1581, 1641, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3191, 642, 672, 1365, 1641, 1701, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3291, 732, 777, 1461, 1761, 1851, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3441, 777, 822, 1521, 1851, 1941, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3591, 822, 867, 1581, 1941, 2031, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3741, 867, 912, 1641, 2031, 2121, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3891, 912, 957, 1701, 2121, 2211, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4041, 1047, 1053, 2301, 2311, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4056, 1053, 1059, 2311, 2321, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4071, 1059, 1065, 2321, 2331, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4086, 1077, 1095, 2301, 2341, 2371, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4131, 1095, 1113, 2311, 2371, 2401, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4176, 1113, 1131, 2321, 2401, 2431, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4221, 1131, 1149, 2331, 2431, 2461, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4266, 1185, 1221, 2371, 2491, 2551, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4356, 1221, 1257, 2401, 2551, 2611, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4446, 1257, 1293, 2431, 2611, 2671, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4536, 1293, 1329, 2461, 2671, 2731, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4626, 1401, 1461, 2551, 2791, 2891, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4776, 1461, 1521, 2611, 2891, 2991, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4926, 1521, 1581, 2671, 2991, 3091, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5076, 1581, 1641, 2731, 3091, 3191, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5226, 1761, 1851, 2891, 3291, 3441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5451, 1851, 1941, 2991, 3441, 3591, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5676, 1941, 2031, 3091, 3591, 3741, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5901, 2031, 2121, 3191, 3741, 3891, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 6126, 2301, 2311, 4041, 4056, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 6147, 2311, 2321, 4056, 4071, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 6168, 2341, 2371, 4041, 4086, 4131, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 6231, 2371, 2401, 4056, 4131, 4176, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 6294, 2401, 2431, 4071, 4176, 4221, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 6357, 2491, 2551, 4131, 4266, 4356, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 6483, 2551, 2611, 4176, 4356, 4446, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 6609, 2611, 2671, 4221, 4446, 4536, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 6735, 2791, 2891, 4356, 4626, 4776, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 6945, 2891, 2991, 4446, 4776, 4926, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 7155, 2991, 3091, 4536, 4926, 5076, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 7365, 3291, 3441, 4776, 5226, 5451, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 7680, 3441, 3591, 4926, 5451, 5676, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 7995, 3591, 3741, 5076, 5676, 5901, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 8310, 4041, 4056, 6126, 6147, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 8338, 4086, 4131, 6126, 6168, 6231, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 8422, 4131, 4176, 6147, 6231, 6294, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 8506, 4266, 4356, 6231, 6357, 6483, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 8674, 4356, 4446, 6294, 6483, 6609, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 8842, 4626, 4776, 6483, 6735, 6945, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 9122, 4776, 4926, 6609, 6945, 7155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 9402, 5226, 5451, 6945, 7365, 7680, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 9822, 5451, 5676, 7155, 7680, 7995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 10242, 6168, 6231, 8310, 8338, 8422, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 10350, 6357, 6483, 8422, 8506, 8674, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 10566, 6735, 6945, 8674, 8842, 9122, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 10926, 7365, 7680, 9122, 9402, 9822, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 2341, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 2491, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 4086, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 135, pbuffer, 4266, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 225, pbuffer, 6168, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 288, pbuffer, 6357, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {2341, 2371});

                pbuffer.scale(2.0 * a_exp, {2491, 2551});

                pbuffer.scale(2.0 * a_exp, {4086, 4131});

                pbuffer.scale(2.0 * a_exp, {4266, 4356});

                pbuffer.scale(2.0 * a_exp, {6168, 6231});

                pbuffer.scale(2.0 * a_exp, {6357, 6483});

                pbuffer.scale(2.0 * a_exp, {8338, 8422});

                pbuffer.scale(2.0 * a_exp, {8506, 8674});

                pbuffer.scale(2.0 * a_exp, {10242, 10350});

                pbuffer.scale(2.0 * a_exp, {10350, 10566});

                t2cfunc::reduce(cbuffer, 1978, pbuffer, 2341, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2008, pbuffer, 2491, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2068, pbuffer, 4086, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2113, pbuffer, 4266, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2203, pbuffer, 6168, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2266, pbuffer, 6357, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2392, pbuffer, 8338, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2476, pbuffer, 8506, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2644, pbuffer, 10242, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2752, pbuffer, 10350, 216, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2341, 2371});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2491, 2551});

                pbuffer.scale(pfactors, 0, 2.0, {2791, 2891});

                pbuffer.scale(pfactors, 0, 2.0, {3291, 3441});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {4086, 4131});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {4266, 4356});

                pbuffer.scale(pfactors, 0, 2.0, {4626, 4776});

                pbuffer.scale(pfactors, 0, 2.0, {5226, 5451});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {6168, 6231});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {6357, 6483});

                pbuffer.scale(pfactors, 0, 2.0, {6735, 6945});

                pbuffer.scale(pfactors, 0, 2.0, {7365, 7680});

                t2cfunc::reduce(cbuffer, 414, pbuffer, 2341, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 444, pbuffer, 2491, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 504, pbuffer, 2791, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 604, pbuffer, 3291, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 754, pbuffer, 4086, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 799, pbuffer, 4266, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 889, pbuffer, 4626, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1039, pbuffer, 5226, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1264, pbuffer, 6168, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1327, pbuffer, 6357, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1453, pbuffer, 6735, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1663, pbuffer, 7365, 315, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {2341, 2371});

                pbuffer.scale(2.0 * a_exp, {2491, 2551});

                pbuffer.scale(2.0 * a_exp, {2791, 2891});

                pbuffer.scale(2.0 * a_exp, {3291, 3441});

                pbuffer.scale(2.0 * a_exp, {4086, 4131});

                pbuffer.scale(2.0 * a_exp, {4266, 4356});

                pbuffer.scale(2.0 * a_exp, {4626, 4776});

                pbuffer.scale(2.0 * a_exp, {5226, 5451});

                pbuffer.scale(2.0 * a_exp, {6168, 6231});

                pbuffer.scale(2.0 * a_exp, {6357, 6483});

                pbuffer.scale(2.0 * a_exp, {6735, 6945});

                pbuffer.scale(2.0 * a_exp, {7365, 7680});

                pbuffer.scale(pfactors, 0, 2.0, {8338, 8422});

                pbuffer.scale(pfactors, 0, 2.0, {8506, 8674});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {8842, 9122});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {9402, 9822});

                pbuffer.scale(pfactors, 0, 2.0, {10242, 10350});

                pbuffer.scale(pfactors, 0, 2.0, {10350, 10566});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {10566, 10926});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {10926, 11466});

                t2cfunc::reduce(cbuffer, 2968, pbuffer, 2341, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2998, pbuffer, 2491, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3058, pbuffer, 2791, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3158, pbuffer, 3291, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3308, pbuffer, 4086, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3353, pbuffer, 4266, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3443, pbuffer, 4626, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3593, pbuffer, 5226, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3818, pbuffer, 6168, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3881, pbuffer, 6357, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4007, pbuffer, 6735, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4217, pbuffer, 7365, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4532, pbuffer, 8338, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4616, pbuffer, 8506, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4784, pbuffer, 8842, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5064, pbuffer, 9402, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5484, pbuffer, 10242, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5592, pbuffer, 10350, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5808, pbuffer, 10566, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6168, pbuffer, 10926, 540, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 570, cbuffer, 0, 30, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2865, cbuffer, 90, 135, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 6222, cbuffer, 225, 288, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 9816, cbuffer, 1978, 2008, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 12111, cbuffer, 2068, 2113, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 15468, cbuffer, 2203, 2266, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 20088, cbuffer, 2392, 2476, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 26172, cbuffer, 2644, 2752, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 414, 444, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 90, cbuffer, 444, 504, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 270, cbuffer, 504, 604, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 660, cbuffer, 0, 0, 90, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 930, cbuffer, 30, 90, 270, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1470, 570, 660, 930, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2010, cbuffer, 754, 799, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2145, cbuffer, 799, 889, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 2415, cbuffer, 889, 1039, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3000, cbuffer, 90, 2010, 2145, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 3405, cbuffer, 135, 2145, 2415, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 4215, 2865, 3000, 3405, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5025, cbuffer, 1264, 1327, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5214, cbuffer, 1327, 1453, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 5592, cbuffer, 1453, 1663, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 6411, cbuffer, 225, 5025, 5214, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 6978, cbuffer, 288, 5214, 5592, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 8112, 6222, 6411, 6978, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9246, cbuffer, 2968, 2998, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9336, cbuffer, 2998, 3058, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 9516, cbuffer, 3058, 3158, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 9906, cbuffer, 1978, 9246, 9336, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 10176, cbuffer, 2008, 9336, 9516, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 10716, 9816, 9906, 10176, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 11256, cbuffer, 3308, 3353, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 11391, cbuffer, 3353, 3443, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 11661, cbuffer, 3443, 3593, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 12246, cbuffer, 2068, 11256, 11391, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 12651, cbuffer, 2113, 11391, 11661, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 13461, 12111, 12246, 12651, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 14271, cbuffer, 3818, 3881, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 14460, cbuffer, 3881, 4007, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 14838, cbuffer, 4007, 4217, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 15657, cbuffer, 2203, 14271, 14460, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 16224, cbuffer, 2266, 14460, 14838, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 17358, 15468, 15657, 16224, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 18492, cbuffer, 4532, 4616, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 18744, cbuffer, 4616, 4784, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 19248, cbuffer, 4784, 5064, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 20340, cbuffer, 2392, 18492, 18744, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 21096, cbuffer, 2476, 18744, 19248, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 22608, 20088, 20340, 21096, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 24120, cbuffer, 5484, 5592, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 24444, cbuffer, 5592, 5808, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 25092, cbuffer, 5808, 6168, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 26496, cbuffer, 2644, 24120, 24444, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 27468, cbuffer, 2752, 24444, 25092, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 29412, 26172, 26496, 27468, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<2, 1>(skbuffer, 0, ckbuffer, 1470, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 150, ckbuffer, 1650, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 300, ckbuffer, 1830, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 1800, ckbuffer, 4215, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 2025, ckbuffer, 4485, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 2250, ckbuffer, 4755, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 4500, ckbuffer, 8112, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 4815, ckbuffer, 8490, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 5130, ckbuffer, 8868, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 70515, ckbuffer, 10716, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 70665, ckbuffer, 10896, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 70815, ckbuffer, 11076, 0, 3);

            t4cfunc::ket_transform<2, 1>(skbuffer, 70965, ckbuffer, 13461, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 71190, ckbuffer, 13731, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 71415, ckbuffer, 14001, 0, 4);

            t4cfunc::ket_transform<2, 1>(skbuffer, 71640, ckbuffer, 17358, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 71955, ckbuffer, 17736, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 72270, ckbuffer, 18114, 0, 5);

            t4cfunc::ket_transform<2, 1>(skbuffer, 72585, ckbuffer, 22608, 0, 6);

            t4cfunc::ket_transform<2, 1>(skbuffer, 73005, ckbuffer, 23112, 0, 6);

            t4cfunc::ket_transform<2, 1>(skbuffer, 73425, ckbuffer, 23616, 0, 6);

            t4cfunc::ket_transform<2, 1>(skbuffer, 73845, ckbuffer, 29412, 0, 7);

            t4cfunc::ket_transform<2, 1>(skbuffer, 74385, ckbuffer, 30060, 0, 7);

            t4cfunc::ket_transform<2, 1>(skbuffer, 74925, ckbuffer, 30708, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 12060, 0, 1800, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 12510, 150, 2025, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 12960, 300, 2250, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 17460, 1800, 4500, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 18135, 2025, 4815, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 18810, 2250, 5130, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 34065, 12060, 17460, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 34965, 12510, 18135, r_ab, 2, 1);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 35865, 12960, 18810, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 450, 70515, 70965, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 2475, 70965, 71640, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 5445, 71640, 72585, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sixx(skbuffer, 8280, 72585, 73845, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 13410, 0, 450, 2475, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 19485, 1800, 2475, 5445, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_phxx(skbuffer, 25560, 4500, 5445, 8280, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 36765, 12060, 13410, 19485, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dgxx(skbuffer, 44865, 17460, 19485, 25560, r_ab, 2, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ffxx(skbuffer, 57015, 34065, 36765, 44865, r_ab, 2, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 57015, 2, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 735, skbuffer, 58515, 2, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1470, skbuffer, 60015, 2, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2205, skbuffer, 61515, 2, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2940, skbuffer, 63015, 2, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 3675, skbuffer, 64515, 2, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 4410, skbuffer, 66015, 2, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 5145, skbuffer, 67515, 2, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 5880, skbuffer, 69015, 2, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 2, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFFDP_hpp */
