#ifndef ElectronRepulsionGeom1010RecFFFS_hpp
#define ElectronRepulsionGeom1010RecFFFS_hpp

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
#include "ElectronRepulsionPrimRecSKSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FF|1/|r-r'||FS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fffs(T& distributor,
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

    CSimdArray<double> pbuffer(11585, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(7020, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(42120, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(35217, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3087, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 315, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 318, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 321, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 324, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 327, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 330, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 333, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 336, 1, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 345, 2, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 354, 3, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 363, 4, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 372, 5, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 381, 6, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 390, 7, 30, 33, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 399, 15, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 417, 18, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 435, 21, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 453, 24, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 471, 27, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 489, 30, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 507, 33, 81, 87, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 525, 51, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 555, 57, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 585, 63, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 615, 69, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 645, 75, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 675, 81, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 705, 87, 165, 175, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 735, 115, 195, 210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 780, 125, 210, 225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 825, 135, 225, 240, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 870, 145, 240, 255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 915, 155, 255, 270, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 960, 165, 270, 285, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 1005, 175, 285, 300, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1050, 0, 1, 315, 318, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1056, 1, 2, 318, 321, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1062, 2, 3, 321, 324, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1068, 3, 4, 324, 327, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1074, 4, 5, 327, 330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1080, 5, 6, 330, 333, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1086, 12, 15, 318, 336, 345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1104, 15, 18, 321, 345, 354, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1122, 18, 21, 324, 354, 363, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1140, 21, 24, 327, 363, 372, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1158, 24, 27, 330, 372, 381, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1176, 27, 30, 333, 381, 390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1194, 45, 51, 345, 399, 417, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1230, 51, 57, 354, 417, 435, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1266, 57, 63, 363, 435, 453, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1302, 63, 69, 372, 453, 471, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1338, 69, 75, 381, 471, 489, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1374, 75, 81, 390, 489, 507, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1410, 105, 115, 417, 525, 555, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1470, 115, 125, 435, 555, 585, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1530, 125, 135, 453, 585, 615, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1590, 135, 145, 471, 615, 645, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1650, 145, 155, 489, 645, 675, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1710, 155, 165, 507, 675, 705, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1770, 195, 210, 555, 735, 780, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1860, 210, 225, 585, 780, 825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1950, 225, 240, 615, 825, 870, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2040, 240, 255, 645, 870, 915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2130, 255, 270, 675, 915, 960, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 2220, 270, 285, 705, 960, 1005, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2310, 315, 318, 1050, 1056, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2320, 318, 321, 1056, 1062, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2330, 321, 324, 1062, 1068, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2340, 324, 327, 1068, 1074, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2350, 327, 330, 1074, 1080, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2360, 336, 345, 1056, 1086, 1104, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2390, 345, 354, 1062, 1104, 1122, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2420, 354, 363, 1068, 1122, 1140, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2450, 363, 372, 1074, 1140, 1158, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2480, 372, 381, 1080, 1158, 1176, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2510, 399, 417, 1104, 1194, 1230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2570, 417, 435, 1122, 1230, 1266, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2630, 435, 453, 1140, 1266, 1302, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2690, 453, 471, 1158, 1302, 1338, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2750, 471, 489, 1176, 1338, 1374, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2810, 525, 555, 1230, 1410, 1470, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2910, 555, 585, 1266, 1470, 1530, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3010, 585, 615, 1302, 1530, 1590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3110, 615, 645, 1338, 1590, 1650, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 3210, 645, 675, 1374, 1650, 1710, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3310, 735, 780, 1470, 1770, 1860, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3460, 780, 825, 1530, 1860, 1950, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3610, 825, 870, 1590, 1950, 2040, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3760, 870, 915, 1650, 2040, 2130, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3910, 915, 960, 1710, 2130, 2220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4060, 1050, 1056, 2310, 2320, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4075, 1056, 1062, 2320, 2330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4090, 1062, 1068, 2330, 2340, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 4105, 1068, 1074, 2340, 2350, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4120, 1086, 1104, 2320, 2360, 2390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4165, 1104, 1122, 2330, 2390, 2420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4210, 1122, 1140, 2340, 2420, 2450, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4255, 1140, 1158, 2350, 2450, 2480, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4300, 1194, 1230, 2390, 2510, 2570, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4390, 1230, 1266, 2420, 2570, 2630, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4480, 1266, 1302, 2450, 2630, 2690, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4570, 1302, 1338, 2480, 2690, 2750, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4660, 1410, 1470, 2570, 2810, 2910, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4810, 1470, 1530, 2630, 2910, 3010, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4960, 1530, 1590, 2690, 3010, 3110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 5110, 1590, 1650, 2750, 3110, 3210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5260, 1770, 1860, 2910, 3310, 3460, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5485, 1860, 1950, 3010, 3460, 3610, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5710, 1950, 2040, 3110, 3610, 3760, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 5935, 2040, 2130, 3210, 3760, 3910, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 6160, 2310, 2320, 4060, 4075, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 6181, 2320, 2330, 4075, 4090, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 6202, 2330, 2340, 4090, 4105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 6223, 2360, 2390, 4075, 4120, 4165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 6286, 2390, 2420, 4090, 4165, 4210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 6349, 2420, 2450, 4105, 4210, 4255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 6412, 2510, 2570, 4165, 4300, 4390, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 6538, 2570, 2630, 4210, 4390, 4480, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 6664, 2630, 2690, 4255, 4480, 4570, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 6790, 2810, 2910, 4390, 4660, 4810, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 7000, 2910, 3010, 4480, 4810, 4960, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 7210, 3010, 3110, 4570, 4960, 5110, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 7420, 3310, 3460, 4810, 5260, 5485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 7735, 3460, 3610, 4960, 5485, 5710, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 8050, 3610, 3760, 5110, 5710, 5935, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 8365, 4060, 4075, 6160, 6181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 8393, 4075, 4090, 6181, 6202, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 8421, 4120, 4165, 6181, 6223, 6286, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 8505, 4165, 4210, 6202, 6286, 6349, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 8589, 4300, 4390, 6286, 6412, 6538, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 8757, 4390, 4480, 6349, 6538, 6664, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 8925, 4660, 4810, 6538, 6790, 7000, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 9205, 4810, 4960, 6664, 7000, 7210, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 9485, 5260, 5485, 7000, 7420, 7735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisg(pbuffer, 9905, 5485, 5710, 7210, 7735, 8050, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_skss(pbuffer, 10325, 6160, 6181, 8365, 8393, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 10361, 6223, 6286, 8393, 8421, 8505, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 10469, 6412, 6538, 8505, 8589, 8757, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 10685, 6790, 7000, 8757, 8925, 9205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksg(pbuffer, 11045, 7420, 7735, 9205, 9485, 9905, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 2310, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 2360, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 40, pbuffer, 2510, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 4060, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 115, pbuffer, 4120, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 4300, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 250, pbuffer, 6160, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 271, pbuffer, 6223, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 334, pbuffer, 6412, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {2310, 2320});

                pbuffer.scale(2.0 * a_exp, {2360, 2390});

                pbuffer.scale(2.0 * a_exp, {2510, 2570});

                pbuffer.scale(2.0 * a_exp, {4060, 4075});

                pbuffer.scale(2.0 * a_exp, {4120, 4165});

                pbuffer.scale(2.0 * a_exp, {4300, 4390});

                pbuffer.scale(2.0 * a_exp, {6160, 6181});

                pbuffer.scale(2.0 * a_exp, {6223, 6286});

                pbuffer.scale(2.0 * a_exp, {6412, 6538});

                pbuffer.scale(2.0 * a_exp, {8365, 8393});

                pbuffer.scale(2.0 * a_exp, {8421, 8505});

                pbuffer.scale(2.0 * a_exp, {8589, 8757});

                pbuffer.scale(2.0 * a_exp, {10325, 10361});

                pbuffer.scale(2.0 * a_exp, {10361, 10469});

                pbuffer.scale(2.0 * a_exp, {10469, 10685});

                t2cfunc::reduce(cbuffer, 2070, pbuffer, 2310, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2080, pbuffer, 2360, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2110, pbuffer, 2510, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2170, pbuffer, 4060, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2185, pbuffer, 4120, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2230, pbuffer, 4300, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2320, pbuffer, 6160, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2341, pbuffer, 6223, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2404, pbuffer, 6412, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2530, pbuffer, 8365, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2558, pbuffer, 8421, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2642, pbuffer, 8589, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2810, pbuffer, 10325, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2846, pbuffer, 10361, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2954, pbuffer, 10469, 216, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2310, 2320});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2360, 2390});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2510, 2570});

                pbuffer.scale(pfactors, 0, 2.0, {2810, 2910});

                pbuffer.scale(pfactors, 0, 2.0, {3310, 3460});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {4060, 4075});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {4120, 4165});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {4300, 4390});

                pbuffer.scale(pfactors, 0, 2.0, {4660, 4810});

                pbuffer.scale(pfactors, 0, 2.0, {5260, 5485});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {6160, 6181});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {6223, 6286});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {6412, 6538});

                pbuffer.scale(pfactors, 0, 2.0, {6790, 7000});

                pbuffer.scale(pfactors, 0, 2.0, {7420, 7735});

                t2cfunc::reduce(cbuffer, 460, pbuffer, 2310, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 470, pbuffer, 2360, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 500, pbuffer, 2510, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 560, pbuffer, 2810, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 660, pbuffer, 3310, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 810, pbuffer, 4060, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 825, pbuffer, 4120, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 870, pbuffer, 4300, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 960, pbuffer, 4660, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1110, pbuffer, 5260, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1335, pbuffer, 6160, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1356, pbuffer, 6223, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1419, pbuffer, 6412, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1545, pbuffer, 6790, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1755, pbuffer, 7420, 315, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {2310, 2320});

                pbuffer.scale(2.0 * a_exp, {2360, 2390});

                pbuffer.scale(2.0 * a_exp, {2510, 2570});

                pbuffer.scale(2.0 * a_exp, {2810, 2910});

                pbuffer.scale(2.0 * a_exp, {3310, 3460});

                pbuffer.scale(2.0 * a_exp, {4060, 4075});

                pbuffer.scale(2.0 * a_exp, {4120, 4165});

                pbuffer.scale(2.0 * a_exp, {4300, 4390});

                pbuffer.scale(2.0 * a_exp, {4660, 4810});

                pbuffer.scale(2.0 * a_exp, {5260, 5485});

                pbuffer.scale(2.0 * a_exp, {6160, 6181});

                pbuffer.scale(2.0 * a_exp, {6223, 6286});

                pbuffer.scale(2.0 * a_exp, {6412, 6538});

                pbuffer.scale(2.0 * a_exp, {6790, 7000});

                pbuffer.scale(2.0 * a_exp, {7420, 7735});

                pbuffer.scale(pfactors, 0, 2.0, {8365, 8393});

                pbuffer.scale(pfactors, 0, 2.0, {8421, 8505});

                pbuffer.scale(pfactors, 0, 2.0, {8589, 8757});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {8925, 9205});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {9485, 9905});

                pbuffer.scale(pfactors, 0, 2.0, {10325, 10361});

                pbuffer.scale(pfactors, 0, 2.0, {10361, 10469});

                pbuffer.scale(pfactors, 0, 2.0, {10469, 10685});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {10685, 11045});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {11045, 11585});

                t2cfunc::reduce(cbuffer, 3170, pbuffer, 2310, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3180, pbuffer, 2360, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3210, pbuffer, 2510, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3270, pbuffer, 2810, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3370, pbuffer, 3310, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3520, pbuffer, 4060, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3535, pbuffer, 4120, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3580, pbuffer, 4300, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3670, pbuffer, 4660, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3820, pbuffer, 5260, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4045, pbuffer, 6160, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4066, pbuffer, 6223, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4129, pbuffer, 6412, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4255, pbuffer, 6790, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4465, pbuffer, 7420, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4780, pbuffer, 8365, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4808, pbuffer, 8421, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4892, pbuffer, 8589, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5060, pbuffer, 8925, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5340, pbuffer, 9485, 420, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5760, pbuffer, 10325, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5796, pbuffer, 10361, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5904, pbuffer, 10469, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6120, pbuffer, 10685, 360, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6480, pbuffer, 11045, 540, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 600, cbuffer, 0, 10, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 720, cbuffer, 10, 40, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 1620, 600, 720, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 3600, cbuffer, 100, 115, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3780, cbuffer, 115, 160, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 5130, 3600, 3780, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 8010, cbuffer, 250, 271, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 8262, cbuffer, 271, 334, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 10152, 8010, 8262, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 13020, cbuffer, 2070, 2080, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 13140, cbuffer, 2080, 2110, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 14040, 13020, 13140, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 16020, cbuffer, 2170, 2185, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 16200, cbuffer, 2185, 2230, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 17550, 16020, 16200, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 20430, cbuffer, 2320, 2341, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 20682, cbuffer, 2341, 2404, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 22572, 20430, 20682, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 26520, cbuffer, 2530, 2558, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 26856, cbuffer, 2558, 2642, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 29376, 26520, 26856, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 34560, cbuffer, 2810, 2846, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 34992, cbuffer, 2846, 2954, cfactors, 6, 0, 7);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 38232, 34560, 34992, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 460, 470, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 30, cbuffer, 470, 500, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 120, cbuffer, 500, 560, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 300, cbuffer, 560, 660, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 630, cbuffer, 0, 0, 30, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 810, cbuffer, 10, 30, 120, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1080, cbuffer, 40, 120, 300, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1680, 600, 630, 810, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1860, 720, 810, 1080, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 2400, 1620, 1680, 1860, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2700, cbuffer, 810, 825, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2745, cbuffer, 825, 870, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2880, cbuffer, 870, 960, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3150, cbuffer, 960, 1110, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 3645, cbuffer, 100, 2700, 2745, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3915, cbuffer, 115, 2745, 2880, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4320, cbuffer, 160, 2880, 3150, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 5220, 3600, 3645, 3915, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 5490, 3780, 3915, 4320, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 6300, 5130, 5220, 5490, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 6750, cbuffer, 1335, 1356, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 6813, cbuffer, 1356, 1419, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 7002, cbuffer, 1419, 1545, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 7380, cbuffer, 1545, 1755, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 8073, cbuffer, 250, 6750, 6813, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 8451, cbuffer, 271, 6813, 7002, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 9018, cbuffer, 334, 7002, 7380, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 10278, 8010, 8073, 8451, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 10656, 8262, 8451, 9018, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 11790, 10152, 10278, 10656, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 12420, cbuffer, 3170, 3180, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 12450, cbuffer, 3180, 3210, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 12540, cbuffer, 3210, 3270, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 12720, cbuffer, 3270, 3370, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 13050, cbuffer, 2070, 12420, 12450, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 13230, cbuffer, 2080, 12450, 12540, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 13500, cbuffer, 2110, 12540, 12720, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 14100, 13020, 13050, 13230, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 14280, 13140, 13230, 13500, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 14820, 14040, 14100, 14280, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 15120, cbuffer, 3520, 3535, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 15165, cbuffer, 3535, 3580, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 15300, cbuffer, 3580, 3670, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 15570, cbuffer, 3670, 3820, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 16065, cbuffer, 2170, 15120, 15165, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 16335, cbuffer, 2185, 15165, 15300, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 16740, cbuffer, 2230, 15300, 15570, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 17640, 16020, 16065, 16335, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 17910, 16200, 16335, 16740, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 18720, 17550, 17640, 17910, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 19170, cbuffer, 4045, 4066, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 19233, cbuffer, 4066, 4129, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 19422, cbuffer, 4129, 4255, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 19800, cbuffer, 4255, 4465, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 20493, cbuffer, 2320, 19170, 19233, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 20871, cbuffer, 2341, 19233, 19422, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 21438, cbuffer, 2404, 19422, 19800, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 22698, 20430, 20493, 20871, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 23076, 20682, 20871, 21438, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 24210, 22572, 22698, 23076, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 24840, cbuffer, 4780, 4808, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 24924, cbuffer, 4808, 4892, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 25176, cbuffer, 4892, 5060, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 25680, cbuffer, 5060, 5340, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 26604, cbuffer, 2530, 24840, 24924, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 27108, cbuffer, 2558, 24924, 25176, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 27864, cbuffer, 2642, 25176, 25680, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 29544, 26520, 26604, 27108, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 30048, 26856, 27108, 27864, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 31560, 29376, 29544, 30048, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 32400, cbuffer, 5760, 5796, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 32508, cbuffer, 5796, 5904, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 32832, cbuffer, 5904, 6120, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 33480, cbuffer, 6120, 6480, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 34668, cbuffer, 2810, 32400, 32508, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 35316, cbuffer, 2846, 32508, 32832, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 36288, cbuffer, 2954, 32832, 33480, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 38448, 34560, 34668, 35316, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 39096, 34992, 35316, 36288, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 41040, 38232, 38448, 39096, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<3, 0>(skbuffer, 0, ckbuffer, 2400, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 70, ckbuffer, 2500, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 140, ckbuffer, 2600, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 840, ckbuffer, 6300, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 945, ckbuffer, 6450, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1050, ckbuffer, 6600, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 2100, ckbuffer, 11790, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 2247, ckbuffer, 12000, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 2394, ckbuffer, 12210, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 32907, ckbuffer, 14820, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 32977, ckbuffer, 14920, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 33047, ckbuffer, 15020, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 33117, ckbuffer, 18720, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 33222, ckbuffer, 18870, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 33327, ckbuffer, 19020, 0, 4);

            t4cfunc::ket_transform<3, 0>(skbuffer, 33432, ckbuffer, 24210, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 33579, ckbuffer, 24420, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 33726, ckbuffer, 24630, 0, 5);

            t4cfunc::ket_transform<3, 0>(skbuffer, 33873, ckbuffer, 31560, 0, 6);

            t4cfunc::ket_transform<3, 0>(skbuffer, 34069, ckbuffer, 31840, 0, 6);

            t4cfunc::ket_transform<3, 0>(skbuffer, 34265, ckbuffer, 32120, 0, 6);

            t4cfunc::ket_transform<3, 0>(skbuffer, 34461, ckbuffer, 41040, 0, 7);

            t4cfunc::ket_transform<3, 0>(skbuffer, 34713, ckbuffer, 41400, 0, 7);

            t4cfunc::ket_transform<3, 0>(skbuffer, 34965, ckbuffer, 41760, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 5628, 0, 840, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 5838, 70, 945, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 6048, 140, 1050, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 8148, 840, 2100, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 8463, 945, 2247, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 8778, 1050, 2394, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 15897, 5628, 8148, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 16317, 5838, 8463, r_ab, 3, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 16737, 6048, 8778, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 210, 32907, 33117, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 1155, 33117, 33432, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 2541, 33432, 33873, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sixx(skbuffer, 3864, 33873, 34461, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 6258, 0, 210, 1155, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 9093, 840, 1155, 2541, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_phxx(skbuffer, 11928, 2100, 2541, 3864, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 17157, 5628, 6258, 9093, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dgxx(skbuffer, 20937, 8148, 9093, 11928, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ffxx(skbuffer, 26607, 15897, 17157, 20937, r_ab, 3, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 26607, 3, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 343, skbuffer, 27307, 3, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 686, skbuffer, 28007, 3, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1029, skbuffer, 28707, 3, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1372, skbuffer, 29407, 3, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1715, skbuffer, 30107, 3, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2058, skbuffer, 30807, 3, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2401, skbuffer, 31507, 3, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2744, skbuffer, 32207, 3, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 3, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFFFS_hpp */
