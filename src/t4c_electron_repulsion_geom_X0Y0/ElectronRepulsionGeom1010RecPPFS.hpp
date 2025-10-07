#ifndef ElectronRepulsionGeom1010RecPPFS_hpp
#define ElectronRepulsionGeom1010RecPPFS_hpp

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
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(PP|1/|r-r'||FS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_ppfs(T& distributor,
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

    CSimdArray<double> pbuffer(1260, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(990, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(5940, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(1596, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(567, 1);

    // setup Boys fuction data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 8, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 8);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 8, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 29, 0, 1, 8, 11, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 35, 1, 2, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 41, 2, 3, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 47, 3, 4, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 4, 5, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 5, 6, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 65, 8, 11, 29, 35, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 75, 11, 14, 35, 41, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 14, 17, 41, 47, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 17, 20, 47, 53, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 20, 23, 53, 59, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 115, 29, 35, 65, 75, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 130, 35, 41, 75, 85, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 145, 41, 47, 85, 95, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 160, 47, 53, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 175, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 178, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 181, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 184, 1, 8, 11, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 193, 2, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 202, 3, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 211, 11, 29, 35, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 229, 14, 35, 41, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 247, 17, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 265, 35, 65, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 295, 41, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 325, 47, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 355, 75, 115, 130, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 400, 85, 130, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 445, 95, 145, 160, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 490, 0, 1, 175, 178, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 496, 1, 2, 178, 181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 502, 8, 11, 178, 184, 193, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 520, 11, 14, 181, 193, 202, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 538, 29, 35, 193, 211, 229, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 574, 35, 41, 202, 229, 247, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 610, 65, 75, 229, 265, 295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 670, 75, 85, 247, 295, 325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 730, 115, 130, 295, 355, 400, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 820, 130, 145, 325, 400, 445, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 910, 175, 178, 490, 496, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 920, 184, 193, 496, 502, 520, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 950, 211, 229, 520, 538, 574, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1010, 265, 295, 574, 610, 670, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1110, 355, 400, 670, 730, 820, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 175, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 184, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12, pbuffer, 211, 18, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {175, 178});

                pbuffer.scale(2.0 * a_exp, {184, 193});

                pbuffer.scale(2.0 * a_exp, {211, 229});

                pbuffer.scale(2.0 * a_exp, {490, 496});

                pbuffer.scale(2.0 * a_exp, {502, 520});

                pbuffer.scale(2.0 * a_exp, {538, 574});

                pbuffer.scale(2.0 * a_exp, {910, 920});

                pbuffer.scale(2.0 * a_exp, {920, 950});

                pbuffer.scale(2.0 * a_exp, {950, 1010});

                t2cfunc::reduce(cbuffer, 135, pbuffer, 175, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 138, pbuffer, 184, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 147, pbuffer, 211, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 165, pbuffer, 490, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 171, pbuffer, 502, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 189, pbuffer, 538, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 225, pbuffer, 910, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 235, pbuffer, 920, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 265, pbuffer, 950, 60, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {175, 178});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {184, 193});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {211, 229});

                pbuffer.scale(pfactors, 0, 2.0, {265, 295});

                pbuffer.scale(pfactors, 0, 2.0, {355, 400});

                t2cfunc::reduce(cbuffer, 30, pbuffer, 175, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 33, pbuffer, 184, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 42, pbuffer, 211, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 265, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 355, 45, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {175, 178});

                pbuffer.scale(2.0 * a_exp, {184, 193});

                pbuffer.scale(2.0 * a_exp, {211, 229});

                pbuffer.scale(2.0 * a_exp, {265, 295});

                pbuffer.scale(2.0 * a_exp, {355, 400});

                pbuffer.scale(pfactors, 0, 2.0, {490, 496});

                pbuffer.scale(pfactors, 0, 2.0, {502, 520});

                pbuffer.scale(pfactors, 0, 2.0, {538, 574});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {610, 670});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {730, 820});

                pbuffer.scale(pfactors, 0, 2.0, {910, 920});

                pbuffer.scale(pfactors, 0, 2.0, {920, 950});

                pbuffer.scale(pfactors, 0, 2.0, {950, 1010});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1010, 1110});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1110, 1260});

                t2cfunc::reduce(cbuffer, 325, pbuffer, 175, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 328, pbuffer, 184, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 337, pbuffer, 211, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 355, pbuffer, 265, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 385, pbuffer, 355, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 430, pbuffer, 490, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 436, pbuffer, 502, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 454, pbuffer, 538, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 490, pbuffer, 610, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 550, pbuffer, 730, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 640, pbuffer, 910, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 650, pbuffer, 920, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 680, pbuffer, 950, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 740, pbuffer, 1010, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 840, pbuffer, 1110, 150, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 180, cbuffer, 0, 3, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 216, cbuffer, 3, 12, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 486, 180, 216, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 990, cbuffer, 135, 138, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1026, cbuffer, 138, 147, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 1296, 990, 1026, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1980, cbuffer, 165, 171, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2052, cbuffer, 171, 189, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 2592, 1980, 2052, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 3840, cbuffer, 225, 235, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 3960, cbuffer, 235, 265, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 4860, 3840, 3960, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 30, 33, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9, cbuffer, 33, 42, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 36, cbuffer, 42, 60, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 90, cbuffer, 60, 90, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 189, cbuffer, 0, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 243, cbuffer, 3, 9, 36, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 324, cbuffer, 12, 36, 90, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 504, 180, 189, 243, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 558, 216, 243, 324, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 720, 486, 504, 558, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 810, cbuffer, 325, 328, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 819, cbuffer, 328, 337, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 846, cbuffer, 337, 355, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 900, cbuffer, 355, 385, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 999, cbuffer, 135, 810, 819, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1053, cbuffer, 138, 819, 846, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1134, cbuffer, 147, 846, 900, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1314, 990, 999, 1053, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 1368, 1026, 1053, 1134, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 1530, 1296, 1314, 1368, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1620, cbuffer, 430, 436, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1638, cbuffer, 436, 454, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1692, cbuffer, 454, 490, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1800, cbuffer, 490, 550, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1998, cbuffer, 165, 1620, 1638, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2106, cbuffer, 171, 1638, 1692, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2268, cbuffer, 189, 1692, 1800, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 2628, 1980, 1998, 2106, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 2736, 2052, 2106, 2268, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 3060, 2592, 2628, 2736, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 3240, cbuffer, 640, 650, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3270, cbuffer, 650, 680, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3360, cbuffer, 680, 740, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3540, cbuffer, 740, 840, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 3870, cbuffer, 225, 3240, 3270, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4050, cbuffer, 235, 3270, 3360, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4320, cbuffer, 265, 3360, 3540, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 4920, 3840, 3870, 4050, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 5100, 3960, 4050, 4320, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 5640, 4860, 4920, 5100, cfactors, 6, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 0, ckbuffer, 720, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 21, ckbuffer, 750, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 42, ckbuffer, 780, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1197, ckbuffer, 1530, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1218, ckbuffer, 1560, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1239, ckbuffer, 1590, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1260, ckbuffer, 3060, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1302, ckbuffer, 3120, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1344, ckbuffer, 3180, 0, 2);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1386, ckbuffer, 5640, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1456, ckbuffer, 5740, 0, 3);

            t4cfunc::ket_transform<3, 0>(skbuffer, 1526, ckbuffer, 5840, 0, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 63, 1197, 1260, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 252, 1260, 1386, r_ab, 3, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 630, 0, 63, 252, r_ab, 3, 0);

            t4cfunc::bra_transform<1, 1>(sbuffer, 0, skbuffer, 630, 3, 0);

            t4cfunc::bra_transform<1, 1>(sbuffer, 63, skbuffer, 693, 3, 0);

            t4cfunc::bra_transform<1, 1>(sbuffer, 126, skbuffer, 756, 3, 0);

            t4cfunc::bra_transform<1, 1>(sbuffer, 189, skbuffer, 819, 3, 0);

            t4cfunc::bra_transform<1, 1>(sbuffer, 252, skbuffer, 882, 3, 0);

            t4cfunc::bra_transform<1, 1>(sbuffer, 315, skbuffer, 945, 3, 0);

            t4cfunc::bra_transform<1, 1>(sbuffer, 378, skbuffer, 1008, 3, 0);

            t4cfunc::bra_transform<1, 1>(sbuffer, 441, skbuffer, 1071, 3, 0);

            t4cfunc::bra_transform<1, 1>(sbuffer, 504, skbuffer, 1134, 3, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 1, 1, 3, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecPPFS_hpp */
