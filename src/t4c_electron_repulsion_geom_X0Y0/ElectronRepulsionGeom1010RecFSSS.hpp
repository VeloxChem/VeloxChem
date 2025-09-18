#ifndef ElectronRepulsionGeom1010RecFSSS_hpp
#define ElectronRepulsionGeom1010RecFSSS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
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
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FS|1/|r-r'||SS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fsss(T& distributor,
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

    CSimdArray<double> pbuffer(281, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(180, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(135, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(945, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(63, 1);

    // setup Boys fuction data

    const CBoysFunc<5> bf_table;

    CSimdArray<double> bf_data(7, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 6, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 6, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 6, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 6);
                }

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 6, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 9, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 12, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 21, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 24, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 27, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 30, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 33, 1, 6, 9, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 42, 2, 9, 12, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 51, 3, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 60, 4, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 69, 0, 1, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 75, 1, 2, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 81, 2, 3, 27, 30, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 87, 6, 9, 24, 33, 42, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 105, 9, 12, 27, 42, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 123, 12, 15, 30, 51, 60, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 141, 21, 24, 69, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 151, 24, 27, 75, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 161, 33, 42, 75, 87, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 191, 42, 51, 81, 105, 123, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 221, 69, 75, 141, 151, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 236, 87, 105, 151, 161, 191, pfactors, 26, r_pb, a_exp, b_exp);

                pbuffer.scale(pfactors, 0, 2.0, {0, 1});

                pbuffer.scale(pfactors, 0, 2.0, {6, 9});

                pbuffer.scale(pfactors, 0, 2.0, {21, 24});

                pbuffer.scale(pfactors, 0, 2.0, {33, 42});

                pbuffer.scale(pfactors, 0, 2.0, {69, 75});

                pbuffer.scale(pfactors, 0, 2.0, {87, 105});

                t2cfunc::reduce(cbuffer, 0, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4, pbuffer, 21, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7, pbuffer, 33, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 69, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 22, pbuffer, 87, 18, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {6, 9});

                pbuffer.scale(2.0 * a_exp, {21, 24});

                pbuffer.scale(2.0 * a_exp, {33, 42});

                pbuffer.scale(2.0 * a_exp, {69, 75});

                pbuffer.scale(2.0 * a_exp, {87, 105});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {141, 151});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {161, 191});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {221, 236});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {236, 281});

                t2cfunc::reduce(cbuffer, 40, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 41, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 44, pbuffer, 21, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 47, pbuffer, 33, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 56, pbuffer, 69, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 62, pbuffer, 87, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 80, pbuffer, 141, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 161, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 120, pbuffer, 221, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 135, pbuffer, 236, 45, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 0, 1, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 3, cbuffer, 4, 7, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 12, cbuffer, 16, 22, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 30, cbuffer, 40, 41, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 33, cbuffer, 44, 47, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 42, cbuffer, 56, 62, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 60, cbuffer, 80, 90, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 90, cbuffer, 120, 135, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 0, ckbuffer, 0, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1, ckbuffer, 1, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2, ckbuffer, 2, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 12, ckbuffer, 3, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 15, ckbuffer, 6, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 18, ckbuffer, 9, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 48, ckbuffer, 12, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 54, ckbuffer, 18, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 60, ckbuffer, 24, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 840, ckbuffer, 30, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 841, ckbuffer, 31, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 842, ckbuffer, 32, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 843, ckbuffer, 33, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 846, ckbuffer, 36, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 849, ckbuffer, 39, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 852, ckbuffer, 42, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 858, ckbuffer, 48, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 864, ckbuffer, 54, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 870, ckbuffer, 60, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 880, ckbuffer, 70, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 890, ckbuffer, 80, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 900, ckbuffer, 90, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 915, ckbuffer, 105, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 930, ckbuffer, 120, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 210, 0, 12, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 213, 1, 15, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 216, 2, 18, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 246, 12, 48, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 255, 15, 54, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 264, 18, 60, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 516, 210, 246, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 522, 213, 255, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 528, 216, 264, r_ab, 0, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 3, 840, 843, r_ab, 0, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 21, 843, 852, r_ab, 0, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 66, 852, 870, r_ab, 0, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 120, 870, 900, r_ab, 0, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 219, 0, 3, 21, r_ab, 0, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 273, 12, 21, 66, r_ab, 0, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 354, 48, 66, 120, r_ab, 0, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 534, 210, 219, 273, r_ab, 0, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 588, 246, 273, 354, r_ab, 0, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fsxx(skbuffer, 750, 516, 534, 588, r_ab, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 750, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 7, skbuffer, 760, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 14, skbuffer, 770, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 21, skbuffer, 780, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 28, skbuffer, 790, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 35, skbuffer, 800, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 42, skbuffer, 810, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 49, skbuffer, 820, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 56, skbuffer, 830, 0, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 0, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFSSS_hpp */
