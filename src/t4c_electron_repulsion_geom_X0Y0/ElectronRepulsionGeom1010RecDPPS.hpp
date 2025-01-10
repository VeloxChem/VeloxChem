#ifndef ElectronRepulsionGeom1010RecDPPS_hpp
#define ElectronRepulsionGeom1010RecDPPS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXPS.hpp"
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DP|1/|r-r'||PS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dpps(T& distributor,
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

    CSimdArray<double> pbuffer(705, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(473, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(903, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(2196, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(405, 1);

    // setup Boys fuction data

    const CBoysFunc<6> bf_table;

    CSimdArray<double> bf_data(8, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 7, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 7);

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 6, pfactors, 16, bf_data, 6);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 7, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 10, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 25, 0, 1, 7, 10, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 31, 1, 2, 10, 13, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 37, 2, 3, 13, 16, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 43, 3, 4, 16, 19, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 49, 4, 5, 19, 22, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 55, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 58, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 61, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 64, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 67, 1, 7, 10, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 76, 2, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 85, 3, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 94, 4, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 103, 10, 25, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 121, 13, 31, 37, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 139, 16, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 157, 19, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 175, 0, 1, 55, 58, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 181, 1, 2, 58, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 187, 2, 3, 61, 64, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 193, 7, 10, 58, 67, 76, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 211, 10, 13, 61, 76, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 229, 13, 16, 64, 85, 94, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 247, 25, 31, 76, 103, 121, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 283, 31, 37, 85, 121, 139, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 319, 37, 43, 94, 139, 157, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 355, 55, 58, 175, 181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 365, 58, 61, 181, 187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 375, 67, 76, 181, 193, 211, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 405, 76, 85, 187, 211, 229, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 435, 103, 121, 211, 247, 283, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 495, 121, 139, 229, 283, 319, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 555, 175, 181, 355, 365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 570, 193, 211, 365, 375, 405, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 615, 247, 283, 405, 435, 495, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 55, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 175, 6, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {55, 58});

                pbuffer.scale(2.0 * a_exp, {175, 181});

                pbuffer.scale(2.0 * a_exp, {355, 365});

                pbuffer.scale(2.0 * a_exp, {555, 570});

                t2cfunc::reduce(cbuffer, 99, pbuffer, 55, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 102, pbuffer, 175, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 108, pbuffer, 355, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 118, pbuffer, 555, 15, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {55, 58});

                pbuffer.scale(pfactors, 0, 2.0, {67, 76});

                pbuffer.scale(pfactors, 0, 2.0, {103, 121});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {175, 181});

                pbuffer.scale(pfactors, 0, 2.0, {193, 211});

                pbuffer.scale(pfactors, 0, 2.0, {247, 283});

                t2cfunc::reduce(cbuffer, 9, pbuffer, 55, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12, pbuffer, 67, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 21, pbuffer, 103, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 39, pbuffer, 175, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 45, pbuffer, 193, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 63, pbuffer, 247, 36, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {55, 58});

                pbuffer.scale(2.0 * a_exp, {67, 76});

                pbuffer.scale(2.0 * a_exp, {103, 121});

                pbuffer.scale(2.0 * a_exp, {175, 181});

                pbuffer.scale(2.0 * a_exp, {193, 211});

                pbuffer.scale(2.0 * a_exp, {247, 283});

                pbuffer.scale(pfactors, 0, 2.0, {355, 365});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {375, 405});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {435, 495});

                pbuffer.scale(pfactors, 0, 2.0, {555, 570});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {570, 615});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {615, 705});

                t2cfunc::reduce(cbuffer, 133, pbuffer, 55, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 136, pbuffer, 67, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 145, pbuffer, 103, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 163, pbuffer, 175, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 169, pbuffer, 193, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 187, pbuffer, 247, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 223, pbuffer, 355, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 233, pbuffer, 375, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 263, pbuffer, 435, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 323, pbuffer, 555, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 338, pbuffer, 570, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 383, pbuffer, 615, 90, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 9, 12, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9, cbuffer, 12, 21, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 36, cbuffer, 0, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 63, cbuffer, 39, 45, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 81, cbuffer, 45, 63, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 135, cbuffer, 3, 63, 81, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 189, cbuffer, 133, 136, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 198, cbuffer, 136, 145, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 225, cbuffer, 99, 189, 198, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 252, cbuffer, 163, 169, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 270, cbuffer, 169, 187, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 324, cbuffer, 102, 252, 270, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 378, cbuffer, 223, 233, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 408, cbuffer, 233, 263, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 498, cbuffer, 108, 378, 408, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 588, cbuffer, 323, 338, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 633, cbuffer, 338, 383, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 768, cbuffer, 118, 588, 633, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 0, ckbuffer, 36, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 9, ckbuffer, 45, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 18, ckbuffer, 54, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 108, ckbuffer, 135, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 126, ckbuffer, 153, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 144, ckbuffer, 171, 0, 2);

            //t4cfunc::ket_transform<1, 0>(skbuffer, 594, ckbuffer, 0, 1, 1);

            //t4cfunc::ket_transform<1, 0>(skbuffer, 621, ckbuffer, 27, 1, 1);

            //t4cfunc::ket_transform<1, 0>(skbuffer, 648, ckbuffer, 54, 1, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 1890, ckbuffer, 225, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 1899, ckbuffer, 234, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 1908, ckbuffer, 243, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 1917, ckbuffer, 324, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 1935, ckbuffer, 342, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 1953, ckbuffer, 360, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 1971, ckbuffer, 498, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2001, ckbuffer, 528, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2031, ckbuffer, 558, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2061, ckbuffer, 768, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2106, ckbuffer, 813, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2151, ckbuffer, 858, 0, 4);
            
            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 594, 0, 108, r_ab, 1, 0);
            
            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 621, 9, 126, r_ab, 1, 0);
            
            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 648, 18, 144, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 27, 1890, 1917, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 162, 1917, 1971, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 324, 1971, 2061, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 675, 0, 27, 162, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 918, 108, 162, 324, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 1404, 594, 675, 918, r_ab, 1, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 1404, 1, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 45, skbuffer, 1458, 1, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 90, skbuffer, 1512, 1, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 135, skbuffer, 1566, 1, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 180, skbuffer, 1620, 1, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 225, skbuffer, 1674, 1, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 270, skbuffer, 1728, 1, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 315, skbuffer, 1782, 1, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 360, skbuffer, 1836, 1, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 1, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDPPS_hpp */
