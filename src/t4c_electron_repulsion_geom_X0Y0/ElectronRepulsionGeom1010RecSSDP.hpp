#ifndef ElectronRepulsionGeom1010RecSSDP_hpp
#define ElectronRepulsionGeom1010RecSSDP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDP.hpp"
#include "ElectronRepulsionGeom1010ContrRecSSXX.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(SS|1/|r-r'||DP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_ssdp(T& distributor,
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

    CSimdArray<double> pbuffer(207, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(172, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(804, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(315, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(135, 1);

    // setup Boys fuction data

    const CBoysFunc<5> bf_table;

    CSimdArray<double> bf_data(7, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 6, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 6);

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

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 21, 0, 1, 6, 9, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 27, 1, 2, 9, 12, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 33, 2, 3, 12, 15, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 39, 3, 4, 15, 18, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 45, 6, 9, 21, 27, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 55, 9, 12, 27, 33, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 65, 12, 15, 33, 39, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 75, 21, 27, 45, 55, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 90, 27, 33, 55, 65, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 105, 1, 6, 9, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 114, 9, 21, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 132, 27, 45, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 162, 55, 75, 90, pfactors, 26, r_pb, a_exp, b_exp);

                pbuffer.scale(2.0 * a_exp, {6, 9});

                pbuffer.scale(2.0 * a_exp, {21, 27});

                pbuffer.scale(2.0 * a_exp, {105, 114});

                pbuffer.scale(2.0 * a_exp, {114, 132});

                t2cfunc::reduce(cbuffer, 0, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 105, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 114, 18, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {6, 9});

                pbuffer.scale(pfactors, 0, 2.0, {21, 27});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {45, 55});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {75, 90});

                pbuffer.scale(pfactors, 0, 2.0, {105, 114});

                pbuffer.scale(pfactors, 0, 2.0, {114, 132});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {132, 162});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {162, 207});

                t2cfunc::reduce(cbuffer, 36, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 39, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 45, pbuffer, 45, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 55, pbuffer, 75, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 70, pbuffer, 105, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 79, pbuffer, 114, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 97, pbuffer, 132, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 127, pbuffer, 162, 45, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 57, cbuffer, 0, 3, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 372, cbuffer, 9, 18, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 36, 39, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 9, cbuffer, 39, 45, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 27, cbuffer, 45, 55, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 66, cbuffer, 0, 0, 9, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 93, cbuffer, 3, 9, 27, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 147, 57, 66, 93, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 201, cbuffer, 70, 79, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 228, cbuffer, 79, 97, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 282, cbuffer, 97, 127, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 399, cbuffer, 9, 201, 228, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 480, cbuffer, 18, 228, 282, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 642, 372, 399, 480, cfactors, 6, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 135, ckbuffer, 147, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 150, ckbuffer, 165, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 165, ckbuffer, 183, 0, 0);

            t4cfunc::ket_transform<2, 1>(skbuffer, 180, ckbuffer, 642, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 225, ckbuffer, 696, 0, 1);

            t4cfunc::ket_transform<2, 1>(skbuffer, 270, ckbuffer, 750, 0, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 0, 135, 180, r_ab, 2, 1);

            t4cfunc::bra_transform<0, 0>(sbuffer, 0, skbuffer, 0, 2, 1);

            t4cfunc::bra_transform<0, 0>(sbuffer, 15, skbuffer, 15, 2, 1);

            t4cfunc::bra_transform<0, 0>(sbuffer, 30, skbuffer, 30, 2, 1);

            t4cfunc::bra_transform<0, 0>(sbuffer, 45, skbuffer, 45, 2, 1);

            t4cfunc::bra_transform<0, 0>(sbuffer, 60, skbuffer, 60, 2, 1);

            t4cfunc::bra_transform<0, 0>(sbuffer, 75, skbuffer, 75, 2, 1);

            t4cfunc::bra_transform<0, 0>(sbuffer, 90, skbuffer, 90, 2, 1);

            t4cfunc::bra_transform<0, 0>(sbuffer, 105, skbuffer, 105, 2, 1);

            t4cfunc::bra_transform<0, 0>(sbuffer, 120, skbuffer, 120, 2, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 0, 0, 2, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecSSDP_hpp */