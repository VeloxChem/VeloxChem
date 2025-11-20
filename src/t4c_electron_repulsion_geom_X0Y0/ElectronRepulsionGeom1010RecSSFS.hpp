#ifndef ElectronRepulsionGeom1010RecSSFS_hpp
#define ElectronRepulsionGeom1010RecSSFS_hpp

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
#include "ElectronRepulsionGeom1010ContrRecSSXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(SS|1/|r-r'||FS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_ssfs(T& distributor,
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

    CSimdArray<double> pbuffer(210, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(180, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(1080, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(147, 1);

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

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 21, 0, 1, 6, 9, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 27, 1, 2, 9, 12, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 33, 2, 3, 12, 15, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 39, 3, 4, 15, 18, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 45, 6, 9, 21, 27, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 55, 9, 12, 27, 33, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 65, 12, 15, 33, 39, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 75, 21, 27, 45, 55, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 90, 27, 33, 55, 65, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 105, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 108, 1, 6, 9, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 117, 9, 21, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 135, 27, 45, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 165, 55, 75, 90, pfactors, 26, r_pb, a_exp, b_exp);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {6, 9});

                pbuffer.scale(2.0 * a_exp, {21, 27});

                pbuffer.scale(2.0 * a_exp, {105, 108});

                pbuffer.scale(2.0 * a_exp, {108, 117});

                pbuffer.scale(2.0 * a_exp, {117, 135});

                t2cfunc::reduce(cbuffer, 0, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 105, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 13, pbuffer, 108, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 22, pbuffer, 117, 18, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {0, 1});

                pbuffer.scale(pfactors, 0, 2.0, {6, 9});

                pbuffer.scale(pfactors, 0, 2.0, {21, 27});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {45, 55});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {75, 90});

                pbuffer.scale(pfactors, 0, 2.0, {105, 108});

                pbuffer.scale(pfactors, 0, 2.0, {108, 117});

                pbuffer.scale(pfactors, 0, 2.0, {117, 135});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {135, 165});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {165, 210});

                t2cfunc::reduce(cbuffer, 40, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 41, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 44, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 50, pbuffer, 45, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 75, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 75, pbuffer, 105, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 78, pbuffer, 108, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 87, pbuffer, 117, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 105, pbuffer, 135, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 135, pbuffer, 165, 45, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 60, cbuffer, 0, 1, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 72, cbuffer, 1, 4, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 162, 60, 72, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 450, cbuffer, 10, 13, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 486, cbuffer, 13, 22, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxds(ckbuffer, 756, 450, 486, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 40, 41, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3, cbuffer, 41, 44, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 12, cbuffer, 44, 50, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 30, cbuffer, 50, 60, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 63, cbuffer, 0, 0, 3, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 81, cbuffer, 1, 3, 12, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 108, cbuffer, 4, 12, 30, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 168, 60, 63, 81, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 186, 72, 81, 108, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 240, 162, 168, 186, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 270, cbuffer, 75, 78, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 279, cbuffer, 78, 87, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 306, cbuffer, 87, 105, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 360, cbuffer, 105, 135, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 459, cbuffer, 10, 270, 279, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 513, cbuffer, 13, 279, 306, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 594, cbuffer, 22, 306, 360, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 774, 450, 459, 513, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 828, 486, 513, 594, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfs(ckbuffer, 990, 756, 774, 828, cfactors, 6, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 63, ckbuffer, 240, 0, 0);

            t4cfunc::ket_transform<3, 0>(skbuffer, 70, ckbuffer, 250, 0, 0);

            t4cfunc::ket_transform<3, 0>(skbuffer, 77, ckbuffer, 260, 0, 0);

            t4cfunc::ket_transform<3, 0>(skbuffer, 84, ckbuffer, 990, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 105, ckbuffer, 1020, 0, 1);

            t4cfunc::ket_transform<3, 0>(skbuffer, 126, ckbuffer, 1050, 0, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 0, 63, 84, r_ab, 3, 0);

            t4cfunc::bra_transform<0, 0>(sbuffer, 0, skbuffer, 0, 3, 0);

            t4cfunc::bra_transform<0, 0>(sbuffer, 7, skbuffer, 7, 3, 0);

            t4cfunc::bra_transform<0, 0>(sbuffer, 14, skbuffer, 14, 3, 0);

            t4cfunc::bra_transform<0, 0>(sbuffer, 21, skbuffer, 21, 3, 0);

            t4cfunc::bra_transform<0, 0>(sbuffer, 28, skbuffer, 28, 3, 0);

            t4cfunc::bra_transform<0, 0>(sbuffer, 35, skbuffer, 35, 3, 0);

            t4cfunc::bra_transform<0, 0>(sbuffer, 42, skbuffer, 42, 3, 0);

            t4cfunc::bra_transform<0, 0>(sbuffer, 49, skbuffer, 49, 3, 0);

            t4cfunc::bra_transform<0, 0>(sbuffer, 56, skbuffer, 56, 3, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 0, 0, 3, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecSSFS_hpp */
