#ifndef ElectronRepulsionGeom1010RecSSSD_hpp
#define ElectronRepulsionGeom1010RecSSSD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom1010ContrRecSSXX.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(SS|1/|r-r'||SD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_sssd(T& distributor,
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

    CSimdArray<double> pbuffer(103, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(64, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(72, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(105, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(45, 1);

    // setup Boys fuction data

    const CBoysFunc<4> bf_table;

    CSimdArray<double> bf_data(6, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 5, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 5, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 5, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 5);
                }

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 5, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 8, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 17, 0, 1, 5, 8, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 23, 1, 2, 8, 11, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 29, 2, 3, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 35, 5, 8, 17, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 45, 8, 11, 23, 29, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 55, 8, 17, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 73, 23, 35, 45, pfactors, 26, r_pb, a_exp, b_exp);

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {17, 23});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {35, 45});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {55, 73});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {73, 103});

                t2cfunc::reduce(cbuffer, 0, pbuffer, 17, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 35, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 55, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 34, pbuffer, 73, 30, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 0, cbuffer, 0, 6, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 18, cbuffer, 16, 34, cfactors, 6, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 45, ckbuffer, 0, 0, 0);

            t4cfunc::ket_transform<0, 2>(skbuffer, 50, ckbuffer, 6, 0, 0);

            t4cfunc::ket_transform<0, 2>(skbuffer, 55, ckbuffer, 12, 0, 0);

            t4cfunc::ket_transform<0, 2>(skbuffer, 60, ckbuffer, 18, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 75, ckbuffer, 36, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 90, ckbuffer, 54, 0, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 0, 45, 60, r_ab, 0, 2);

            t4cfunc::bra_transform<0, 0>(sbuffer, 0, skbuffer, 0, 0, 2);

            t4cfunc::bra_transform<0, 0>(sbuffer, 5, skbuffer, 5, 0, 2);

            t4cfunc::bra_transform<0, 0>(sbuffer, 10, skbuffer, 10, 0, 2);

            t4cfunc::bra_transform<0, 0>(sbuffer, 15, skbuffer, 15, 0, 2);

            t4cfunc::bra_transform<0, 0>(sbuffer, 20, skbuffer, 20, 0, 2);

            t4cfunc::bra_transform<0, 0>(sbuffer, 25, skbuffer, 25, 0, 2);

            t4cfunc::bra_transform<0, 0>(sbuffer, 30, skbuffer, 30, 0, 2);

            t4cfunc::bra_transform<0, 0>(sbuffer, 35, skbuffer, 35, 0, 2);

            t4cfunc::bra_transform<0, 0>(sbuffer, 40, skbuffer, 40, 0, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 0, 0, 0, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecSSSD_hpp */
