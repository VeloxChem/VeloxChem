#ifndef ElectronRepulsionGeom1100RecDPSS_hpp
#define ElectronRepulsionGeom1100RecDPSS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DP|1/|r-r'||SS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dpss(T& distributor,
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

    CSimdArray<double> pfactors(23, ket_npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(126, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(92, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(878, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(135, 1);

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

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // loop over basis function pairs on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            // zero integral buffers

            cbuffer.zero();

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

                t4cfunc::comp_distances_wp(pfactors, 20, 17, r_p);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 6, 0, 1, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 9, 1, 2, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 12, 2, 3, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 15, 3, 4, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 18, 4, 5, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 21, 0, 1, 6, 9, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 27, 1, 2, 9, 12, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 33, 2, 3, 12, 15, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 39, 3, 4, 15, 18, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 45, 6, 9, 21, 27, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 55, 9, 12, 27, 33, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 65, 12, 15, 33, 39, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 75, 21, 27, 45, 55, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 90, 27, 33, 55, 65, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 105, 45, 55, 75, 90, pfactors, 20, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1, pbuffer, 6, 3, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {21, 27});

                pbuffer.scale(2.0 * b_exp, {45, 55});

                t2cfunc::reduce(cbuffer, 4, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 45, 10, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {6, 9});

                pbuffer.scale(a_exp / b_exp, {21, 27});

                pbuffer.scale(a_exp / b_exp, {45, 55});

                t2cfunc::reduce(cbuffer, 20, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 21, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 24, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 45, 10, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {21, 27});

                pbuffer.scale(2.0 * b_exp, {45, 55});

                pbuffer.scale(4.0 * a_exp * b_exp, {75, 90});

                pbuffer.scale(4.0 * a_exp * b_exp, {105, 126});

                t2cfunc::reduce(cbuffer, 40, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 46, pbuffer, 45, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 56, pbuffer, 75, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 71, pbuffer, 105, 21, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 0>(skbuffer, 0, cbuffer, 0, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1, cbuffer, 1, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 661, cbuffer, 4, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 667, cbuffer, 10, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 677, cbuffer, 20, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 678, cbuffer, 21, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 690, cbuffer, 24, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 714, cbuffer, 30, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 826, cbuffer, 40, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 832, cbuffer, 46, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 842, cbuffer, 56, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 857, cbuffer, 71, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 799, 678, 690, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 808, 690, 714, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 229, 1, 799, 808, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 4, 0, 661, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 40, 1, 667, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 202, 1, 4, 40, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 681, 677, 826, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 696, 678, 832, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 724, 690, 842, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 754, 714, 857, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 13, 678, 681, 696, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 58, 690, 696, 724, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 112, 714, 724, 754, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 256, 4, 799, 13, 58, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 337, 40, 808, 58, 112, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 499, 202, 229, 256, 337, r_ab, 0, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 499, 0, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 15, skbuffer, 517, 0, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 30, skbuffer, 535, 0, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 45, skbuffer, 553, 0, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 60, skbuffer, 571, 0, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 75, skbuffer, 589, 0, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 90, skbuffer, 607, 0, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 105, skbuffer, 625, 0, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 120, skbuffer, 643, 0, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 0, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDPSS_hpp */
