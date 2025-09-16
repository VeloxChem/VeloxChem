#ifndef ElectronRepulsionGeom1100RecDFSS_hpp
#define ElectronRepulsionGeom1100RecDFSS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSIXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSHXX.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSKSS.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DF|1/|r-r'||SS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dfss(T& distributor,
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

    CSimdArray<double> pbuffer(330, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(204, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(2385, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(315, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 8, 0, 1, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 11, 1, 2, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 14, 2, 3, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 17, 3, 4, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 20, 4, 5, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 23, 5, 6, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 26, 6, 7, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 29, 0, 1, 8, 11, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 35, 1, 2, 11, 14, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 41, 2, 3, 14, 17, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 47, 3, 4, 17, 20, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 53, 4, 5, 20, 23, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 59, 5, 6, 23, 26, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 65, 8, 11, 29, 35, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 75, 11, 14, 35, 41, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 85, 14, 17, 41, 47, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 95, 17, 20, 47, 53, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 105, 20, 23, 53, 59, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 115, 29, 35, 65, 75, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 130, 35, 41, 75, 85, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 145, 41, 47, 85, 95, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 160, 47, 53, 95, 105, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 175, 65, 75, 115, 130, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 196, 75, 85, 130, 145, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 217, 85, 95, 145, 160, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 238, 115, 130, 175, 196, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 266, 130, 145, 196, 217, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_skss(pbuffer, 294, 175, 196, 238, 266, pfactors, 20, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 65, 10, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {115, 130});

                pbuffer.scale(2.0 * b_exp, {175, 196});

                t2cfunc::reduce(cbuffer, 16, pbuffer, 115, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 31, pbuffer, 175, 21, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {29, 35});

                pbuffer.scale(2.0 * a_exp, {65, 75});

                pbuffer.scale(a_exp / b_exp, {115, 130});

                pbuffer.scale(a_exp / b_exp, {175, 196});

                t2cfunc::reduce(cbuffer, 52, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 58, pbuffer, 65, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 68, pbuffer, 115, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 83, pbuffer, 175, 21, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {115, 130});

                pbuffer.scale(2.0 * b_exp, {175, 196});

                pbuffer.scale(4.0 * a_exp * b_exp, {238, 266});

                pbuffer.scale(4.0 * a_exp * b_exp, {294, 330});

                t2cfunc::reduce(cbuffer, 104, pbuffer, 115, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 119, pbuffer, 175, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 140, pbuffer, 238, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 168, pbuffer, 294, 36, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 0>(skbuffer, 0, cbuffer, 0, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 6, cbuffer, 6, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1900, cbuffer, 16, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1915, cbuffer, 31, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1936, cbuffer, 52, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1942, cbuffer, 58, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1982, cbuffer, 68, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2042, cbuffer, 83, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2285, cbuffer, 104, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2300, cbuffer, 119, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2321, cbuffer, 140, 0, 6);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2349, cbuffer, 168, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 2210, 1942, 1982, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 2240, 1982, 2042, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 595, 6, 2210, 2240, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 16, 0, 1900, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 136, 6, 1915, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 505, 6, 16, 136, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 1952, 1936, 2285, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 1997, 1942, 2300, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 2063, 1982, 2321, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 2126, 2042, 2349, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 46, 1942, 1952, 1997, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 181, 1982, 1997, 2063, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 316, 2042, 2063, 2126, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 685, 16, 2210, 46, 181, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 955, 136, 2240, 181, 316, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 1360, 505, 595, 685, 955, r_ab, 0, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 0, skbuffer, 1360, 0, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 35, skbuffer, 1420, 0, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 70, skbuffer, 1480, 0, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 105, skbuffer, 1540, 0, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 140, skbuffer, 1600, 0, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 175, skbuffer, 1660, 0, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 210, skbuffer, 1720, 0, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 245, skbuffer, 1780, 0, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 280, skbuffer, 1840, 0, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 3, 0, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDFSS_hpp */
