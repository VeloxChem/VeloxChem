#ifndef ElectronRepulsionGeom1100RecFPSS_hpp
#define ElectronRepulsionGeom1100RecFPSS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecFPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FP|1/|r-r'||SS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_fpss(T& distributor,
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

    CSimdArray<double> pbuffer(210, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(156, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(2289, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(189, 1);

    // setup Boys fuction data

    const CBoysFunc<6> bf_table;

    CSimdArray<double> bf_data(8, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 7, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 7, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 7, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 7);
                }

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 6, pfactors, 16, bf_data, 6);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 7, 0, 1, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 10, 1, 2, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 13, 2, 3, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 16, 3, 4, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 19, 4, 5, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 22, 5, 6, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 25, 0, 1, 7, 10, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 31, 1, 2, 10, 13, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 37, 2, 3, 13, 16, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 43, 3, 4, 16, 19, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 49, 4, 5, 19, 22, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 55, 7, 10, 25, 31, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 65, 10, 13, 31, 37, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 75, 13, 16, 37, 43, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 85, 16, 19, 43, 49, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 95, 25, 31, 55, 65, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 110, 31, 37, 65, 75, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 125, 37, 43, 75, 85, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 140, 55, 65, 95, 110, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 161, 65, 75, 110, 125, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 182, 95, 110, 140, 161, pfactors, 20, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4, pbuffer, 25, 6, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {25, 31});

                pbuffer.scale(2.0 * b_exp, {55, 65});

                pbuffer.scale(2.0 * b_exp, {95, 110});

                t2cfunc::reduce(cbuffer, 10, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 26, pbuffer, 95, 15, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {7, 10});

                pbuffer.scale(a_exp / b_exp, {25, 31});

                pbuffer.scale(a_exp / b_exp, {55, 65});

                pbuffer.scale(a_exp / b_exp, {95, 110});

                t2cfunc::reduce(cbuffer, 41, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 42, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 45, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 51, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 61, pbuffer, 95, 15, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {25, 31});

                pbuffer.scale(2.0 * b_exp, {55, 65});

                pbuffer.scale(2.0 * b_exp, {95, 110});

                pbuffer.scale(4.0 * a_exp * b_exp, {140, 161});

                pbuffer.scale(4.0 * a_exp * b_exp, {182, 210});

                t2cfunc::reduce(cbuffer, 76, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 82, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 92, pbuffer, 95, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 107, pbuffer, 140, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 128, pbuffer, 182, 28, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 0>(skbuffer, 0, cbuffer, 0, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1, cbuffer, 1, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 40, cbuffer, 4, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1921, cbuffer, 10, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1927, cbuffer, 16, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1937, cbuffer, 26, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1952, cbuffer, 41, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1953, cbuffer, 42, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1965, cbuffer, 45, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1989, cbuffer, 51, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2029, cbuffer, 61, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2209, cbuffer, 76, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2215, cbuffer, 82, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2225, cbuffer, 92, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2240, cbuffer, 107, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 2261, cbuffer, 128, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 373, 1, 40, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 2152, 1953, 1965, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2161, 1965, 1989, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 2179, 1989, 2029, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 409, 1, 2152, 2161, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 571, 40, 2161, 2179, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dpxx(skbuffer, 1111, 373, 409, 571, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 4, 0, 1921, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 46, 1, 1927, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 118, 40, 1937, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 382, 1, 4, 46, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 517, 40, 46, 118, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dpxx(skbuffer, 1057, 373, 382, 517, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 1956, 1952, 2209, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 1971, 1953, 2215, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 1999, 1965, 2225, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 2044, 1989, 2240, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 2089, 2029, 2261, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 13, 1953, 1956, 1971, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 64, 1965, 1971, 1999, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 148, 1989, 1999, 2044, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 238, 2029, 2044, 2089, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 436, 4, 2152, 13, 64, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 625, 46, 2161, 64, 148, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 787, 118, 2179, 148, 238, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 1165, 382, 409, 436, 625, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 1327, 517, 571, 625, 787, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_fpxx(skbuffer, 1651, 1057, 1111, 1165, 1327, r_ab, 0, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 1651, 0, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 21, skbuffer, 1681, 0, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 42, skbuffer, 1711, 0, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 63, skbuffer, 1741, 0, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 84, skbuffer, 1771, 0, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 105, skbuffer, 1801, 0, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 126, skbuffer, 1831, 0, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 147, skbuffer, 1861, 0, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 168, skbuffer, 1891, 0, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 0, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFPSS_hpp */
