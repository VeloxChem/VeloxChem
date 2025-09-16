#ifndef ElectronRepulsionGeom1100RecDPSP_hpp
#define ElectronRepulsionGeom1100RecDPSP_hpp

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
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DP|1/|r-r'||SP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dpsp(T& distributor,
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

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(450, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(276, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(2634, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(405, 1);

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

                t4cfunc::comp_distances_qd(pfactors, 20, 10, 7);

                t4cfunc::comp_distances_wq(pfactors, 23, 17, 10);

                t4cfunc::comp_distances_wp(pfactors, 26, 17, r_p);

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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 7, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 10, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 25, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 28, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 31, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 34, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 37, 1, 7, 10, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 46, 2, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 55, 3, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 64, 4, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 73, 5, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 82, 1, 2, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 88, 2, 3, 28, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 94, 3, 4, 31, 34, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 100, 7, 10, 25, 37, 46, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 118, 10, 13, 28, 46, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 136, 13, 16, 31, 55, 64, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 154, 16, 19, 34, 64, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 172, 25, 28, 82, 88, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 182, 28, 31, 88, 94, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 192, 37, 46, 82, 100, 118, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 222, 46, 55, 88, 118, 136, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 252, 55, 64, 94, 136, 154, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 282, 82, 88, 172, 182, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 297, 100, 118, 172, 192, 222, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 342, 118, 136, 182, 222, 252, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 387, 192, 222, 282, 297, 342, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 37, 9, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {100, 118});

                pbuffer.scale(2.0 * b_exp, {192, 222});

                t2cfunc::reduce(cbuffer, 12, pbuffer, 100, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 192, 30, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {7, 10});

                pbuffer.scale(2.0 * a_exp, {37, 46});

                pbuffer.scale(a_exp / b_exp, {100, 118});

                pbuffer.scale(a_exp / b_exp, {192, 222});

                t2cfunc::reduce(cbuffer, 60, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 63, pbuffer, 37, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 72, pbuffer, 100, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 192, 30, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {100, 118});

                pbuffer.scale(2.0 * b_exp, {192, 222});

                pbuffer.scale(4.0 * a_exp * b_exp, {297, 342});

                pbuffer.scale(4.0 * a_exp * b_exp, {387, 450});

                t2cfunc::reduce(cbuffer, 120, pbuffer, 100, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 138, pbuffer, 192, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 168, pbuffer, 297, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 213, pbuffer, 387, 63, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 1>(skbuffer, 0, cbuffer, 0, 0, 0);

            t4cfunc::ket_transform<0, 1>(skbuffer, 3, cbuffer, 3, 0, 1);

            t4cfunc::ket_transform<0, 1>(skbuffer, 1983, cbuffer, 12, 0, 2);

            t4cfunc::ket_transform<0, 1>(skbuffer, 2001, cbuffer, 30, 0, 3);

            t4cfunc::ket_transform<0, 1>(skbuffer, 2031, cbuffer, 60, 0, 0);

            t4cfunc::ket_transform<0, 1>(skbuffer, 2034, cbuffer, 63, 0, 1);

            t4cfunc::ket_transform<0, 1>(skbuffer, 2070, cbuffer, 72, 0, 2);

            t4cfunc::ket_transform<0, 1>(skbuffer, 2142, cbuffer, 90, 0, 3);

            t4cfunc::ket_transform<0, 1>(skbuffer, 2478, cbuffer, 120, 0, 2);

            t4cfunc::ket_transform<0, 1>(skbuffer, 2496, cbuffer, 138, 0, 3);

            t4cfunc::ket_transform<0, 1>(skbuffer, 2526, cbuffer, 168, 0, 4);

            t4cfunc::ket_transform<0, 1>(skbuffer, 2571, cbuffer, 213, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 2397, 2034, 2070, r_ab, 0, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2424, 2070, 2142, r_ab, 0, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 687, 3, 2397, 2424, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 12, 0, 1983, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 120, 3, 2001, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 606, 3, 12, 120, r_ab, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 2043, 2031, 2478, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 2088, 2034, 2496, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 2172, 2070, 2526, 0, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 2262, 2142, 2571, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 39, 2034, 2043, 2088, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 174, 2070, 2088, 2172, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 336, 2142, 2172, 2262, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 768, 12, 2397, 39, 174, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 1011, 120, 2424, 174, 336, r_ab, 0, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 1497, 606, 687, 768, 1011, r_ab, 0, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 1497, 0, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 45, skbuffer, 1551, 0, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 90, skbuffer, 1605, 0, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 135, skbuffer, 1659, 0, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 180, skbuffer, 1713, 0, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 225, skbuffer, 1767, 0, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 270, skbuffer, 1821, 0, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 315, skbuffer, 1875, 0, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 360, skbuffer, 1929, 0, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 0, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDPSP_hpp */
