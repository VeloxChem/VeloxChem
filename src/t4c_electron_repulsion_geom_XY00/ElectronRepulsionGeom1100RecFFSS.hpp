#ifndef ElectronRepulsionGeom1100RecFFSS_hpp
#define ElectronRepulsionGeom1100RecFFSS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSIXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSKXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecFFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPHXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSIXX.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSKSS.hpp"
#include "ElectronRepulsionPrimRecSLSS.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FF|1/|r-r'||SS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_ffss(T& distributor,
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

    CSimdArray<double> pbuffer(495, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(320, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(5924, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(441, 1);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 9, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 9);
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

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 8, pfactors, 16, bf_data, 8);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 9, 0, 1, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 12, 1, 2, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 15, 2, 3, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 18, 3, 4, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 21, 4, 5, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 24, 5, 6, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 27, 6, 7, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 30, 7, 8, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 33, 0, 1, 9, 12, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 39, 1, 2, 12, 15, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 45, 2, 3, 15, 18, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 51, 3, 4, 18, 21, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 57, 4, 5, 21, 24, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 63, 5, 6, 24, 27, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 69, 6, 7, 27, 30, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 75, 9, 12, 33, 39, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 85, 12, 15, 39, 45, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 95, 15, 18, 45, 51, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 105, 18, 21, 51, 57, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 115, 21, 24, 57, 63, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 125, 24, 27, 63, 69, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 135, 33, 39, 75, 85, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 150, 39, 45, 85, 95, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 165, 45, 51, 95, 105, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 180, 51, 57, 105, 115, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 195, 57, 63, 115, 125, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 210, 75, 85, 135, 150, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 231, 85, 95, 150, 165, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 252, 95, 105, 165, 180, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 273, 105, 115, 180, 195, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 294, 135, 150, 210, 231, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 322, 150, 165, 231, 252, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 350, 165, 180, 252, 273, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_skss(pbuffer, 378, 210, 231, 294, 322, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_skss(pbuffer, 414, 231, 252, 322, 350, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_slss(pbuffer, 450, 294, 322, 378, 414, pfactors, 20, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 135, 15, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {135, 150});

                pbuffer.scale(2.0 * b_exp, {210, 231});

                pbuffer.scale(2.0 * b_exp, {294, 322});

                t2cfunc::reduce(cbuffer, 31, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 46, pbuffer, 210, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 67, pbuffer, 294, 28, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {33, 39});

                pbuffer.scale(2.0 * a_exp, {75, 85});

                pbuffer.scale(a_exp / b_exp, {135, 150});

                pbuffer.scale(a_exp / b_exp, {210, 231});

                pbuffer.scale(a_exp / b_exp, {294, 322});

                t2cfunc::reduce(cbuffer, 95, pbuffer, 33, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 101, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 111, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 126, pbuffer, 210, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 147, pbuffer, 294, 28, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {135, 150});

                pbuffer.scale(2.0 * b_exp, {210, 231});

                pbuffer.scale(2.0 * b_exp, {294, 322});

                pbuffer.scale(4.0 * a_exp * b_exp, {378, 414});

                pbuffer.scale(4.0 * a_exp * b_exp, {450, 495});

                t2cfunc::reduce(cbuffer, 175, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 190, pbuffer, 210, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 211, pbuffer, 294, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 239, pbuffer, 378, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 275, pbuffer, 450, 45, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 0>(skbuffer, 0, cbuffer, 0, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 6, cbuffer, 6, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 136, cbuffer, 16, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5167, cbuffer, 31, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5182, cbuffer, 46, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5203, cbuffer, 67, 0, 6);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5231, cbuffer, 95, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5237, cbuffer, 101, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5277, cbuffer, 111, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5337, cbuffer, 126, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5421, cbuffer, 147, 0, 6);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5779, cbuffer, 175, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5794, cbuffer, 190, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5815, cbuffer, 211, 0, 6);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5843, cbuffer, 239, 0, 7);

            t4cfunc::ket_transform<0, 0>(skbuffer, 5879, cbuffer, 275, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 835, 6, 136, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 5641, 5237, 5277, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 5671, 5277, 5337, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 5716, 5337, 5421, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 955, 6, 5641, 5671, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 1450, 136, 5671, 5716, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 2737, 835, 955, 1450, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 16, 0, 5167, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 151, 6, 5182, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 331, 136, 5203, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pfxx(skbuffer, 865, 6, 16, 151, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pgxx(skbuffer, 1315, 136, 151, 331, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dfxx(skbuffer, 2557, 835, 865, 1315, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 5247, 5231, 5779, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 5292, 5237, 5794, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 5358, 5277, 5815, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sixx(skbuffer, 5449, 5337, 5843, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_skxx(skbuffer, 5533, 5421, 5879, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 46, 5237, 5247, 5292, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 196, 5277, 5292, 5358, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_shxx(skbuffer, 394, 5337, 5358, 5449, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sixx(skbuffer, 583, 5421, 5449, 5533, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 1045, 16, 5641, 46, 196, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pgxx(skbuffer, 1585, 151, 5671, 196, 394, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_phxx(skbuffer, 1990, 331, 5716, 394, 583, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dfxx(skbuffer, 2917, 865, 955, 1045, 1585, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dgxx(skbuffer, 3457, 1315, 1450, 1585, 1990, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ffxx(skbuffer, 4267, 2557, 2737, 2917, 3457, r_ab, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 4267, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 49, skbuffer, 4367, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 98, skbuffer, 4467, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 147, skbuffer, 4567, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 196, skbuffer, 4667, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 245, skbuffer, 4767, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 294, skbuffer, 4867, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 343, skbuffer, 4967, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 392, skbuffer, 5067, 0, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 0, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFFSS_hpp */
