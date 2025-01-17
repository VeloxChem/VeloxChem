#ifndef ElectronRepulsionGeom1100RecDDSD_hpp
#define ElectronRepulsionGeom1100RecDDSD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DD|1/|r-r'||SD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_ddsd(T& distributor,
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

    CSimdArray<double> pbuffer(1718, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(852, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(7685, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1125, 1);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_data(10, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 9, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 9);

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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 9, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 12, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 21, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 24, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 27, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 30, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 33, 0, 1, 9, 12, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 39, 1, 2, 12, 15, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 45, 2, 3, 15, 18, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 51, 3, 4, 18, 21, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 57, 4, 5, 21, 24, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 63, 5, 6, 24, 27, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 69, 6, 7, 27, 30, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 75, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 78, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 81, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 84, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 87, 2, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 96, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 105, 4, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 114, 5, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 123, 6, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 132, 12, 33, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 150, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 168, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 186, 21, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 204, 24, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 222, 27, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 240, 2, 3, 75, 78, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 246, 3, 4, 78, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 252, 4, 5, 81, 84, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 258, 12, 15, 75, 87, 96, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 276, 15, 18, 78, 96, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 294, 18, 21, 81, 105, 114, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 312, 21, 24, 84, 114, 123, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 330, 33, 39, 87, 132, 150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 366, 39, 45, 96, 150, 168, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 402, 45, 51, 105, 168, 186, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 438, 51, 57, 114, 186, 204, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 474, 57, 63, 123, 204, 222, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 510, 75, 78, 240, 246, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 520, 78, 81, 246, 252, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 530, 87, 96, 240, 258, 276, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 560, 96, 105, 246, 276, 294, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 590, 105, 114, 252, 294, 312, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 620, 132, 150, 258, 330, 366, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 680, 150, 168, 276, 366, 402, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 740, 168, 186, 294, 402, 438, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 800, 186, 204, 312, 438, 474, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 860, 240, 246, 510, 520, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 875, 258, 276, 510, 530, 560, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 920, 276, 294, 520, 560, 590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 965, 330, 366, 530, 620, 680, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1055, 366, 402, 560, 680, 740, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1145, 402, 438, 590, 740, 800, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1235, 530, 560, 860, 875, 920, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 1298, 620, 680, 875, 965, 1055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 1424, 680, 740, 920, 1055, 1145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 1550, 965, 1055, 1235, 1298, 1424, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 132, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 330, 36, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {620, 680});

                pbuffer.scale(2.0 * b_exp, {965, 1055});

                t2cfunc::reduce(cbuffer, 54, pbuffer, 620, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 114, pbuffer, 965, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {132, 150});

                pbuffer.scale(2.0 * a_exp, {330, 366});

                pbuffer.scale(a_exp / b_exp, {620, 680});

                pbuffer.scale(a_exp / b_exp, {965, 1055});

                t2cfunc::reduce(cbuffer, 204, pbuffer, 132, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 222, pbuffer, 330, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 258, pbuffer, 620, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 318, pbuffer, 965, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {620, 680});

                pbuffer.scale(2.0 * b_exp, {965, 1055});

                pbuffer.scale(4.0 * a_exp * b_exp, {1298, 1424});

                pbuffer.scale(4.0 * a_exp * b_exp, {1550, 1718});

                t2cfunc::reduce(cbuffer, 408, pbuffer, 620, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 468, pbuffer, 965, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 558, pbuffer, 1298, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 684, pbuffer, 1550, 168, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 2>(skbuffer, 0, cbuffer, 0, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 15, cbuffer, 18, 0, 2);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6000, cbuffer, 54, 0, 3);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6050, cbuffer, 114, 0, 4);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6125, cbuffer, 204, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6140, cbuffer, 222, 0, 2);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6260, cbuffer, 258, 0, 3);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6460, cbuffer, 318, 0, 4);

            t4cfunc::ket_transform<0, 2>(skbuffer, 7315, cbuffer, 408, 0, 3);

            t4cfunc::ket_transform<0, 2>(skbuffer, 7365, cbuffer, 468, 0, 4);

            t4cfunc::ket_transform<0, 2>(skbuffer, 7440, cbuffer, 558, 0, 5);

            t4cfunc::ket_transform<0, 2>(skbuffer, 7545, cbuffer, 684, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 7075, 6140, 6260, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 7165, 6260, 6460, r_ab, 0, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 1950, 15, 7075, 7165, r_ab, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 45, 0, 6000, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 405, 15, 6050, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 1680, 15, 45, 405, r_ab, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 6170, 6125, 7315, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 6310, 6140, 7365, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 6535, 6260, 7440, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 6760, 6460, 7545, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 135, 6140, 6170, 6310, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 555, 6260, 6310, 6535, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 1005, 6460, 6535, 6760, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 2220, 45, 7075, 135, 555, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 3030, 405, 7165, 555, 1005, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 4380, 1680, 1950, 2220, 3030, r_ab, 0, 2);

            t4cfunc::bra_transform<2, 2>(sbuffer, 0, skbuffer, 4380, 0, 2);

            t4cfunc::bra_transform<2, 2>(sbuffer, 125, skbuffer, 4560, 0, 2);

            t4cfunc::bra_transform<2, 2>(sbuffer, 250, skbuffer, 4740, 0, 2);

            t4cfunc::bra_transform<2, 2>(sbuffer, 375, skbuffer, 4920, 0, 2);

            t4cfunc::bra_transform<2, 2>(sbuffer, 500, skbuffer, 5100, 0, 2);

            t4cfunc::bra_transform<2, 2>(sbuffer, 625, skbuffer, 5280, 0, 2);

            t4cfunc::bra_transform<2, 2>(sbuffer, 750, skbuffer, 5460, 0, 2);

            t4cfunc::bra_transform<2, 2>(sbuffer, 875, skbuffer, 5640, 0, 2);

            t4cfunc::bra_transform<2, 2>(sbuffer, 1000, skbuffer, 5820, 0, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 2, 0, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDDSD_hpp */
