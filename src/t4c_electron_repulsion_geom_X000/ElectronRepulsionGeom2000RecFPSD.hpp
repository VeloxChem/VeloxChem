#ifndef ElectronRepulsionGeom2000RecFPSD_hpp
#define ElectronRepulsionGeom2000RecFPSD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDDXX.hpp"
#include "ElectronRepulsionContrRecDFXX.hpp"
#include "ElectronRepulsionContrRecDGXX.hpp"
#include "ElectronRepulsionContrRecDPXX.hpp"
#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecFPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"
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

/// @brief Computes d^(2)/dA^(2)(FP|1/|r-r'||SD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_fpsd(T& distributor,
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

    CSimdArray<double> cbuffer(756, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(8730, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(630, 1);

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

                t4cfunc::comp_distances_qd(pfactors, 20, 10, 7);

                t4cfunc::comp_distances_wq(pfactors, 23, 17, 10);

                t4cfunc::comp_distances_wp(pfactors, 26, 17, r_p);

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

                pbuffer.scale(2.0 * a_exp, {132, 150});

                pbuffer.scale(2.0 * a_exp, {330, 366});

                pbuffer.scale(2.0 * a_exp, {620, 680});

                pbuffer.scale(2.0 * a_exp, {965, 1055});

                t2cfunc::reduce(cbuffer, 54, pbuffer, 132, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 72, pbuffer, 330, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 108, pbuffer, 620, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 168, pbuffer, 965, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {132, 150});

                pbuffer.scale(2.0 * a_exp, {330, 366});

                pbuffer.scale(2.0 * a_exp, {620, 680});

                pbuffer.scale(2.0 * a_exp, {965, 1055});

                pbuffer.scale(4.0 * a_exp * a_exp, {1298, 1424});

                pbuffer.scale(4.0 * a_exp * a_exp, {1550, 1718});

                t2cfunc::reduce(cbuffer, 258, pbuffer, 132, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 276, pbuffer, 330, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 312, pbuffer, 620, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 372, pbuffer, 965, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 462, pbuffer, 1298, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 588, pbuffer, 1550, 168, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 2>(skbuffer, 0, cbuffer, 0, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 105, cbuffer, 18, 0, 2);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6015, cbuffer, 54, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6030, cbuffer, 72, 0, 2);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6060, cbuffer, 108, 0, 3);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6110, cbuffer, 168, 0, 4);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6470, cbuffer, 258, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6485, cbuffer, 276, 0, 2);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6515, cbuffer, 312, 0, 3);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6565, cbuffer, 372, 0, 4);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6640, cbuffer, 462, 0, 5);

            t4cfunc::ket_transform<0, 2>(skbuffer, 6745, cbuffer, 588, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1065, 0, 105, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 6185, 6015, 6030, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 6230, 6030, 6060, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 6320, 6060, 6110, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 6885, 6470, 6485, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 6930, 6485, 6515, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 7020, 6515, 6565, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 7170, 6565, 6640, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 7395, 6640, 6745, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 7710, 6885, 6930, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 7800, 6930, 7020, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 7980, 7020, 7170, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_dgxx(skbuffer, 8280, 7170, 7395, r_ab, 0, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 1110, 0, 6185, 6230, r_ab, 0, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 1515, 105, 6230, 6320, r_ab, 0, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dpxx(skbuffer, 3225, 1065, 1110, 1515, r_ab, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 15, 6015, 7710, 1, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 135, 6030, 7800, 2, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 315, 6060, 7980, 3, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 615, 6110, 8280, 4, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_ppxx(skbuffer, 1245, 6185, 15, 135, r_ab, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_pdxx(skbuffer, 1785, 6230, 135, 315, r_ab, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_pfxx(skbuffer, 2325, 6320, 315, 615, r_ab, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dpxx(skbuffer, 3495, 1110, 1245, 1785, r_ab, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_ddxx(skbuffer, 4035, 1515, 1785, 2325, r_ab, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_fpxx(skbuffer, 5115, 3225, 3495, 4035, r_ab, 0, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 5115, 0, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 105, skbuffer, 5265, 0, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 210, skbuffer, 5415, 0, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 315, skbuffer, 5565, 0, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 420, skbuffer, 5715, 0, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 525, skbuffer, 5865, 0, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 0, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecFPSD_hpp */
