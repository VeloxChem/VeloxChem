#ifndef ElectronRepulsionGeom2000RecFSSS_hpp
#define ElectronRepulsionGeom2000RecFSSS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDDXX.hpp"
#include "ElectronRepulsionContrRecDFXX.hpp"
#include "ElectronRepulsionContrRecDPXX.hpp"
#include "ElectronRepulsionContrRecDSXX.hpp"
#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecFSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"
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

/// @brief Computes d^(2)/dA^(2)(FS|1/|r-r'||SS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_fsss(T& distributor,
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

    CSimdArray<double> cbuffer(80, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(896, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(42, 1);

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

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {6, 9});

                pbuffer.scale(2.0 * a_exp, {21, 27});

                pbuffer.scale(2.0 * a_exp, {45, 55});

                t2cfunc::reduce(cbuffer, 4, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 5, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 8, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14, pbuffer, 45, 10, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {6, 9});

                pbuffer.scale(2.0 * a_exp, {21, 27});

                pbuffer.scale(2.0 * a_exp, {45, 55});

                pbuffer.scale(4.0 * a_exp * a_exp, {75, 90});

                pbuffer.scale(4.0 * a_exp * a_exp, {105, 126});

                t2cfunc::reduce(cbuffer, 24, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 25, pbuffer, 6, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 28, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 34, pbuffer, 45, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 44, pbuffer, 75, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 59, pbuffer, 105, 21, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 0>(skbuffer, 0, cbuffer, 0, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 7, cbuffer, 1, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 565, cbuffer, 4, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 566, cbuffer, 5, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 569, cbuffer, 8, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 575, cbuffer, 14, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 615, cbuffer, 24, 0, 0);

            t4cfunc::ket_transform<0, 0>(skbuffer, 616, cbuffer, 25, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 619, cbuffer, 28, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 625, cbuffer, 34, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 635, cbuffer, 44, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 650, cbuffer, 59, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 124, 0, 7, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 585, 565, 566, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 588, 566, 569, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 597, 569, 575, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 671, 615, 616, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 674, 616, 619, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 683, 619, 625, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 701, 625, 635, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 731, 635, 650, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 776, 671, 674, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 782, 674, 683, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 800, 683, 701, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 836, 701, 731, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_psxx(skbuffer, 127, 0, 585, 588, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 154, 7, 588, 597, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dsxx(skbuffer, 343, 124, 127, 154, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 1, 565, 776, 0, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 10, 566, 782, 1, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 28, 569, 800, 2, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 64, 575, 836, 3, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_psxx(skbuffer, 136, 585, 1, 10, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_ppxx(skbuffer, 181, 588, 10, 28, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_pdxx(skbuffer, 235, 597, 28, 64, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dsxx(skbuffer, 361, 127, 136, 181, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dpxx(skbuffer, 397, 154, 181, 235, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_fsxx(skbuffer, 505, 343, 361, 397, r_ab, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 505, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 7, skbuffer, 515, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 14, skbuffer, 525, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 21, skbuffer, 535, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 28, skbuffer, 545, 0, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 35, skbuffer, 555, 0, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 0, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecFSSS_hpp */
