#ifndef ElectronRepulsionGeom1000RecFFSS_hpp
#define ElectronRepulsionGeom1000RecFFSS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDFXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecPIXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDGXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecFFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPHXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)(FF|1/|r-r'||SS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1000_ffss(T& distributor,
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

    CSimdArray<double> cbuffer(156, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(1677, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(147, 1);

    // setup Boys fuction data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 8);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 65, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 115, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 25, pbuffer, 175, 21, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {65, 75});

                pbuffer.scale(2.0 * a_exp, {115, 130});

                pbuffer.scale(2.0 * a_exp, {175, 196});

                pbuffer.scale(2.0 * a_exp, {238, 266});

                pbuffer.scale(2.0 * a_exp, {294, 330});

                t2cfunc::reduce(cbuffer, 46, pbuffer, 65, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 56, pbuffer, 115, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 71, pbuffer, 175, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 92, pbuffer, 238, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 120, pbuffer, 294, 36, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 0>(skbuffer, 0, cbuffer, 0, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 10, cbuffer, 10, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 25, cbuffer, 25, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 181, cbuffer, 46, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 191, cbuffer, 56, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 206, cbuffer, 71, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 227, cbuffer, 92, 0, 6);

            t4cfunc::ket_transform<0, 0>(skbuffer, 255, cbuffer, 120, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 46, 0, 10, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 76, 10, 25, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 121, 46, 76, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 291, 181, 191, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 321, 191, 206, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 366, 206, 227, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pixx(skbuffer, 429, 227, 255, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 513, 0, 291, 321, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 603, 10, 321, 366, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_phxx(skbuffer, 738, 25, 366, 429, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 927, 46, 513, 603, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dgxx(skbuffer, 1107, 76, 603, 738, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ffxx(skbuffer, 1377, 121, 927, 1107, r_ab, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 1377, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 49, skbuffer, 1477, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 98, skbuffer, 1577, 0, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 0, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1000RecFFSS_hpp */