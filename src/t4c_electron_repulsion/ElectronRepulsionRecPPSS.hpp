#ifndef ElectronRepulsionRecPPSS_hpp
#define ElectronRepulsionRecPPSS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes (PP|1/|r-r'||SS)  integrals for two basis function pairs blocks.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
/// @param bra_eq_ket True if basis function pairs blocks on bra and ket are the same, False otherwise.
template <class T, int N>
inline auto
comp_electron_repulsion_ppss(T& distributor,
                             const CGtoPairBlock& bra_gto_pair_block,
                             const CGtoPairBlock& ket_gto_pair_block,
                             const std::pair<size_t, size_t>& bra_indices,
                             const std::pair<size_t, size_t>& ket_indices,
                             const bool bra_eq_ket) -> void
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

    if constexpr (N == 1) CSimdArray<double> pbuffer(15, ket_npgtos);

    if constexpr (N == 2) CSimdArray<double> pbuffer(15, ket_npgtos);

    if constexpr (N == 3) CSimdArray<double> pbuffer(30, ket_npgtos);

    // allocate aligned Cartesian integrals

    if constexpr (N == 1) CSimdArray<double> cbuffer(9, 1);

    if constexpr (N == 2) CSimdArray<double> cbuffer(9, 1);

    if constexpr (N == 3) CSimdArray<double> cbuffer(18, 1);

    // allocate aligned half transformed integrals

    if constexpr (N == 1) CSimdArray<double> skbuffer(18, 1);

    if constexpr (N == 2) CSimdArray<double> skbuffer(18, 1);

    if constexpr (N == 3) CSimdArray<double> skbuffer(36, 1);

    // allocate aligned spherical integrals

    if constexpr (N == 1) CSimdArray<double> sbuffer(9, 1);

    if constexpr (N == 2) CSimdArray<double> sbuffer(9, 1);

    if constexpr (N == 3) CSimdArray<double> sbuffer(18, 1);

    // setup Boys fuction data

    const CBoysFunc<2> bf_table;

    if constexpr (N == 1) CSimdArray<double> bf_data(4, ket_npgtos);

    if constexpr (N == 2) CSimdArray<double> bf_data(4, ket_npgtos);

    if constexpr (N == 3) CSimdArray<double> bf_data(8, ket_npgtos);

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

                if constexpr (N == 1) t4cfunc::comp_boys_args(bf_data, 3, pfactors, 13, a_exp, b_exp);

                if constexpr (N == 2) t4cfunc::comp_boys_args(bf_data, 3, pfactors, 13, a_exp, b_exp, omega);

                if constexpr (N == 3)
                {
                    t4cfunc::comp_boys_args(bf_data, 3, pfactors, 13, a_exp, b_exp);

                    t4cfunc::comp_boys_args(bf_data, 7, pfactors, 13, a_exp, b_exp, omega);
                }

                if constexpr (N == 1) bf_table.compute(bf_data, 0, 3);

                if constexpr (N == 2) bf_table.compute(bf_data, 0, 3, pfactors, a_exp, b_exp, omega);

                if constexpr (N == 3)
                {
                    bf_table.compute(bf_data, 0, 3);

                    bf_table.compute(bf_data, 4, 7, pfactors, a_exp, b_exp, omega);
                }

                t4cfunc::comp_ovl_factors(pfactors, 16, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(buffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(buffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(buffer, 2, pfactors, 16, bf_data, 2);

                if constexpr (N == 3)
                {
                    erirec::comp_prim_electron_repulsion_ssss(buffer, 15, pfactors, 16, bf_data, 4);

                    erirec::comp_prim_electron_repulsion_ssss(buffer, 16, pfactors, 16, bf_data, 5);

                    erirec::comp_prim_electron_repulsion_ssss(buffer, 17, pfactors, 16, bf_data, 6);
                }

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 3, 0, 1, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 6, 1, 2, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 9, 0, 1, 3, 6, pfactors, 20, r_pb, a_exp, b_exp);

                if constexpr (N == 3)
                {
                    erirec::comp_prim_electron_repulsion_spss(pbuffer, 18, 15, 16, pfactors, 20, r_pb);

                    erirec::comp_prim_electron_repulsion_spss(pbuffer, 21, 16, 17, pfactors, 20, r_pb);

                    erirec::comp_prim_electron_repulsion_sdss(pbuffer, 24, 15, 16, 18, 21, pfactors, 20, r_pb, a_exp, b_exp);
                }

                t2cfunc::reduce(cbuffer, 0, pbuffer, 3, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 9, 6, ket_width, ket_npgtos);

                if constexpr (N == 3)
                {
                    t2cfunc::reduce(cbuffer, 9, pbuffer, 18, 3, ket_width, ket_npgtos);

                    t2cfunc::reduce(cbuffer, 12, pbuffer, 24, 6, ket_width, ket_npgtos);
                }
            }

            t4cfunc::ket_transform<0, 0>(skbuffer, 0, cbuffer, 0, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3, cbuffer, 3, 0, 2);

            if constexpr (N == 3)
            {
                t4cfunc::ket_transform<0, 0>(skbuffer, 18, cbuffer, 9, 0, 1);

                t4cfunc::ket_transform<0, 0>(skbuffer, 21, cbuffer, 12, 0, 2);
            }

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 9, 0, 3, r_ab, 0, 0);

            if constexpr (N == 3)
            {
                erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 27, 18, 21, r_ab, 0, 0);
            }

            t4cfunc::bra_transform<1, 1>(sbuffer, 0, skbuffer, 9, 0, 0);

            if constexpr (N == 3)
            {
                t4cfunc::bra_transform<1, 1>(sbuffer, 9, skbuffer, 27, 0, 0);
            }
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionRecPPSS_hpp */
