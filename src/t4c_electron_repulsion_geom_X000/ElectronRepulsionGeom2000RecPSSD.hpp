#ifndef ElectronRepulsionGeom2000RecPSSD_hpp
#define ElectronRepulsionGeom2000RecPSSD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDPXX.hpp"
#include "ElectronRepulsionContrRecDSXX.hpp"
#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
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

/// @brief Computes d^(2)/dA^(2)(PS|1/|r-r'||SD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_pssd(T& distributor,
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

    CSimdArray<double> pbuffer(270, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(144, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(615, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(90, 1);

    // setup Boys fuction data

    const CBoysFunc<5> bf_table;

    CSimdArray<double> bf_data(7, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 6, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 6);

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 6, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 9, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 12, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 21, 0, 1, 6, 9, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 27, 1, 2, 9, 12, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 33, 2, 3, 12, 15, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 39, 3, 4, 15, 18, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 45, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 48, 2, 9, 12, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 57, 3, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 66, 9, 21, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 84, 12, 27, 33, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 102, 15, 33, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 120, 9, 12, 45, 48, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 138, 21, 27, 48, 66, 84, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 174, 27, 33, 57, 84, 102, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 210, 66, 84, 120, 138, 174, pfactors, 26, r_pb, a_exp, b_exp);

                pbuffer.scale(2.0 * a_exp, {21, 27});

                pbuffer.scale(2.0 * a_exp, {66, 84});

                t2cfunc::reduce(cbuffer, 0, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 66, 18, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {21, 27});

                pbuffer.scale(2.0 * a_exp, {66, 84});

                pbuffer.scale(4.0 * a_exp * a_exp, {138, 174});

                pbuffer.scale(4.0 * a_exp * a_exp, {210, 270});

                t2cfunc::reduce(cbuffer, 24, pbuffer, 21, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 66, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 48, pbuffer, 138, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 84, pbuffer, 210, 60, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 2>(skbuffer, 210, cbuffer, 0, 0, 0);

            t4cfunc::ket_transform<0, 2>(skbuffer, 215, cbuffer, 6, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 245, cbuffer, 24, 0, 0);

            t4cfunc::ket_transform<0, 2>(skbuffer, 250, cbuffer, 30, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 265, cbuffer, 48, 0, 2);

            t4cfunc::ket_transform<0, 2>(skbuffer, 295, cbuffer, 84, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 230, 210, 215, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 345, 245, 250, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 360, 250, 265, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 405, 265, 295, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 495, 345, 360, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 525, 360, 405, r_ab, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 0, 210, 495, 0, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 30, 215, 525, 1, 0, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_psxx(skbuffer, 120, 230, 0, 30, r_ab, 0, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 0, skbuffer, 120, 0, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 15, skbuffer, 135, 0, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 30, skbuffer, 150, 0, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 45, skbuffer, 165, 0, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 60, skbuffer, 180, 0, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 75, skbuffer, 195, 0, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 1, 0, 0, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecPSSD_hpp */
