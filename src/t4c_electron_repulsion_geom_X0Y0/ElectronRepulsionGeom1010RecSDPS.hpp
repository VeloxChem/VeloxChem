#ifndef ElectronRepulsionGeom1010RecSDPS_hpp
#define ElectronRepulsionGeom1010RecSDPS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(SD|1/|r-r'||PS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_sdps(T& distributor,
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

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(355, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(176, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(336, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(306, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(135, 1);

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

        cfactors.replicate_points(c_coords, ket_range, 0, 1);

        cfactors.replicate_points(d_coords, ket_range, 3, 1);

        t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        ckbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // loop over basis function pairs on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            // zero integral buffers

            cbuffer.zero();

            ckbuffer.zero();

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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 6, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 9, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 12, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 15, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 18, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 21, 0, 1, 6, 9, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 27, 1, 2, 9, 12, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 33, 2, 3, 12, 15, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 39, 3, 4, 15, 18, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 45, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 48, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 51, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 54, 1, 6, 9, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 63, 2, 9, 12, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 72, 3, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 81, 9, 21, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 99, 12, 27, 33, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 117, 15, 33, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 135, 0, 1, 45, 48, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 141, 1, 2, 48, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 147, 6, 9, 48, 54, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 165, 9, 12, 51, 63, 72, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 183, 21, 27, 63, 81, 99, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 219, 27, 33, 72, 99, 117, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 255, 45, 48, 135, 141, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 265, 54, 63, 141, 147, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 295, 81, 99, 165, 183, 219, pfactors, 26, r_pb, a_exp, b_exp);

                pbuffer.scale(2.0 * a_exp, {135, 141});

                pbuffer.scale(2.0 * a_exp, {255, 265});

                t2cfunc::reduce(cbuffer, 0, pbuffer, 135, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 255, 10, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 2.0, {135, 141});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {147, 165});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {183, 219});

                pbuffer.scale(pfactors, 0, 2.0, {255, 265});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {265, 295});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {295, 355});

                t2cfunc::reduce(cbuffer, 16, pbuffer, 135, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 22, pbuffer, 147, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 40, pbuffer, 183, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 76, pbuffer, 255, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 86, pbuffer, 265, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 116, pbuffer, 295, 60, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 16, 22, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 18, cbuffer, 22, 40, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 72, cbuffer, 0, 0, 18, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 126, cbuffer, 76, 86, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 156, cbuffer, 86, 116, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 246, cbuffer, 6, 126, 156, cfactors, 6, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 162, ckbuffer, 72, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 180, ckbuffer, 90, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 198, ckbuffer, 108, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 216, ckbuffer, 246, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 246, ckbuffer, 276, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 276, ckbuffer, 306, 0, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 0, 162, 216, r_ab, 1, 0);

            t4cfunc::bra_transform<0, 2>(sbuffer, 0, skbuffer, 0, 1, 0);

            t4cfunc::bra_transform<0, 2>(sbuffer, 15, skbuffer, 18, 1, 0);

            t4cfunc::bra_transform<0, 2>(sbuffer, 30, skbuffer, 36, 1, 0);

            t4cfunc::bra_transform<0, 2>(sbuffer, 45, skbuffer, 54, 1, 0);

            t4cfunc::bra_transform<0, 2>(sbuffer, 60, skbuffer, 72, 1, 0);

            t4cfunc::bra_transform<0, 2>(sbuffer, 75, skbuffer, 90, 1, 0);

            t4cfunc::bra_transform<0, 2>(sbuffer, 90, skbuffer, 108, 1, 0);

            t4cfunc::bra_transform<0, 2>(sbuffer, 105, skbuffer, 126, 1, 0);

            t4cfunc::bra_transform<0, 2>(sbuffer, 120, skbuffer, 144, 1, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 0, 2, 1, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecSDPS_hpp */
