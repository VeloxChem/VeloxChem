#ifndef ElectronRepulsionGeom1100RecPSPD_hpp
#define ElectronRepulsionGeom1100RecPSPD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPSXX.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(PS|1/|r-r'||PD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_pspd(T& distributor,
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

    CSimdArray<double> pbuffer(630, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(416, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(954, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(1785, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(405, 1);

    // setup Boys fuction data

    const CBoysFunc<6> bf_table;

    CSimdArray<double> bf_data(8, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 7, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 7);

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

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 25, 0, 1, 7, 10, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 31, 1, 2, 10, 13, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 37, 2, 3, 13, 16, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 43, 3, 4, 16, 19, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 49, 4, 5, 19, 22, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 55, 7, 10, 25, 31, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 65, 10, 13, 31, 37, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 75, 13, 16, 37, 43, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 16, 19, 43, 49, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 95, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 98, 2, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 107, 3, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 116, 10, 25, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 134, 13, 31, 37, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 152, 16, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 170, 31, 55, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 200, 37, 65, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 230, 43, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 260, 10, 13, 95, 98, 107, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 278, 25, 31, 98, 116, 134, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 314, 31, 37, 107, 134, 152, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 350, 55, 65, 134, 170, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 410, 65, 75, 152, 200, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 470, 116, 134, 260, 278, 314, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 530, 170, 200, 314, 350, 410, pfactors, 26, r_pb, a_exp, b_exp);

                pbuffer.scale(2.0 * b_exp, {116, 134});

                pbuffer.scale(2.0 * b_exp, {170, 200});

                t2cfunc::reduce(cbuffer, 0, pbuffer, 116, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 170, 30, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {25, 31});

                pbuffer.scale(2.0 * a_exp, {55, 65});

                pbuffer.scale(a_exp / b_exp, {116, 134});

                pbuffer.scale(a_exp / b_exp, {170, 200});

                t2cfunc::reduce(cbuffer, 48, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 54, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 64, pbuffer, 116, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 82, pbuffer, 170, 30, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {116, 134});

                pbuffer.scale(2.0 * b_exp, {170, 200});

                pbuffer.scale(4.0 * a_exp * b_exp, {278, 314});

                pbuffer.scale(4.0 * a_exp * b_exp, {350, 410});

                pbuffer.scale(4.0 * a_exp * b_exp, {470, 530});

                pbuffer.scale(4.0 * a_exp * b_exp, {530, 630});

                t2cfunc::reduce(cbuffer, 112, pbuffer, 116, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 130, pbuffer, 170, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 160, pbuffer, 278, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 196, pbuffer, 350, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 256, pbuffer, 470, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 316, pbuffer, 530, 100, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 18, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 54, cbuffer, 48, 54, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 72, cbuffer, 64, 82, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 612, cbuffer, 112, 130, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 666, cbuffer, 160, 196, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 774, cbuffer, 256, 316, cfactors, 6, 0, 3);

            t4cfunc::ket_transform<1, 2>(skbuffer, 810, ckbuffer, 0, 0, 1);

            t4cfunc::ket_transform<1, 2>(skbuffer, 855, ckbuffer, 54, 0, 0);

            t4cfunc::ket_transform<1, 2>(skbuffer, 870, ckbuffer, 72, 0, 1);

            t4cfunc::ket_transform<1, 2>(skbuffer, 1365, ckbuffer, 612, 0, 1);

            t4cfunc::ket_transform<1, 2>(skbuffer, 1410, ckbuffer, 666, 0, 2);

            t4cfunc::ket_transform<1, 2>(skbuffer, 1500, ckbuffer, 774, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 1320, 855, 870, r_ab, 1, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1650, 1365, 1410, r_ab, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 915, 855, 1410, 1, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 1050, 870, 1500, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 0, 870, 915, 1050, r_ab, 1, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_psxx(skbuffer, 405, 810, 1320, 1650, 0, r_ab, 1, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 0, skbuffer, 405, 1, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 45, skbuffer, 450, 1, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 90, skbuffer, 495, 1, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 135, skbuffer, 540, 1, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 180, skbuffer, 585, 1, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 225, skbuffer, 630, 1, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 270, skbuffer, 675, 1, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 315, skbuffer, 720, 1, 2);

            t4cfunc::bra_transform<1, 0>(sbuffer, 360, skbuffer, 765, 1, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 1, 0, 1, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecPSPD_hpp */
