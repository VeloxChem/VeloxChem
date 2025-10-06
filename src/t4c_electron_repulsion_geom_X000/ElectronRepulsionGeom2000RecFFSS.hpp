#ifndef ElectronRepulsionGeom2000RecFFSS_hpp
#define ElectronRepulsionGeom2000RecFFSS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDFXX.hpp"
#include "ElectronRepulsionContrRecDGXX.hpp"
#include "ElectronRepulsionContrRecDHXX.hpp"
#include "ElectronRepulsionContrRecDIXX.hpp"
#include "ElectronRepulsionContrRecPFXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecPIXX.hpp"
#include "ElectronRepulsionContrRecPKXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecDGXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecFFXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecPHXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"
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

/// @brief Computes d^(2)/dA^(2)(FF|1/|r-r'||SS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_ffss(T& distributor,
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

    CSimdArray<double> cbuffer(254, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(4373, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(294, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 135, 15, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {75, 85});

                pbuffer.scale(2.0 * a_exp, {135, 150});

                pbuffer.scale(2.0 * a_exp, {210, 231});

                pbuffer.scale(2.0 * a_exp, {294, 322});

                t2cfunc::reduce(cbuffer, 25, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 35, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 50, pbuffer, 210, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 71, pbuffer, 294, 28, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {75, 85});

                pbuffer.scale(2.0 * a_exp, {135, 150});

                pbuffer.scale(2.0 * a_exp, {210, 231});

                pbuffer.scale(2.0 * a_exp, {294, 322});

                pbuffer.scale(4.0 * a_exp * a_exp, {378, 414});

                pbuffer.scale(4.0 * a_exp * a_exp, {450, 495});

                t2cfunc::reduce(cbuffer, 99, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 109, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 124, pbuffer, 210, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 145, pbuffer, 294, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 173, pbuffer, 378, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 209, pbuffer, 450, 45, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 0>(skbuffer, 0, cbuffer, 0, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 70, cbuffer, 10, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3232, cbuffer, 25, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3242, cbuffer, 35, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3257, cbuffer, 50, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3278, cbuffer, 71, 0, 6);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3444, cbuffer, 99, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3454, cbuffer, 109, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3469, cbuffer, 124, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3490, cbuffer, 145, 0, 6);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3518, cbuffer, 173, 0, 7);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3554, cbuffer, 209, 0, 8);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 469, 0, 70, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 3306, 3232, 3242, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 3336, 3242, 3257, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 3381, 3257, 3278, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 3599, 3444, 3454, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 3629, 3454, 3469, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_phxx(skbuffer, 3674, 3469, 3490, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pixx(skbuffer, 3737, 3490, 3518, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pkxx(skbuffer, 3821, 3518, 3554, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 3929, 3599, 3629, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_dgxx(skbuffer, 3989, 3629, 3674, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_dhxx(skbuffer, 4079, 3674, 3737, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_dixx(skbuffer, 4205, 3737, 3821, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pfxx(skbuffer, 499, 0, 3306, 3336, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pgxx(skbuffer, 769, 70, 3336, 3381, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dfxx(skbuffer, 1552, 469, 499, 769, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 10, 3232, 3929, 3, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 85, 3242, 3989, 4, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 175, 3257, 4079, 5, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 301, 3278, 4205, 6, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_pfxx(skbuffer, 589, 3306, 10, 85, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_pgxx(skbuffer, 904, 3336, 85, 175, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_phxx(skbuffer, 1174, 3381, 175, 301, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dfxx(skbuffer, 1732, 499, 589, 904, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_dgxx(skbuffer, 2092, 769, 904, 1174, r_ab, 0, 0);

            erirec::comp_bra_geom20_hrr_electron_repulsion_ffxx(skbuffer, 2632, 1552, 1732, 2092, r_ab, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 2632, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 49, skbuffer, 2732, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 98, skbuffer, 2832, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 147, skbuffer, 2932, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 196, skbuffer, 3032, 0, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 245, skbuffer, 3132, 0, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 0, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecFFSS_hpp */
