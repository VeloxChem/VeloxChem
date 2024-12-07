#ifndef ElectronRepulsionGeom1100RecDDSS_hpp
#define ElectronRepulsionGeom1100RecDDSS_hpp

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
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DD|1/|r-r'||SS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_ddss(T& distributor,
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

    CSimdArray<double> pbuffer(210, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(142, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(1537, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(225, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 7, 0, 1, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 10, 1, 2, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 13, 2, 3, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 16, 3, 4, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 19, 4, 5, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 22, 5, 6, pfactors, 20, r_pb);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 25, 0, 1, 7, 10, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 31, 1, 2, 10, 13, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 37, 2, 3, 13, 16, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 43, 3, 4, 16, 19, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 49, 4, 5, 19, 22, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 55, 7, 10, 25, 31, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 65, 10, 13, 31, 37, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 75, 13, 16, 37, 43, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 85, 16, 19, 43, 49, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 95, 25, 31, 55, 65, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 110, 31, 37, 65, 75, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 125, 37, 43, 75, 85, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 140, 55, 65, 95, 110, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 161, 65, 75, 110, 125, pfactors, 20, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 182, 95, 110, 140, 161, pfactors, 20, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 25, 6, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {55, 65});

                pbuffer.scale(2.0 * b_exp, {95, 110});

                t2cfunc::reduce(cbuffer, 9, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 19, pbuffer, 95, 15, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {7, 10});

                pbuffer.scale(2.0 * a_exp, {25, 31});

                pbuffer.scale(a_exp / b_exp, {55, 65});

                pbuffer.scale(a_exp / b_exp, {95, 110});

                t2cfunc::reduce(cbuffer, 34, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 37, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 43, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 53, pbuffer, 95, 15, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {55, 65});

                pbuffer.scale(2.0 * b_exp, {95, 110});

                pbuffer.scale(4.0 * a_exp * b_exp, {140, 161});

                pbuffer.scale(4.0 * a_exp * b_exp, {182, 210});

                t2cfunc::reduce(cbuffer, 68, pbuffer, 55, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 78, pbuffer, 95, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 93, pbuffer, 140, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 114, pbuffer, 182, 28, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 0>(skbuffer, 0, cbuffer, 0, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 3, cbuffer, 3, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1200, cbuffer, 9, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1210, cbuffer, 19, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1225, cbuffer, 34, 0, 1);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1228, cbuffer, 37, 0, 2);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1252, cbuffer, 43, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1292, cbuffer, 53, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1463, cbuffer, 68, 0, 3);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1473, cbuffer, 78, 0, 4);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1488, cbuffer, 93, 0, 5);

            t4cfunc::ket_transform<0, 0>(skbuffer, 1509, cbuffer, 114, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 1415, 1228, 1252, r_ab, 0, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 1433, 1252, 1292, r_ab, 0, 0);

            erirec::comp_bra_geom10_hrr_electron_repulsion_pdxx(skbuffer, 390, 3, 1415, 1433, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 9, 0, 1200, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 81, 3, 1210, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_pdxx(skbuffer, 336, 3, 9, 81, r_ab, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 1234, 1225, 1463, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 1262, 1228, 1473, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 1307, 1252, 1488, 0, 0);

            erirec::comp_bra_geom01_hrr_electron_repulsion_shxx(skbuffer, 1352, 1292, 1509, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 27, 1228, 1234, 1262, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 111, 1252, 1262, 1307, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sgxx(skbuffer, 201, 1292, 1307, 1352, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 444, 9, 1415, 27, 111, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pfxx(skbuffer, 606, 81, 1433, 111, 201, r_ab, 0, 0);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ddxx(skbuffer, 876, 336, 390, 444, 606, r_ab, 0, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 0, skbuffer, 876, 0, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 25, skbuffer, 912, 0, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 50, skbuffer, 948, 0, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 75, skbuffer, 984, 0, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 100, skbuffer, 1020, 0, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 125, skbuffer, 1056, 0, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 150, skbuffer, 1092, 0, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 175, skbuffer, 1128, 0, 0);

            t4cfunc::bra_transform<2, 2>(sbuffer, 200, skbuffer, 1164, 0, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 2, 0, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDDSS_hpp */
