#ifndef ElectronRepulsionGeom1100RecFSSD_hpp
#define ElectronRepulsionGeom1100RecFSSD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPSXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecFSXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(FS|1/|r-r'||SD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_fssd(T& distributor,
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

    CSimdArray<double> pbuffer(1011, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(588, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(5620, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(315, 1);

    // setup Boys fuction data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 8, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 8);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 8, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 29, 0, 1, 8, 11, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 35, 1, 2, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 41, 2, 3, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 47, 3, 4, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 4, 5, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 5, 6, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 65, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 68, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 71, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 74, 2, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 83, 3, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 92, 4, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 101, 5, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 110, 11, 29, 35, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 128, 14, 35, 41, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 146, 17, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 164, 20, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 182, 23, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 200, 2, 3, 65, 68, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 206, 3, 4, 68, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 212, 11, 14, 65, 74, 83, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 230, 14, 17, 68, 83, 92, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 248, 17, 20, 71, 92, 101, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 266, 29, 35, 74, 110, 128, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 302, 35, 41, 83, 128, 146, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 338, 41, 47, 92, 146, 164, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 374, 47, 53, 101, 164, 182, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 410, 65, 68, 200, 206, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 420, 74, 83, 200, 212, 230, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 450, 83, 92, 206, 230, 248, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 480, 110, 128, 212, 266, 302, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 540, 128, 146, 230, 302, 338, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 600, 146, 164, 248, 338, 374, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 660, 212, 230, 410, 420, 450, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 705, 266, 302, 420, 480, 540, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 795, 302, 338, 450, 540, 600, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 885, 480, 540, 660, 705, 795, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 110, 18, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {110, 128});

                pbuffer.scale(2.0 * b_exp, {266, 302});

                pbuffer.scale(2.0 * b_exp, {480, 540});

                t2cfunc::reduce(cbuffer, 24, pbuffer, 110, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 42, pbuffer, 266, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 78, pbuffer, 480, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {29, 35});

                pbuffer.scale(a_exp / b_exp, {110, 128});

                pbuffer.scale(a_exp / b_exp, {266, 302});

                pbuffer.scale(a_exp / b_exp, {480, 540});

                t2cfunc::reduce(cbuffer, 138, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 144, pbuffer, 110, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 162, pbuffer, 266, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 198, pbuffer, 480, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {110, 128});

                pbuffer.scale(2.0 * b_exp, {266, 302});

                pbuffer.scale(2.0 * b_exp, {480, 540});

                pbuffer.scale(4.0 * a_exp * b_exp, {705, 795});

                pbuffer.scale(4.0 * a_exp * b_exp, {885, 1011});

                t2cfunc::reduce(cbuffer, 258, pbuffer, 110, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 276, pbuffer, 266, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 312, pbuffer, 480, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 372, pbuffer, 705, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 462, pbuffer, 885, 126, ket_width, ket_npgtos);

            }

            t4cfunc::ket_transform<0, 2>(skbuffer, 0, cbuffer, 0, 0, 0);

            t4cfunc::ket_transform<0, 2>(skbuffer, 5, cbuffer, 6, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 4445, cbuffer, 24, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 4460, cbuffer, 42, 0, 2);

            t4cfunc::ket_transform<0, 2>(skbuffer, 4490, cbuffer, 78, 0, 3);

            t4cfunc::ket_transform<0, 2>(skbuffer, 4540, cbuffer, 138, 0, 0);

            t4cfunc::ket_transform<0, 2>(skbuffer, 4545, cbuffer, 144, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 4605, cbuffer, 162, 0, 2);

            t4cfunc::ket_transform<0, 2>(skbuffer, 4725, cbuffer, 198, 0, 3);

            t4cfunc::ket_transform<0, 2>(skbuffer, 5300, cbuffer, 258, 0, 1);

            t4cfunc::ket_transform<0, 2>(skbuffer, 5315, cbuffer, 276, 0, 2);

            t4cfunc::ket_transform<0, 2>(skbuffer, 5345, cbuffer, 312, 0, 3);

            t4cfunc::ket_transform<0, 2>(skbuffer, 5395, cbuffer, 372, 0, 4);

            t4cfunc::ket_transform<0, 2>(skbuffer, 5470, cbuffer, 462, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 1010, 0, 5, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 5150, 4540, 4545, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 5165, 4545, 4605, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 5210, 4605, 4725, r_ab, 0, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 5575, 5300, 5315, r_ab, 0, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_psxx(skbuffer, 1070, 0, 5150, 5165, r_ab, 0, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 1385, 5, 5165, 5210, r_ab, 0, 2);

            erirec::comp_bra_geom10_hrr_electron_repulsion_dsxx(skbuffer, 2825, 1010, 1070, 1385, r_ab, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 20, 0, 4460, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 200, 5, 4490, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_psxx(skbuffer, 1025, 0, 4445, 20, r_ab, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 1250, 5, 20, 200, r_ab, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_dsxx(skbuffer, 2735, 1010, 1025, 1250, r_ab, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 4560, 4540, 5315, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 4635, 4545, 5345, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 4775, 4605, 5395, 0, 2);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 4925, 4725, 5470, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 65, 4545, 4560, 4635, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 290, 4605, 4635, 4775, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 560, 4725, 4775, 4925, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_psxx(skbuffer, 1115, 4445, 5150, 5575, 65, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 1520, 20, 5165, 65, 290, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 1925, 200, 5210, 290, 560, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dsxx(skbuffer, 2915, 1025, 1070, 1115, 1520, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 3185, 1250, 1385, 1520, 1925, r_ab, 0, 2);

            erirec::comp_bra_geom11_hrr_electron_repulsion_fsxx(skbuffer, 3995, 2735, 2825, 2915, 3185, r_ab, 0, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 3995, 0, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 35, skbuffer, 4045, 0, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 70, skbuffer, 4095, 0, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 105, skbuffer, 4145, 0, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 140, skbuffer, 4195, 0, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 175, skbuffer, 4245, 0, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 210, skbuffer, 4295, 0, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 245, skbuffer, 4345, 0, 2);

            t4cfunc::bra_transform<3, 0>(sbuffer, 280, skbuffer, 4395, 0, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 0, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecFSSD_hpp */
