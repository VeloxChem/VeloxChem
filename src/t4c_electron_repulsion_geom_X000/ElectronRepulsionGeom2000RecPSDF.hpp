#ifndef ElectronRepulsionGeom2000RecPSDF_hpp
#define ElectronRepulsionGeom2000RecPSDF_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDPXX.hpp"
#include "ElectronRepulsionContrRecDSXX.hpp"
#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecPSXX.hpp"
#include "ElectronRepulsionContrRecXXDF.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionGeom2000ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSH.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(2)/dA^(2)(PS|1/|r-r'||DF)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_psdf(T& distributor,
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

    CSimdArray<double> pbuffer(1801, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1104, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(3240, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(4305, 1);

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

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 75, 9, 12, 33, 39, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 12, 15, 39, 45, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 15, 18, 45, 51, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 18, 21, 51, 57, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 21, 24, 57, 63, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 24, 27, 63, 69, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 135, 33, 39, 75, 85, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 150, 39, 45, 85, 95, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 165, 45, 51, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 180, 51, 57, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 195, 57, 63, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 210, 75, 85, 135, 150, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 231, 85, 95, 150, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 252, 95, 105, 165, 180, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 273, 105, 115, 180, 195, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 294, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 303, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 321, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 339, 39, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 369, 45, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 399, 51, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 429, 85, 135, 150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 474, 95, 150, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 519, 105, 165, 180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 564, 150, 210, 231, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 627, 165, 231, 252, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 690, 180, 252, 273, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 753, 39, 45, 294, 303, 321, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 789, 75, 85, 303, 339, 369, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 849, 85, 95, 321, 369, 399, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 909, 135, 150, 369, 429, 474, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 999, 150, 165, 399, 474, 519, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1089, 210, 231, 474, 564, 627, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1215, 231, 252, 519, 627, 690, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1341, 339, 369, 753, 789, 849, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1441, 429, 474, 849, 909, 999, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 1591, 564, 627, 999, 1089, 1215, pfactors, 26, r_pb, a_exp, b_exp);

                pbuffer.scale(2.0 * a_exp, {75, 85});

                pbuffer.scale(2.0 * a_exp, {135, 150});

                pbuffer.scale(2.0 * a_exp, {210, 231});

                pbuffer.scale(2.0 * a_exp, {339, 369});

                pbuffer.scale(2.0 * a_exp, {429, 474});

                pbuffer.scale(2.0 * a_exp, {564, 627});

                t2cfunc::reduce(cbuffer, 0, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 25, pbuffer, 210, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 46, pbuffer, 339, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 76, pbuffer, 429, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 121, pbuffer, 564, 63, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {75, 85});

                pbuffer.scale(2.0 * a_exp, {135, 150});

                pbuffer.scale(2.0 * a_exp, {210, 231});

                pbuffer.scale(2.0 * a_exp, {339, 369});

                pbuffer.scale(2.0 * a_exp, {429, 474});

                pbuffer.scale(2.0 * a_exp, {564, 627});

                pbuffer.scale(4.0 * a_exp * a_exp, {789, 849});

                pbuffer.scale(4.0 * a_exp * a_exp, {909, 999});

                pbuffer.scale(4.0 * a_exp * a_exp, {1089, 1215});

                pbuffer.scale(4.0 * a_exp * a_exp, {1341, 1441});

                pbuffer.scale(4.0 * a_exp * a_exp, {1441, 1591});

                pbuffer.scale(4.0 * a_exp * a_exp, {1591, 1801});

                t2cfunc::reduce(cbuffer, 184, pbuffer, 75, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 194, pbuffer, 135, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 209, pbuffer, 210, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 230, pbuffer, 339, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 260, pbuffer, 429, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 305, pbuffer, 564, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 368, pbuffer, 789, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 428, pbuffer, 909, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 518, pbuffer, 1089, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 644, pbuffer, 1341, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 744, pbuffer, 1441, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 894, pbuffer, 1591, 210, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 0, cbuffer, 0, 10, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 30, cbuffer, 10, 25, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 75, 0, 30, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 135, cbuffer, 46, 76, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 225, cbuffer, 76, 121, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 360, 135, 225, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 540, cbuffer, 184, 194, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 570, cbuffer, 194, 209, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 615, 540, 570, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 675, cbuffer, 230, 260, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 765, cbuffer, 260, 305, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 900, 675, 765, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1080, cbuffer, 368, 428, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 1260, cbuffer, 428, 518, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 1530, 1080, 1260, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1890, cbuffer, 644, 744, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 2190, cbuffer, 744, 894, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdf(ckbuffer, 2640, 1890, 2190, cfactors, 6, 0, 3);

            t4cfunc::ket_transform<2, 3>(skbuffer, 1470, ckbuffer, 75, 0, 0);

            t4cfunc::ket_transform<2, 3>(skbuffer, 1505, ckbuffer, 360, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 1715, ckbuffer, 615, 0, 0);

            t4cfunc::ket_transform<2, 3>(skbuffer, 1750, ckbuffer, 900, 0, 1);

            t4cfunc::ket_transform<2, 3>(skbuffer, 1855, ckbuffer, 1530, 0, 2);

            t4cfunc::ket_transform<2, 3>(skbuffer, 2065, ckbuffer, 2640, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 1610, 1470, 1505, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 2415, 1715, 1750, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 2520, 1750, 1855, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2835, 1855, 2065, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 3465, 2415, 2520, r_ab, 2, 3);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 3675, 2520, 2835, r_ab, 2, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 0, 1470, 3465, 0, 2, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 210, 1505, 3675, 1, 2, 3);

            erirec::comp_bra_geom20_hrr_electron_repulsion_psxx(skbuffer, 840, 1610, 0, 210, r_ab, 2, 3);

            t4cfunc::bra_transform<1, 0>(sbuffer, 0, skbuffer, 840, 2, 3);

            t4cfunc::bra_transform<1, 0>(sbuffer, 105, skbuffer, 945, 2, 3);

            t4cfunc::bra_transform<1, 0>(sbuffer, 210, skbuffer, 1050, 2, 3);

            t4cfunc::bra_transform<1, 0>(sbuffer, 315, skbuffer, 1155, 2, 3);

            t4cfunc::bra_transform<1, 0>(sbuffer, 420, skbuffer, 1260, 2, 3);

            t4cfunc::bra_transform<1, 0>(sbuffer, 525, skbuffer, 1365, 2, 3);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 1, 0, 2, 3, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecPSDF_hpp */
