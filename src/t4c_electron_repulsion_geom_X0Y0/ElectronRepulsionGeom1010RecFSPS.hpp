#ifndef ElectronRepulsionGeom1010RecFSPS_hpp
#define ElectronRepulsionGeom1010RecFSPS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSSXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FS|1/|r-r'||PS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fsps(T& distributor,
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

    CSimdArray<double> pbuffer(705, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(495, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(945, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(2835, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(189, 1);

    // setup Boys fuction data

    const CBoysFunc<6> bf_table;

    CSimdArray<double> bf_data(8, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 7, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 7, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 7, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 7);
                }

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 55, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 58, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 61, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 64, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 67, 1, 7, 10, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 76, 2, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 85, 3, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 94, 4, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 103, 10, 25, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 121, 13, 31, 37, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 139, 16, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 157, 19, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 175, 0, 1, 55, 58, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 181, 1, 2, 58, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 187, 2, 3, 61, 64, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 193, 7, 10, 58, 67, 76, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 211, 10, 13, 61, 76, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 229, 13, 16, 64, 85, 94, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 247, 25, 31, 76, 103, 121, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 283, 31, 37, 85, 121, 139, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 319, 37, 43, 94, 139, 157, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 355, 55, 58, 175, 181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 365, 58, 61, 181, 187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 375, 67, 76, 181, 193, 211, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 405, 76, 85, 187, 211, 229, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 435, 103, 121, 211, 247, 283, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 495, 121, 139, 229, 283, 319, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 555, 175, 181, 355, 365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 570, 193, 211, 365, 375, 405, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 615, 247, 283, 405, 435, 495, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1, pbuffer, 55, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4, pbuffer, 175, 6, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {55, 58});

                pbuffer.scale(2.0 * a_exp, {175, 181});

                pbuffer.scale(2.0 * a_exp, {355, 365});

                pbuffer.scale(2.0 * a_exp, {555, 570});

                t2cfunc::reduce(cbuffer, 110, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 111, pbuffer, 55, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 114, pbuffer, 175, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 120, pbuffer, 355, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 130, pbuffer, 555, 15, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {0, 1});

                pbuffer.scale(pfactors, 0, 2.0, {7, 10});

                pbuffer.scale(pfactors, 0, 2.0, {25, 31});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {55, 58});

                pbuffer.scale(pfactors, 0, 2.0, {67, 76});

                pbuffer.scale(pfactors, 0, 2.0, {103, 121});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {175, 181});

                pbuffer.scale(pfactors, 0, 2.0, {193, 211});

                pbuffer.scale(pfactors, 0, 2.0, {247, 283});

                t2cfunc::reduce(cbuffer, 10, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 11, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 14, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 20, pbuffer, 55, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 23, pbuffer, 67, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 32, pbuffer, 103, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 50, pbuffer, 175, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 56, pbuffer, 193, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 74, pbuffer, 247, 36, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {7, 10});

                pbuffer.scale(2.0 * a_exp, {25, 31});

                pbuffer.scale(2.0 * a_exp, {55, 58});

                pbuffer.scale(2.0 * a_exp, {67, 76});

                pbuffer.scale(2.0 * a_exp, {103, 121});

                pbuffer.scale(2.0 * a_exp, {175, 181});

                pbuffer.scale(2.0 * a_exp, {193, 211});

                pbuffer.scale(2.0 * a_exp, {247, 283});

                pbuffer.scale(pfactors, 0, 2.0, {355, 365});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {375, 405});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {435, 495});

                pbuffer.scale(pfactors, 0, 2.0, {555, 570});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {570, 615});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {615, 705});

                t2cfunc::reduce(cbuffer, 145, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 146, pbuffer, 7, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 149, pbuffer, 25, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 155, pbuffer, 55, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 158, pbuffer, 67, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 167, pbuffer, 103, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 185, pbuffer, 175, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 191, pbuffer, 193, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 209, pbuffer, 247, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 245, pbuffer, 355, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 255, pbuffer, 375, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 285, pbuffer, 435, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 345, pbuffer, 555, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 360, pbuffer, 570, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 405, pbuffer, 615, 90, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 10, 11, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3, cbuffer, 11, 14, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 12, cbuffer, 0, 0, 3, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 21, cbuffer, 20, 23, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 30, cbuffer, 23, 32, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 57, cbuffer, 1, 21, 30, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 84, cbuffer, 50, 56, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 102, cbuffer, 56, 74, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 156, cbuffer, 4, 84, 102, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 210, cbuffer, 145, 146, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 213, cbuffer, 146, 149, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 222, cbuffer, 110, 210, 213, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 231, cbuffer, 155, 158, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 240, cbuffer, 158, 167, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 267, cbuffer, 111, 231, 240, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 294, cbuffer, 185, 191, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 312, cbuffer, 191, 209, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 366, cbuffer, 114, 294, 312, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 420, cbuffer, 245, 255, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 450, cbuffer, 255, 285, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 540, cbuffer, 120, 420, 450, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 630, cbuffer, 345, 360, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 675, cbuffer, 360, 405, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 810, cbuffer, 130, 630, 675, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 0, ckbuffer, 12, 0, 0);

            t4cfunc::ket_transform<1, 0>(skbuffer, 3, ckbuffer, 15, 0, 0);

            t4cfunc::ket_transform<1, 0>(skbuffer, 6, ckbuffer, 18, 0, 0);

            t4cfunc::ket_transform<1, 0>(skbuffer, 36, ckbuffer, 57, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 45, ckbuffer, 66, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 54, ckbuffer, 75, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 144, ckbuffer, 156, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 162, ckbuffer, 174, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 180, ckbuffer, 192, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2520, ckbuffer, 222, 0, 0);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2523, ckbuffer, 225, 0, 0);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2526, ckbuffer, 228, 0, 0);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2529, ckbuffer, 267, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2538, ckbuffer, 276, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2547, ckbuffer, 285, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2556, ckbuffer, 366, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2574, ckbuffer, 384, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2592, ckbuffer, 402, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2610, ckbuffer, 540, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2640, ckbuffer, 570, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2670, ckbuffer, 600, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2700, ckbuffer, 810, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2745, ckbuffer, 855, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 2790, ckbuffer, 900, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 630, 0, 36, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 639, 3, 45, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 648, 6, 54, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 738, 36, 144, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 765, 45, 162, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 792, 54, 180, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 1548, 630, 738, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 1566, 639, 765, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 1584, 648, 792, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 9, 2520, 2529, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 63, 2529, 2556, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 198, 2556, 2610, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 360, 2610, 2700, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 657, 0, 9, 63, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 819, 36, 63, 198, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 1062, 144, 198, 360, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 1602, 630, 657, 819, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 1764, 738, 819, 1062, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fsxx(skbuffer, 2250, 1548, 1602, 1764, r_ab, 1, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 2250, 1, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 21, skbuffer, 2280, 1, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 42, skbuffer, 2310, 1, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 63, skbuffer, 2340, 1, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 84, skbuffer, 2370, 1, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 105, skbuffer, 2400, 1, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 126, skbuffer, 2430, 1, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 147, skbuffer, 2460, 1, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 168, skbuffer, 2490, 1, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 1, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFSPS_hpp */
