#ifndef ElectronRepulsionGeom1010RecFFPP_hpp
#define ElectronRepulsionGeom1010RecFFPP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPHXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSIXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSF.hpp"
#include "ElectronRepulsionPrimRecSKSP.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FF|1/|r-r'||PP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_ffpp(T& distributor,
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

    CSimdArray<double> pbuffer(6496, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(3432, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(8424, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(45279, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(3969, 1);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 11, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 11);
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

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 9, pfactors, 16, bf_data, 9);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 10, pfactors, 16, bf_data, 10);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 29, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 32, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 35, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 38, 9, 10, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 41, 0, 1, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 47, 1, 2, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 2, 3, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 3, 4, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 65, 4, 5, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 71, 5, 6, 26, 29, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 77, 6, 7, 29, 32, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 83, 7, 8, 32, 35, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 89, 8, 9, 35, 38, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 11, 14, 41, 47, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 14, 17, 47, 53, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 17, 20, 53, 59, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 20, 23, 59, 65, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 23, 26, 65, 71, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 26, 29, 71, 77, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 155, 29, 32, 77, 83, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 165, 32, 35, 83, 89, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 175, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 178, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 181, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 184, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 187, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 190, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 193, 1, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 202, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 211, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 220, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 229, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 238, 6, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 247, 7, 29, 32, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 256, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 274, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 292, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 310, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 328, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 346, 29, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 364, 32, 77, 83, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 382, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 412, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 442, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 472, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 502, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 532, 77, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 562, 83, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 592, 1, 2, 175, 178, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 598, 2, 3, 178, 181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 604, 3, 4, 181, 184, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 610, 4, 5, 184, 187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 616, 5, 6, 187, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 622, 11, 14, 175, 193, 202, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 640, 14, 17, 178, 202, 211, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 658, 17, 20, 181, 211, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 676, 20, 23, 184, 220, 229, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 694, 23, 26, 187, 229, 238, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 712, 26, 29, 190, 238, 247, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 730, 41, 47, 202, 256, 274, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 766, 47, 53, 211, 274, 292, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 802, 53, 59, 220, 292, 310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 838, 59, 65, 229, 310, 328, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 874, 65, 71, 238, 328, 346, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 910, 71, 77, 247, 346, 364, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 946, 95, 105, 274, 382, 412, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1006, 105, 115, 292, 412, 442, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1066, 115, 125, 310, 442, 472, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1126, 125, 135, 328, 472, 502, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1186, 135, 145, 346, 502, 532, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1246, 145, 155, 364, 532, 562, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1306, 175, 178, 592, 598, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1316, 178, 181, 598, 604, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1326, 181, 184, 604, 610, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1336, 184, 187, 610, 616, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1346, 193, 202, 592, 622, 640, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1376, 202, 211, 598, 640, 658, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1406, 211, 220, 604, 658, 676, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1436, 220, 229, 610, 676, 694, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1466, 229, 238, 616, 694, 712, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1496, 256, 274, 640, 730, 766, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1556, 274, 292, 658, 766, 802, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1616, 292, 310, 676, 802, 838, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1676, 310, 328, 694, 838, 874, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1736, 328, 346, 712, 874, 910, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1796, 382, 412, 766, 946, 1006, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1896, 412, 442, 802, 1006, 1066, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1996, 442, 472, 838, 1066, 1126, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2096, 472, 502, 874, 1126, 1186, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2196, 502, 532, 910, 1186, 1246, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2296, 592, 598, 1306, 1316, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2311, 598, 604, 1316, 1326, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2326, 604, 610, 1326, 1336, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2341, 622, 640, 1306, 1346, 1376, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2386, 640, 658, 1316, 1376, 1406, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2431, 658, 676, 1326, 1406, 1436, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2476, 676, 694, 1336, 1436, 1466, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2521, 730, 766, 1376, 1496, 1556, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2611, 766, 802, 1406, 1556, 1616, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2701, 802, 838, 1436, 1616, 1676, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2791, 838, 874, 1466, 1676, 1736, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2881, 946, 1006, 1556, 1796, 1896, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3031, 1006, 1066, 1616, 1896, 1996, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3181, 1066, 1126, 1676, 1996, 2096, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3331, 1126, 1186, 1736, 2096, 2196, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 3481, 1306, 1316, 2296, 2311, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 3502, 1316, 1326, 2311, 2326, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3523, 1346, 1376, 2296, 2341, 2386, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3586, 1376, 1406, 2311, 2386, 2431, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3649, 1406, 1436, 2326, 2431, 2476, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3712, 1496, 1556, 2386, 2521, 2611, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3838, 1556, 1616, 2431, 2611, 2701, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3964, 1616, 1676, 2476, 2701, 2791, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4090, 1796, 1896, 2611, 2881, 3031, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4300, 1896, 1996, 2701, 3031, 3181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4510, 1996, 2096, 2791, 3181, 3331, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 4720, 2296, 2311, 3481, 3502, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 4748, 2341, 2386, 3481, 3523, 3586, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 4832, 2386, 2431, 3502, 3586, 3649, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 4916, 2521, 2611, 3586, 3712, 3838, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 5084, 2611, 2701, 3649, 3838, 3964, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 5252, 2881, 3031, 3838, 4090, 4300, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 5532, 3031, 3181, 3964, 4300, 4510, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 5812, 3523, 3586, 4720, 4748, 4832, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 5920, 3712, 3838, 4832, 4916, 5084, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 6136, 4090, 4300, 5084, 5252, 5532, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1346, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 30, pbuffer, 2341, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 75, pbuffer, 3523, 63, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1346, 1376});

                pbuffer.scale(2.0 * a_exp, {2341, 2386});

                pbuffer.scale(2.0 * a_exp, {3523, 3586});

                pbuffer.scale(2.0 * a_exp, {4748, 4832});

                pbuffer.scale(2.0 * a_exp, {5812, 5920});

                t2cfunc::reduce(cbuffer, 1012, pbuffer, 1346, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1042, pbuffer, 2341, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1087, pbuffer, 3523, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1150, pbuffer, 4748, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1234, pbuffer, 5812, 108, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1346, 1376});

                pbuffer.scale(pfactors, 0, 2.0, {1496, 1556});

                pbuffer.scale(pfactors, 0, 2.0, {1796, 1896});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2341, 2386});

                pbuffer.scale(pfactors, 0, 2.0, {2521, 2611});

                pbuffer.scale(pfactors, 0, 2.0, {2881, 3031});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3523, 3586});

                pbuffer.scale(pfactors, 0, 2.0, {3712, 3838});

                pbuffer.scale(pfactors, 0, 2.0, {4090, 4300});

                t2cfunc::reduce(cbuffer, 138, pbuffer, 1346, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 168, pbuffer, 1496, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 228, pbuffer, 1796, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 328, pbuffer, 2341, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 373, pbuffer, 2521, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 463, pbuffer, 2881, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 613, pbuffer, 3523, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 676, pbuffer, 3712, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 802, pbuffer, 4090, 210, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1346, 1376});

                pbuffer.scale(2.0 * a_exp, {1496, 1556});

                pbuffer.scale(2.0 * a_exp, {1796, 1896});

                pbuffer.scale(2.0 * a_exp, {2341, 2386});

                pbuffer.scale(2.0 * a_exp, {2521, 2611});

                pbuffer.scale(2.0 * a_exp, {2881, 3031});

                pbuffer.scale(2.0 * a_exp, {3523, 3586});

                pbuffer.scale(2.0 * a_exp, {3712, 3838});

                pbuffer.scale(2.0 * a_exp, {4090, 4300});

                pbuffer.scale(pfactors, 0, 2.0, {4748, 4832});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {4916, 5084});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5252, 5532});

                pbuffer.scale(pfactors, 0, 2.0, {5812, 5920});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5920, 6136});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {6136, 6496});

                t2cfunc::reduce(cbuffer, 1342, pbuffer, 1346, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1372, pbuffer, 1496, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1432, pbuffer, 1796, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1532, pbuffer, 2341, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1577, pbuffer, 2521, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1667, pbuffer, 2881, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1817, pbuffer, 3523, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1880, pbuffer, 3712, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2006, pbuffer, 4090, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2216, pbuffer, 4748, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2300, pbuffer, 4916, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2468, pbuffer, 5252, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2748, pbuffer, 5812, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2856, pbuffer, 5920, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3072, pbuffer, 6136, 360, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 138, 168, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 90, cbuffer, 168, 228, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 270, cbuffer, 0, 0, 90, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 540, cbuffer, 328, 373, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 675, cbuffer, 373, 463, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 945, cbuffer, 30, 540, 675, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1350, cbuffer, 613, 676, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1539, cbuffer, 676, 802, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1917, cbuffer, 75, 1350, 1539, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2484, cbuffer, 1342, 1372, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2574, cbuffer, 1372, 1432, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2754, cbuffer, 1012, 2484, 2574, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3024, cbuffer, 1532, 1577, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3159, cbuffer, 1577, 1667, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3429, cbuffer, 1042, 3024, 3159, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3834, cbuffer, 1817, 1880, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4023, cbuffer, 1880, 2006, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4401, cbuffer, 1087, 3834, 4023, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4968, cbuffer, 2216, 2300, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5220, cbuffer, 2300, 2468, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 5724, cbuffer, 1150, 4968, 5220, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 6480, cbuffer, 2748, 2856, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6804, cbuffer, 2856, 3072, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 7452, cbuffer, 1234, 6480, 6804, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<1, 1>(skbuffer, 0, ckbuffer, 270, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 90, ckbuffer, 360, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 180, ckbuffer, 450, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 1080, ckbuffer, 945, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 1215, ckbuffer, 1080, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 1350, ckbuffer, 1215, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2700, ckbuffer, 1917, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 2889, ckbuffer, 2106, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 3078, ckbuffer, 2295, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 42309, ckbuffer, 2754, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 42399, ckbuffer, 2844, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 42489, ckbuffer, 2934, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 42579, ckbuffer, 3429, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 42714, ckbuffer, 3564, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 42849, ckbuffer, 3699, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 42984, ckbuffer, 4401, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 43173, ckbuffer, 4590, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 43362, ckbuffer, 4779, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 43551, ckbuffer, 5724, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 43803, ckbuffer, 5976, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 44055, ckbuffer, 6228, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 44307, ckbuffer, 7452, 0, 7);

            t4cfunc::ket_transform<1, 1>(skbuffer, 44631, ckbuffer, 7776, 0, 7);

            t4cfunc::ket_transform<1, 1>(skbuffer, 44955, ckbuffer, 8100, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 7236, 0, 1080, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 7506, 90, 1215, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 7776, 180, 1350, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 10476, 1080, 2700, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 10881, 1215, 2889, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 11286, 1350, 3078, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 20439, 7236, 10476, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 20979, 7506, 10881, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 21519, 7776, 11286, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 270, 42309, 42579, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 1485, 42579, 42984, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 3267, 42984, 43551, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sixx(skbuffer, 4968, 43551, 44307, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 8046, 0, 270, 1485, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 11691, 1080, 1485, 3267, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_phxx(skbuffer, 15336, 2700, 3267, 4968, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 22059, 7236, 8046, 11691, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dgxx(skbuffer, 26919, 10476, 11691, 15336, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ffxx(skbuffer, 34209, 20439, 22059, 26919, r_ab, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 34209, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 441, skbuffer, 35109, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 882, skbuffer, 36009, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1323, skbuffer, 36909, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1764, skbuffer, 37809, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2205, skbuffer, 38709, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 2646, skbuffer, 39609, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 3087, skbuffer, 40509, 1, 1);

            t4cfunc::bra_transform<3, 3>(sbuffer, 3528, skbuffer, 41409, 1, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 1, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFFPP_hpp */
