#ifndef ElectronRepulsionGeom1010RecFDPP_hpp
#define ElectronRepulsionGeom1010RecFDPP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom1010ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FD|1/|r-r'||PP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fdpp(T& distributor,
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

    CSimdArray<double> pbuffer(4132, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2442, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(5994, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(29646, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2835, 1);

    // setup Boys fuction data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_data(11, ket_npgtos);

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
                    t4cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 10, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 10);
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

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 10, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 25, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 28, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 31, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 34, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 37, 0, 1, 10, 13, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 43, 1, 2, 13, 16, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 49, 2, 3, 16, 19, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 55, 3, 4, 19, 22, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 61, 4, 5, 22, 25, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 67, 5, 6, 25, 28, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 73, 6, 7, 28, 31, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 79, 7, 8, 31, 34, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 10, 13, 37, 43, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 13, 16, 43, 49, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 16, 19, 49, 55, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 19, 22, 55, 61, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 22, 25, 61, 67, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 25, 28, 67, 73, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 28, 31, 73, 79, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 155, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 158, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 161, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 164, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 167, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 170, 1, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 179, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 188, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 197, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 206, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 215, 6, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 224, 13, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 242, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 260, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 278, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 296, 25, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 314, 28, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 332, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 362, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 392, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 422, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 452, 67, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 482, 73, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 512, 1, 2, 155, 158, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 518, 2, 3, 158, 161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 524, 3, 4, 161, 164, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 530, 4, 5, 164, 167, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 536, 10, 13, 155, 170, 179, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 554, 13, 16, 158, 179, 188, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 572, 16, 19, 161, 188, 197, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 590, 19, 22, 164, 197, 206, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 608, 22, 25, 167, 206, 215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 626, 37, 43, 179, 224, 242, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 662, 43, 49, 188, 242, 260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 698, 49, 55, 197, 260, 278, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 734, 55, 61, 206, 278, 296, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 770, 61, 67, 215, 296, 314, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 806, 85, 95, 242, 332, 362, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 866, 95, 105, 260, 362, 392, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 926, 105, 115, 278, 392, 422, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 986, 115, 125, 296, 422, 452, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1046, 125, 135, 314, 452, 482, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1106, 155, 158, 512, 518, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1116, 158, 161, 518, 524, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1126, 161, 164, 524, 530, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1136, 170, 179, 512, 536, 554, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1166, 179, 188, 518, 554, 572, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1196, 188, 197, 524, 572, 590, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1226, 197, 206, 530, 590, 608, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1256, 224, 242, 554, 626, 662, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1316, 242, 260, 572, 662, 698, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1376, 260, 278, 590, 698, 734, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1436, 278, 296, 608, 734, 770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1496, 332, 362, 662, 806, 866, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1596, 362, 392, 698, 866, 926, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1696, 392, 422, 734, 926, 986, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1796, 422, 452, 770, 986, 1046, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1896, 512, 518, 1106, 1116, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1911, 518, 524, 1116, 1126, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1926, 536, 554, 1106, 1136, 1166, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1971, 554, 572, 1116, 1166, 1196, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2016, 572, 590, 1126, 1196, 1226, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2061, 626, 662, 1166, 1256, 1316, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2151, 662, 698, 1196, 1316, 1376, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2241, 698, 734, 1226, 1376, 1436, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2331, 806, 866, 1316, 1496, 1596, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2481, 866, 926, 1376, 1596, 1696, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2631, 926, 986, 1436, 1696, 1796, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 2781, 1106, 1116, 1896, 1911, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2802, 1136, 1166, 1896, 1926, 1971, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2865, 1166, 1196, 1911, 1971, 2016, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2928, 1256, 1316, 1971, 2061, 2151, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3054, 1316, 1376, 2016, 2151, 2241, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3180, 1496, 1596, 2151, 2331, 2481, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3390, 1596, 1696, 2241, 2481, 2631, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 3600, 1926, 1971, 2781, 2802, 2865, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 3684, 2061, 2151, 2865, 2928, 3054, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 3852, 2331, 2481, 3054, 3180, 3390, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 536, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 1136, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 48, pbuffer, 1926, 45, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {536, 554});

                pbuffer.scale(2.0 * a_exp, {1136, 1166});

                pbuffer.scale(2.0 * a_exp, {1926, 1971});

                pbuffer.scale(2.0 * a_exp, {2802, 2865});

                pbuffer.scale(2.0 * a_exp, {3600, 3684});

                t2cfunc::reduce(cbuffer, 682, pbuffer, 536, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 700, pbuffer, 1136, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 730, pbuffer, 1926, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 775, pbuffer, 2802, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 838, pbuffer, 3600, 84, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {536, 554});

                pbuffer.scale(pfactors, 0, 2.0, {626, 662});

                pbuffer.scale(pfactors, 0, 2.0, {806, 866});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1136, 1166});

                pbuffer.scale(pfactors, 0, 2.0, {1256, 1316});

                pbuffer.scale(pfactors, 0, 2.0, {1496, 1596});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1926, 1971});

                pbuffer.scale(pfactors, 0, 2.0, {2061, 2151});

                pbuffer.scale(pfactors, 0, 2.0, {2331, 2481});

                t2cfunc::reduce(cbuffer, 93, pbuffer, 536, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 111, pbuffer, 626, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 147, pbuffer, 806, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 207, pbuffer, 1136, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 237, pbuffer, 1256, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 297, pbuffer, 1496, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 397, pbuffer, 1926, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 442, pbuffer, 2061, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 532, pbuffer, 2331, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {536, 554});

                pbuffer.scale(2.0 * a_exp, {626, 662});

                pbuffer.scale(2.0 * a_exp, {806, 866});

                pbuffer.scale(2.0 * a_exp, {1136, 1166});

                pbuffer.scale(2.0 * a_exp, {1256, 1316});

                pbuffer.scale(2.0 * a_exp, {1496, 1596});

                pbuffer.scale(2.0 * a_exp, {1926, 1971});

                pbuffer.scale(2.0 * a_exp, {2061, 2151});

                pbuffer.scale(2.0 * a_exp, {2331, 2481});

                pbuffer.scale(pfactors, 0, 2.0, {2802, 2865});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2928, 3054});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3180, 3390});

                pbuffer.scale(pfactors, 0, 2.0, {3600, 3684});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3684, 3852});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3852, 4132});

                t2cfunc::reduce(cbuffer, 922, pbuffer, 536, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 940, pbuffer, 626, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 976, pbuffer, 806, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1036, pbuffer, 1136, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1066, pbuffer, 1256, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1126, pbuffer, 1496, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1226, pbuffer, 1926, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1271, pbuffer, 2061, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1361, pbuffer, 2331, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1511, pbuffer, 2802, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1574, pbuffer, 2928, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1700, pbuffer, 3180, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1910, pbuffer, 3600, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1994, pbuffer, 3684, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2162, pbuffer, 3852, 280, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 93, 111, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 54, cbuffer, 111, 147, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 162, cbuffer, 0, 0, 54, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 324, cbuffer, 207, 237, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 414, cbuffer, 237, 297, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 594, cbuffer, 18, 324, 414, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 864, cbuffer, 397, 442, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 999, cbuffer, 442, 532, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1269, cbuffer, 48, 864, 999, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1674, cbuffer, 922, 940, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1728, cbuffer, 940, 976, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1836, cbuffer, 682, 1674, 1728, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1998, cbuffer, 1036, 1066, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2088, cbuffer, 1066, 1126, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2268, cbuffer, 700, 1998, 2088, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2538, cbuffer, 1226, 1271, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2673, cbuffer, 1271, 1361, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2943, cbuffer, 730, 2538, 2673, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3348, cbuffer, 1511, 1574, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3537, cbuffer, 1574, 1700, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3915, cbuffer, 775, 3348, 3537, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4482, cbuffer, 1910, 1994, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4734, cbuffer, 1994, 2162, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 5238, cbuffer, 838, 4482, 4734, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 0, ckbuffer, 162, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 54, ckbuffer, 216, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 108, ckbuffer, 270, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 648, ckbuffer, 594, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 738, ckbuffer, 684, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 828, ckbuffer, 774, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 1728, ckbuffer, 1269, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 1863, ckbuffer, 1404, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 1998, ckbuffer, 1539, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 27486, ckbuffer, 1836, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 27540, ckbuffer, 1890, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 27594, ckbuffer, 1944, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 27648, ckbuffer, 2268, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 27738, ckbuffer, 2358, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 27828, ckbuffer, 2448, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 27918, ckbuffer, 2943, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 28053, ckbuffer, 3078, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 28188, ckbuffer, 3213, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 28323, ckbuffer, 3915, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 28512, ckbuffer, 4104, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 28701, ckbuffer, 4293, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 28890, ckbuffer, 5238, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 29142, ckbuffer, 5490, 0, 6);

            t4cfunc::ket_transform<1, 1>(skbuffer, 29394, ckbuffer, 5742, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 5049, 0, 648, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 5211, 54, 738, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 5373, 108, 828, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 6993, 648, 1728, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 7263, 738, 1863, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 7533, 828, 1998, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 13878, 5049, 6993, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 14202, 5211, 7263, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 14526, 5373, 7533, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 162, 27486, 27648, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 918, 27648, 27918, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 2133, 27918, 28323, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 3348, 28323, 28890, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 5535, 0, 162, 918, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 7803, 648, 918, 2133, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 10233, 1728, 2133, 3348, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 14850, 5049, 5535, 7803, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 17766, 6993, 7803, 10233, r_ab, 1, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fdxx(skbuffer, 22626, 13878, 14850, 17766, r_ab, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 22626, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 315, skbuffer, 23166, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 630, skbuffer, 23706, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 945, skbuffer, 24246, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1260, skbuffer, 24786, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1575, skbuffer, 25326, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1890, skbuffer, 25866, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2205, skbuffer, 26406, 1, 1);

            t4cfunc::bra_transform<3, 2>(sbuffer, 2520, skbuffer, 26946, 1, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 1, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFDPP_hpp */
