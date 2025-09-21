#ifndef ElectronRepulsionGeom1010RecFPPS_hpp
#define ElectronRepulsionGeom1010RecFPPS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
#include "ElectronRepulsionGeom1010ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FP|1/|r-r'||PS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fpps(T& distributor,
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

    CSimdArray<double> pbuffer(1265, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(814, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(1554, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(5796, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(567, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 65, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 68, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 71, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 74, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 77, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 80, 1, 8, 11, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 89, 2, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 98, 3, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 107, 4, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 116, 5, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 125, 11, 29, 35, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 143, 14, 35, 41, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 161, 17, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 179, 20, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 197, 23, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 215, 0, 1, 65, 68, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 221, 1, 2, 68, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 227, 2, 3, 71, 74, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 233, 3, 4, 74, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 239, 8, 11, 68, 80, 89, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 257, 11, 14, 71, 89, 98, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 275, 14, 17, 74, 98, 107, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 293, 17, 20, 77, 107, 116, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 311, 29, 35, 89, 125, 143, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 347, 35, 41, 98, 143, 161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 383, 41, 47, 107, 161, 179, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 419, 47, 53, 116, 179, 197, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 455, 65, 68, 215, 221, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 465, 68, 71, 221, 227, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 475, 71, 74, 227, 233, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 485, 80, 89, 221, 239, 257, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 515, 89, 98, 227, 257, 275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 545, 98, 107, 233, 275, 293, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 575, 125, 143, 257, 311, 347, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 635, 143, 161, 275, 347, 383, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 695, 161, 179, 293, 383, 419, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 755, 215, 221, 455, 465, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 770, 221, 227, 465, 475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 785, 239, 257, 465, 485, 515, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 830, 257, 275, 475, 515, 545, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 875, 311, 347, 515, 575, 635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 965, 347, 383, 545, 635, 695, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 1055, 455, 465, 755, 770, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1076, 485, 515, 770, 785, 830, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 1139, 575, 635, 830, 875, 965, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 65, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 215, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 455, 10, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {65, 68});

                pbuffer.scale(2.0 * a_exp, {215, 221});

                pbuffer.scale(2.0 * a_exp, {455, 465});

                pbuffer.scale(2.0 * a_exp, {755, 770});

                pbuffer.scale(2.0 * a_exp, {1055, 1076});

                t2cfunc::reduce(cbuffer, 209, pbuffer, 65, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 212, pbuffer, 215, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 218, pbuffer, 455, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 228, pbuffer, 755, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 243, pbuffer, 1055, 21, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {65, 68});

                pbuffer.scale(pfactors, 0, 2.0, {80, 89});

                pbuffer.scale(pfactors, 0, 2.0, {125, 143});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {215, 221});

                pbuffer.scale(pfactors, 0, 2.0, {239, 257});

                pbuffer.scale(pfactors, 0, 2.0, {311, 347});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {455, 465});

                pbuffer.scale(pfactors, 0, 2.0, {485, 515});

                pbuffer.scale(pfactors, 0, 2.0, {575, 635});

                t2cfunc::reduce(cbuffer, 19, pbuffer, 65, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 22, pbuffer, 80, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 31, pbuffer, 125, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 49, pbuffer, 215, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 55, pbuffer, 239, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 73, pbuffer, 311, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 109, pbuffer, 455, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 119, pbuffer, 485, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 149, pbuffer, 575, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {65, 68});

                pbuffer.scale(2.0 * a_exp, {80, 89});

                pbuffer.scale(2.0 * a_exp, {125, 143});

                pbuffer.scale(2.0 * a_exp, {215, 221});

                pbuffer.scale(2.0 * a_exp, {239, 257});

                pbuffer.scale(2.0 * a_exp, {311, 347});

                pbuffer.scale(2.0 * a_exp, {455, 465});

                pbuffer.scale(2.0 * a_exp, {485, 515});

                pbuffer.scale(2.0 * a_exp, {575, 635});

                pbuffer.scale(pfactors, 0, 2.0, {755, 770});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {785, 830});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {875, 965});

                pbuffer.scale(pfactors, 0, 2.0, {1055, 1076});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1076, 1139});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1139, 1265});

                t2cfunc::reduce(cbuffer, 264, pbuffer, 65, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 267, pbuffer, 80, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 276, pbuffer, 125, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 294, pbuffer, 215, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 300, pbuffer, 239, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 318, pbuffer, 311, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 354, pbuffer, 455, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 364, pbuffer, 485, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 394, pbuffer, 575, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 454, pbuffer, 755, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 469, pbuffer, 785, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 514, pbuffer, 875, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 604, pbuffer, 1055, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 625, pbuffer, 1076, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 688, pbuffer, 1139, 126, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 19, 22, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9, cbuffer, 22, 31, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 36, cbuffer, 0, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 63, cbuffer, 49, 55, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 81, cbuffer, 55, 73, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 135, cbuffer, 3, 63, 81, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 189, cbuffer, 109, 119, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 219, cbuffer, 119, 149, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 309, cbuffer, 9, 189, 219, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 399, cbuffer, 264, 267, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 408, cbuffer, 267, 276, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 435, cbuffer, 209, 399, 408, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 462, cbuffer, 294, 300, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 480, cbuffer, 300, 318, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 534, cbuffer, 212, 462, 480, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 588, cbuffer, 354, 364, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 618, cbuffer, 364, 394, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 708, cbuffer, 218, 588, 618, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 798, cbuffer, 454, 469, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 843, cbuffer, 469, 514, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 978, cbuffer, 228, 798, 843, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1113, cbuffer, 604, 625, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1176, cbuffer, 625, 688, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1365, cbuffer, 243, 1113, 1176, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 0, ckbuffer, 36, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 9, ckbuffer, 45, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 18, ckbuffer, 54, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 108, ckbuffer, 135, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 126, ckbuffer, 153, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 144, ckbuffer, 171, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 324, ckbuffer, 309, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 354, ckbuffer, 339, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 384, ckbuffer, 369, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5301, ckbuffer, 435, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5310, ckbuffer, 444, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5319, ckbuffer, 453, 0, 1);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5328, ckbuffer, 534, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5346, ckbuffer, 552, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5364, ckbuffer, 570, 0, 2);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5382, ckbuffer, 708, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5412, ckbuffer, 738, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5442, ckbuffer, 768, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5472, ckbuffer, 978, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5517, ckbuffer, 1023, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5562, ckbuffer, 1068, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5607, ckbuffer, 1365, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5670, ckbuffer, 1428, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5733, ckbuffer, 1491, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1089, 0, 108, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1116, 9, 126, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1143, 18, 144, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 1413, 108, 324, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 1467, 126, 354, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 1521, 144, 384, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 2871, 1089, 1413, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 2925, 1116, 1467, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 2979, 1143, 1521, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 27, 5301, 5328, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 162, 5328, 5382, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 414, 5382, 5472, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 684, 5472, 5607, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 1170, 0, 27, 162, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 1575, 108, 162, 414, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 2061, 324, 414, 684, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 3033, 1089, 1170, 1575, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 3519, 1413, 1575, 2061, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fpxx(skbuffer, 4491, 2871, 3033, 3519, r_ab, 1, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 4491, 1, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 63, skbuffer, 4581, 1, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 126, skbuffer, 4671, 1, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 189, skbuffer, 4761, 1, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 252, skbuffer, 4851, 1, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 315, skbuffer, 4941, 1, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 378, skbuffer, 5031, 1, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 441, skbuffer, 5121, 1, 0);

            t4cfunc::bra_transform<3, 1>(sbuffer, 504, skbuffer, 5211, 1, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 1, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFPPS_hpp */
