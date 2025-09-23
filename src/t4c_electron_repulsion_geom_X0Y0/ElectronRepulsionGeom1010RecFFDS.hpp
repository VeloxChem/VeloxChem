#ifndef ElectronRepulsionGeom1010RecFFDS_hpp
#define ElectronRepulsionGeom1010RecFFDS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
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
#include "ElectronRepulsionPrimRecSKSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FF|1/|r-r'||DS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_ffds(T& distributor,
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

    CSimdArray<double> pbuffer(6615, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(3744, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(13572, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(25155, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(2205, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 175, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 178, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 181, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 184, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 187, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 190, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 193, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 196, 1, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 205, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 214, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 223, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 232, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 241, 6, 26, 29, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 250, 7, 29, 32, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 259, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 277, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 295, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 313, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 331, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 349, 29, 71, 77, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 367, 32, 77, 83, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 385, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 415, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 445, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 475, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 505, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 535, 77, 145, 155, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 565, 83, 155, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 595, 0, 1, 175, 178, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 601, 1, 2, 178, 181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 607, 2, 3, 181, 184, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 613, 3, 4, 184, 187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 619, 4, 5, 187, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 625, 5, 6, 190, 193, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 631, 11, 14, 178, 196, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 649, 14, 17, 181, 205, 214, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 667, 17, 20, 184, 214, 223, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 685, 20, 23, 187, 223, 232, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 703, 23, 26, 190, 232, 241, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 721, 26, 29, 193, 241, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 739, 41, 47, 205, 259, 277, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 775, 47, 53, 214, 277, 295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 811, 53, 59, 223, 295, 313, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 847, 59, 65, 232, 313, 331, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 883, 65, 71, 241, 331, 349, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 919, 71, 77, 250, 349, 367, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 955, 95, 105, 277, 385, 415, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1015, 105, 115, 295, 415, 445, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1075, 115, 125, 313, 445, 475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1135, 125, 135, 331, 475, 505, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1195, 135, 145, 349, 505, 535, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1255, 145, 155, 367, 535, 565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1315, 175, 178, 595, 601, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1325, 178, 181, 601, 607, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1335, 181, 184, 607, 613, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1345, 184, 187, 613, 619, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1355, 187, 190, 619, 625, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1365, 196, 205, 601, 631, 649, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1395, 205, 214, 607, 649, 667, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1425, 214, 223, 613, 667, 685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1455, 223, 232, 619, 685, 703, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1485, 232, 241, 625, 703, 721, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1515, 259, 277, 649, 739, 775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1575, 277, 295, 667, 775, 811, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1635, 295, 313, 685, 811, 847, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1695, 313, 331, 703, 847, 883, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1755, 331, 349, 721, 883, 919, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1815, 385, 415, 775, 955, 1015, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1915, 415, 445, 811, 1015, 1075, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2015, 445, 475, 847, 1075, 1135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2115, 475, 505, 883, 1135, 1195, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2215, 505, 535, 919, 1195, 1255, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2315, 595, 601, 1315, 1325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2330, 601, 607, 1325, 1335, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2345, 607, 613, 1335, 1345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 2360, 613, 619, 1345, 1355, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2375, 631, 649, 1325, 1365, 1395, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2420, 649, 667, 1335, 1395, 1425, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2465, 667, 685, 1345, 1425, 1455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2510, 685, 703, 1355, 1455, 1485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2555, 739, 775, 1395, 1515, 1575, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2645, 775, 811, 1425, 1575, 1635, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2735, 811, 847, 1455, 1635, 1695, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2825, 847, 883, 1485, 1695, 1755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2915, 955, 1015, 1575, 1815, 1915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3065, 1015, 1075, 1635, 1915, 2015, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3215, 1075, 1135, 1695, 2015, 2115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3365, 1135, 1195, 1755, 2115, 2215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 3515, 1315, 1325, 2315, 2330, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 3536, 1325, 1335, 2330, 2345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 3557, 1335, 1345, 2345, 2360, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3578, 1365, 1395, 2330, 2375, 2420, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3641, 1395, 1425, 2345, 2420, 2465, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 3704, 1425, 1455, 2360, 2465, 2510, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3767, 1515, 1575, 2420, 2555, 2645, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3893, 1575, 1635, 2465, 2645, 2735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 4019, 1635, 1695, 2510, 2735, 2825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4145, 1815, 1915, 2645, 2915, 3065, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4355, 1915, 2015, 2735, 3065, 3215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 4565, 2015, 2115, 2825, 3215, 3365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 4775, 2315, 2330, 3515, 3536, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 4803, 2330, 2345, 3536, 3557, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 4831, 2375, 2420, 3536, 3578, 3641, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 4915, 2420, 2465, 3557, 3641, 3704, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 4999, 2555, 2645, 3641, 3767, 3893, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 5167, 2645, 2735, 3704, 3893, 4019, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 5335, 2915, 3065, 3893, 4145, 4355, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 5615, 3065, 3215, 4019, 4355, 4565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_skss(pbuffer, 5895, 3515, 3536, 4775, 4803, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 5931, 3578, 3641, 4803, 4831, 4915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 6039, 3767, 3893, 4915, 4999, 5167, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksf(pbuffer, 6255, 4145, 4355, 5167, 5335, 5615, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1315, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 1365, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 40, pbuffer, 2315, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 55, pbuffer, 2375, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 100, pbuffer, 3515, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 121, pbuffer, 3578, 63, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1315, 1325});

                pbuffer.scale(2.0 * a_exp, {1365, 1395});

                pbuffer.scale(2.0 * a_exp, {2315, 2330});

                pbuffer.scale(2.0 * a_exp, {2375, 2420});

                pbuffer.scale(2.0 * a_exp, {3515, 3536});

                pbuffer.scale(2.0 * a_exp, {3578, 3641});

                pbuffer.scale(2.0 * a_exp, {4775, 4803});

                pbuffer.scale(2.0 * a_exp, {4831, 4915});

                pbuffer.scale(2.0 * a_exp, {5895, 5931});

                pbuffer.scale(2.0 * a_exp, {5931, 6039});

                t2cfunc::reduce(cbuffer, 1104, pbuffer, 1315, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1114, pbuffer, 1365, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1144, pbuffer, 2315, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1159, pbuffer, 2375, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1204, pbuffer, 3515, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1225, pbuffer, 3578, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1288, pbuffer, 4775, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1316, pbuffer, 4831, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1400, pbuffer, 5895, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1436, pbuffer, 5931, 108, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1315, 1325});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1365, 1395});

                pbuffer.scale(pfactors, 0, 2.0, {1515, 1575});

                pbuffer.scale(pfactors, 0, 2.0, {1815, 1915});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2315, 2330});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2375, 2420});

                pbuffer.scale(pfactors, 0, 2.0, {2555, 2645});

                pbuffer.scale(pfactors, 0, 2.0, {2915, 3065});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3515, 3536});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {3578, 3641});

                pbuffer.scale(pfactors, 0, 2.0, {3767, 3893});

                pbuffer.scale(pfactors, 0, 2.0, {4145, 4355});

                t2cfunc::reduce(cbuffer, 184, pbuffer, 1315, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 194, pbuffer, 1365, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 224, pbuffer, 1515, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 284, pbuffer, 1815, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 384, pbuffer, 2315, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 399, pbuffer, 2375, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 444, pbuffer, 2555, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 534, pbuffer, 2915, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 684, pbuffer, 3515, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 705, pbuffer, 3578, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 768, pbuffer, 3767, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 894, pbuffer, 4145, 210, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1315, 1325});

                pbuffer.scale(2.0 * a_exp, {1365, 1395});

                pbuffer.scale(2.0 * a_exp, {1515, 1575});

                pbuffer.scale(2.0 * a_exp, {1815, 1915});

                pbuffer.scale(2.0 * a_exp, {2315, 2330});

                pbuffer.scale(2.0 * a_exp, {2375, 2420});

                pbuffer.scale(2.0 * a_exp, {2555, 2645});

                pbuffer.scale(2.0 * a_exp, {2915, 3065});

                pbuffer.scale(2.0 * a_exp, {3515, 3536});

                pbuffer.scale(2.0 * a_exp, {3578, 3641});

                pbuffer.scale(2.0 * a_exp, {3767, 3893});

                pbuffer.scale(2.0 * a_exp, {4145, 4355});

                pbuffer.scale(pfactors, 0, 2.0, {4775, 4803});

                pbuffer.scale(pfactors, 0, 2.0, {4831, 4915});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {4999, 5167});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5335, 5615});

                pbuffer.scale(pfactors, 0, 2.0, {5895, 5931});

                pbuffer.scale(pfactors, 0, 2.0, {5931, 6039});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {6039, 6255});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {6255, 6615});

                t2cfunc::reduce(cbuffer, 1544, pbuffer, 1315, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1554, pbuffer, 1365, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1584, pbuffer, 1515, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1644, pbuffer, 1815, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1744, pbuffer, 2315, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1759, pbuffer, 2375, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1804, pbuffer, 2555, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1894, pbuffer, 2915, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2044, pbuffer, 3515, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2065, pbuffer, 3578, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2128, pbuffer, 3767, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2254, pbuffer, 4145, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2464, pbuffer, 4775, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2492, pbuffer, 4831, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2576, pbuffer, 4999, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2744, pbuffer, 5335, 280, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3024, pbuffer, 5895, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3060, pbuffer, 5931, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3168, pbuffer, 6039, 216, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3384, pbuffer, 6255, 360, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 300, cbuffer, 0, 10, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1320, cbuffer, 40, 55, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2805, cbuffer, 100, 121, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 4302, cbuffer, 1104, 1114, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 5322, cbuffer, 1144, 1159, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 6807, cbuffer, 1204, 1225, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 8844, cbuffer, 1288, 1316, cfactors, 6, 0, 6);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 11520, cbuffer, 1400, 1436, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 184, 194, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 30, cbuffer, 194, 224, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 120, cbuffer, 224, 284, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 330, cbuffer, 0, 0, 30, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 420, cbuffer, 10, 30, 120, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 690, 300, 330, 420, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 870, cbuffer, 384, 399, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 915, cbuffer, 399, 444, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1050, cbuffer, 444, 534, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1365, cbuffer, 40, 870, 915, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1500, cbuffer, 55, 915, 1050, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1905, 1320, 1365, 1500, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2175, cbuffer, 684, 705, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2238, cbuffer, 705, 768, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2427, cbuffer, 768, 894, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2868, cbuffer, 100, 2175, 2238, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3057, cbuffer, 121, 2238, 2427, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 3624, 2805, 2868, 3057, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 4002, cbuffer, 1544, 1554, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4032, cbuffer, 1554, 1584, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4122, cbuffer, 1584, 1644, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 4332, cbuffer, 1104, 4002, 4032, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4422, cbuffer, 1114, 4032, 4122, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 4692, 4302, 4332, 4422, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 4872, cbuffer, 1744, 1759, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4917, cbuffer, 1759, 1804, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5052, cbuffer, 1804, 1894, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 5367, cbuffer, 1144, 4872, 4917, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 5502, cbuffer, 1159, 4917, 5052, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 5907, 5322, 5367, 5502, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 6177, cbuffer, 2044, 2065, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 6240, cbuffer, 2065, 2128, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6429, cbuffer, 2128, 2254, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 6870, cbuffer, 1204, 6177, 6240, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 7059, cbuffer, 1225, 6240, 6429, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 7626, 6807, 6870, 7059, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 8004, cbuffer, 2464, 2492, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 8088, cbuffer, 2492, 2576, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 8340, cbuffer, 2576, 2744, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 8928, cbuffer, 1288, 8004, 8088, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 9180, cbuffer, 1316, 8088, 8340, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 9936, 8844, 8928, 9180, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 10440, cbuffer, 3024, 3060, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 10548, cbuffer, 3060, 3168, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 10872, cbuffer, 3168, 3384, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 11628, cbuffer, 1400, 10440, 10548, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 11952, cbuffer, 1436, 10548, 10872, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 12924, 11520, 11628, 11952, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<2, 0>(skbuffer, 0, ckbuffer, 690, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 50, ckbuffer, 750, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 100, ckbuffer, 810, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 600, ckbuffer, 1905, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 675, ckbuffer, 1995, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 750, ckbuffer, 2085, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1500, ckbuffer, 3624, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1605, ckbuffer, 3750, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1710, ckbuffer, 3876, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 23505, ckbuffer, 4692, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 23555, ckbuffer, 4752, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 23605, ckbuffer, 4812, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 23655, ckbuffer, 5907, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 23730, ckbuffer, 5997, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 23805, ckbuffer, 6087, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 23880, ckbuffer, 7626, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 23985, ckbuffer, 7752, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 24090, ckbuffer, 7878, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 24195, ckbuffer, 9936, 0, 6);

            t4cfunc::ket_transform<2, 0>(skbuffer, 24335, ckbuffer, 10104, 0, 6);

            t4cfunc::ket_transform<2, 0>(skbuffer, 24475, ckbuffer, 10272, 0, 6);

            t4cfunc::ket_transform<2, 0>(skbuffer, 24615, ckbuffer, 12924, 0, 7);

            t4cfunc::ket_transform<2, 0>(skbuffer, 24795, ckbuffer, 13140, 0, 7);

            t4cfunc::ket_transform<2, 0>(skbuffer, 24975, ckbuffer, 13356, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 4020, 0, 600, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 4170, 50, 675, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 4320, 100, 750, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 5820, 600, 1500, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 6045, 675, 1605, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 6270, 750, 1710, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 11355, 4020, 5820, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 11655, 4170, 6045, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 11955, 4320, 6270, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 150, 23505, 23655, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 825, 23655, 23880, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 1815, 23880, 24195, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sixx(skbuffer, 2760, 24195, 24615, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 4470, 0, 150, 825, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 6495, 600, 825, 1815, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_phxx(skbuffer, 8520, 1500, 1815, 2760, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 12255, 4020, 4470, 6495, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dgxx(skbuffer, 14955, 5820, 6495, 8520, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ffxx(skbuffer, 19005, 11355, 12255, 14955, r_ab, 2, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 19005, 2, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 245, skbuffer, 19505, 2, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 490, skbuffer, 20005, 2, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 735, skbuffer, 20505, 2, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 980, skbuffer, 21005, 2, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1225, skbuffer, 21505, 2, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1470, skbuffer, 22005, 2, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1715, skbuffer, 22505, 2, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1960, skbuffer, 23005, 2, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 2, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFFDS_hpp */
