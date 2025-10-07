#ifndef ElectronRepulsionGeom1010RecDFPS_hpp
#define ElectronRepulsionGeom1010RecDFPS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
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
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DF|1/|r-r'||PS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dfps(T& distributor,
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

    CSimdArray<double> pbuffer(2105, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1089, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(2079, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(6048, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(945, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 75, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 78, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 81, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 84, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 87, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 90, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 93, 1, 9, 12, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 102, 2, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 111, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 120, 4, 18, 21, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 129, 5, 21, 24, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 138, 6, 24, 27, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 147, 12, 33, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 165, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 183, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 201, 21, 51, 57, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 219, 24, 57, 63, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 237, 27, 63, 69, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 255, 0, 1, 75, 78, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 261, 1, 2, 78, 81, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 267, 2, 3, 81, 84, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 273, 3, 4, 84, 87, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 279, 4, 5, 87, 90, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 285, 9, 12, 78, 93, 102, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 303, 12, 15, 81, 102, 111, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 321, 15, 18, 84, 111, 120, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 339, 18, 21, 87, 120, 129, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 357, 21, 24, 90, 129, 138, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 375, 33, 39, 102, 147, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 411, 39, 45, 111, 165, 183, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 447, 45, 51, 120, 183, 201, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 483, 51, 57, 129, 201, 219, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 519, 57, 63, 138, 219, 237, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 555, 75, 78, 255, 261, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 565, 78, 81, 261, 267, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 575, 81, 84, 267, 273, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 585, 84, 87, 273, 279, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 595, 93, 102, 261, 285, 303, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 625, 102, 111, 267, 303, 321, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 655, 111, 120, 273, 321, 339, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 685, 120, 129, 279, 339, 357, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 715, 147, 165, 303, 375, 411, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 775, 165, 183, 321, 411, 447, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 835, 183, 201, 339, 447, 483, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 895, 201, 219, 357, 483, 519, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 955, 255, 261, 555, 565, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 970, 261, 267, 565, 575, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 985, 267, 273, 575, 585, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1000, 285, 303, 565, 595, 625, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1045, 303, 321, 575, 625, 655, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1090, 321, 339, 585, 655, 685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1135, 375, 411, 625, 715, 775, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1225, 411, 447, 655, 775, 835, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1315, 447, 483, 685, 835, 895, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 1405, 555, 565, 955, 970, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 1426, 565, 575, 970, 985, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1447, 595, 625, 970, 1000, 1045, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1510, 625, 655, 985, 1045, 1090, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 1573, 715, 775, 1045, 1135, 1225, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 1699, 775, 835, 1090, 1225, 1315, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 1825, 955, 970, 1405, 1426, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 1853, 1000, 1045, 1426, 1447, 1510, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 1937, 1135, 1225, 1510, 1573, 1699, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 555, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 955, 15, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {555, 565});

                pbuffer.scale(2.0 * a_exp, {955, 970});

                pbuffer.scale(2.0 * a_exp, {1405, 1426});

                pbuffer.scale(2.0 * a_exp, {1825, 1853});

                t2cfunc::reduce(cbuffer, 275, pbuffer, 555, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 285, pbuffer, 955, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 300, pbuffer, 1405, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 321, pbuffer, 1825, 28, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {555, 565});

                pbuffer.scale(pfactors, 0, 2.0, {595, 625});

                pbuffer.scale(pfactors, 0, 2.0, {715, 775});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {955, 970});

                pbuffer.scale(pfactors, 0, 2.0, {1000, 1045});

                pbuffer.scale(pfactors, 0, 2.0, {1135, 1225});

                t2cfunc::reduce(cbuffer, 25, pbuffer, 555, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 35, pbuffer, 595, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 65, pbuffer, 715, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 125, pbuffer, 955, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 140, pbuffer, 1000, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 185, pbuffer, 1135, 90, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {555, 565});

                pbuffer.scale(2.0 * a_exp, {595, 625});

                pbuffer.scale(2.0 * a_exp, {715, 775});

                pbuffer.scale(2.0 * a_exp, {955, 970});

                pbuffer.scale(2.0 * a_exp, {1000, 1045});

                pbuffer.scale(2.0 * a_exp, {1135, 1225});

                pbuffer.scale(pfactors, 0, 2.0, {1405, 1426});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1447, 1510});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1573, 1699});

                pbuffer.scale(pfactors, 0, 2.0, {1825, 1853});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1853, 1937});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1937, 2105});

                t2cfunc::reduce(cbuffer, 349, pbuffer, 555, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 359, pbuffer, 595, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 389, pbuffer, 715, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 449, pbuffer, 955, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 464, pbuffer, 1000, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 509, pbuffer, 1135, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 599, pbuffer, 1405, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 620, pbuffer, 1447, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 683, pbuffer, 1573, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 809, pbuffer, 1825, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 837, pbuffer, 1853, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 921, pbuffer, 1937, 168, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 25, 35, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 30, cbuffer, 35, 65, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 120, cbuffer, 0, 0, 30, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 210, cbuffer, 125, 140, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 255, cbuffer, 140, 185, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 390, cbuffer, 10, 210, 255, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 525, cbuffer, 349, 359, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 555, cbuffer, 359, 389, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 645, cbuffer, 275, 525, 555, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 735, cbuffer, 449, 464, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 780, cbuffer, 464, 509, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 915, cbuffer, 285, 735, 780, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1050, cbuffer, 599, 620, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1113, cbuffer, 620, 683, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1302, cbuffer, 300, 1050, 1113, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1491, cbuffer, 809, 837, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1575, cbuffer, 837, 921, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1827, cbuffer, 321, 1491, 1575, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<1, 0>(skbuffer, 0, ckbuffer, 120, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 30, ckbuffer, 150, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 60, ckbuffer, 180, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 360, ckbuffer, 390, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 405, ckbuffer, 435, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 450, ckbuffer, 480, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5382, ckbuffer, 645, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5412, ckbuffer, 675, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5442, ckbuffer, 705, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5472, ckbuffer, 915, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5517, ckbuffer, 960, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5562, ckbuffer, 1005, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5607, ckbuffer, 1302, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5670, ckbuffer, 1365, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5733, ckbuffer, 1428, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5796, ckbuffer, 1827, 0, 6);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5880, ckbuffer, 1911, 0, 6);

            t4cfunc::ket_transform<1, 0>(skbuffer, 5964, ckbuffer, 1995, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 1467, 0, 360, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 1557, 30, 405, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 1647, 60, 450, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 90, 5382, 5472, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 495, 5472, 5607, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 900, 5607, 5796, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 1737, 0, 90, 495, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 2547, 360, 495, 900, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 3762, 1467, 1737, 2547, r_ab, 1, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 0, skbuffer, 3762, 1, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 105, skbuffer, 3942, 1, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 210, skbuffer, 4122, 1, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 315, skbuffer, 4302, 1, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 420, skbuffer, 4482, 1, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 525, skbuffer, 4662, 1, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 630, skbuffer, 4842, 1, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 735, skbuffer, 5022, 1, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 840, skbuffer, 5202, 1, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 3, 1, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDFPS_hpp */
