#ifndef ElectronRepulsionGeom1010RecDPDS_hpp
#define ElectronRepulsionGeom1010RecDPDS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDS.hpp"
#include "ElectronRepulsionGeom1010ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DP|1/|r-r'||DS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dpds(T& distributor,
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

    CSimdArray<double> pbuffer(1415, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1032, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(3741, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(3660, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(675, 1);

    // setup Boys fuction data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 8);

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

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 65, 8, 11, 29, 35, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 75, 11, 14, 35, 41, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 14, 17, 41, 47, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 17, 20, 47, 53, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 20, 23, 53, 59, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 115, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 118, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 121, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 124, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 127, 1, 8, 11, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 136, 2, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 145, 3, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 154, 4, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 163, 11, 29, 35, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 181, 14, 35, 41, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 199, 17, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 217, 20, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 235, 35, 65, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 265, 41, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 295, 47, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 325, 53, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 355, 0, 1, 115, 118, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 361, 1, 2, 118, 121, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 367, 2, 3, 121, 124, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 373, 8, 11, 118, 127, 136, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 391, 11, 14, 121, 136, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 409, 14, 17, 124, 145, 154, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 427, 29, 35, 136, 163, 181, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 463, 35, 41, 145, 181, 199, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 499, 41, 47, 154, 199, 217, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 535, 65, 75, 181, 235, 265, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 595, 75, 85, 199, 265, 295, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 655, 85, 95, 217, 295, 325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 715, 115, 118, 355, 361, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 725, 118, 121, 361, 367, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 735, 127, 136, 361, 373, 391, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 765, 136, 145, 367, 391, 409, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 795, 163, 181, 391, 427, 463, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 855, 181, 199, 409, 463, 499, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 915, 235, 265, 463, 535, 595, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1015, 265, 295, 499, 595, 655, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1115, 355, 361, 715, 725, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1130, 373, 391, 725, 735, 765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1175, 427, 463, 765, 795, 855, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 1265, 535, 595, 855, 915, 1015, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 115, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 127, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 12, pbuffer, 355, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 373, 18, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {115, 118});

                pbuffer.scale(2.0 * a_exp, {127, 136});

                pbuffer.scale(2.0 * a_exp, {355, 361});

                pbuffer.scale(2.0 * a_exp, {373, 391});

                pbuffer.scale(2.0 * a_exp, {715, 725});

                pbuffer.scale(2.0 * a_exp, {735, 765});

                pbuffer.scale(2.0 * a_exp, {1115, 1130});

                pbuffer.scale(2.0 * a_exp, {1130, 1175});

                t2cfunc::reduce(cbuffer, 216, pbuffer, 115, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 219, pbuffer, 127, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 228, pbuffer, 355, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 234, pbuffer, 373, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 252, pbuffer, 715, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 262, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 292, pbuffer, 1115, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 307, pbuffer, 1130, 45, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {115, 118});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {127, 136});

                pbuffer.scale(pfactors, 0, 2.0, {163, 181});

                pbuffer.scale(pfactors, 0, 2.0, {235, 265});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {355, 361});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {373, 391});

                pbuffer.scale(pfactors, 0, 2.0, {427, 463});

                pbuffer.scale(pfactors, 0, 2.0, {535, 595});

                t2cfunc::reduce(cbuffer, 36, pbuffer, 115, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 39, pbuffer, 127, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 48, pbuffer, 163, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 66, pbuffer, 235, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 96, pbuffer, 355, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 102, pbuffer, 373, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 120, pbuffer, 427, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 156, pbuffer, 535, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {115, 118});

                pbuffer.scale(2.0 * a_exp, {127, 136});

                pbuffer.scale(2.0 * a_exp, {163, 181});

                pbuffer.scale(2.0 * a_exp, {235, 265});

                pbuffer.scale(2.0 * a_exp, {355, 361});

                pbuffer.scale(2.0 * a_exp, {373, 391});

                pbuffer.scale(2.0 * a_exp, {427, 463});

                pbuffer.scale(2.0 * a_exp, {535, 595});

                pbuffer.scale(pfactors, 0, 2.0, {715, 725});

                pbuffer.scale(pfactors, 0, 2.0, {735, 765});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {795, 855});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {915, 1015});

                pbuffer.scale(pfactors, 0, 2.0, {1115, 1130});

                pbuffer.scale(pfactors, 0, 2.0, {1130, 1175});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1175, 1265});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1265, 1415});

                t2cfunc::reduce(cbuffer, 352, pbuffer, 115, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 355, pbuffer, 127, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 364, pbuffer, 163, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 382, pbuffer, 235, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 412, pbuffer, 355, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 418, pbuffer, 373, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 436, pbuffer, 427, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 472, pbuffer, 535, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 532, pbuffer, 715, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 542, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 572, pbuffer, 795, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 632, pbuffer, 915, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 732, pbuffer, 1115, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 747, pbuffer, 1130, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 792, pbuffer, 1175, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 882, pbuffer, 1265, 150, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 90, cbuffer, 0, 3, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 441, cbuffer, 12, 18, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 873, cbuffer, 216, 219, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1224, cbuffer, 228, 234, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1866, cbuffer, 252, 262, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2886, cbuffer, 292, 307, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 36, 39, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 9, cbuffer, 39, 48, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 36, cbuffer, 48, 66, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 99, cbuffer, 0, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 126, cbuffer, 3, 9, 36, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 207, 90, 99, 126, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 261, cbuffer, 96, 102, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 279, cbuffer, 102, 120, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 333, cbuffer, 120, 156, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 459, cbuffer, 12, 261, 279, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 513, cbuffer, 18, 279, 333, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 675, 441, 459, 513, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 783, cbuffer, 352, 355, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 792, cbuffer, 355, 364, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 819, cbuffer, 364, 382, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 882, cbuffer, 216, 783, 792, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 909, cbuffer, 219, 792, 819, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 990, 873, 882, 909, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1044, cbuffer, 412, 418, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1062, cbuffer, 418, 436, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1116, cbuffer, 436, 472, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1242, cbuffer, 228, 1044, 1062, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1296, cbuffer, 234, 1062, 1116, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1458, 1224, 1242, 1296, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1566, cbuffer, 532, 542, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1596, cbuffer, 542, 572, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1686, cbuffer, 572, 632, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1896, cbuffer, 252, 1566, 1596, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1986, cbuffer, 262, 1596, 1686, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 2256, 1866, 1896, 1986, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2436, cbuffer, 732, 747, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2481, cbuffer, 747, 792, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2616, cbuffer, 792, 882, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2931, cbuffer, 292, 2436, 2481, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3066, cbuffer, 307, 2481, 2616, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 3471, 2886, 2931, 3066, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 0, ckbuffer, 207, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15, ckbuffer, 225, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 30, ckbuffer, 243, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 180, ckbuffer, 675, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 210, ckbuffer, 711, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 240, ckbuffer, 747, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 990, ckbuffer, 0, 1, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1035, ckbuffer, 54, 1, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1080, ckbuffer, 108, 1, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3150, ckbuffer, 990, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3165, ckbuffer, 1008, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3180, ckbuffer, 1026, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3195, ckbuffer, 1458, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3225, ckbuffer, 1494, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3255, ckbuffer, 1530, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3285, ckbuffer, 2256, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3335, ckbuffer, 2316, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3385, ckbuffer, 2376, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3435, ckbuffer, 3471, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3510, ckbuffer, 3561, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 3585, ckbuffer, 3651, 0, 4);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 45, 3150, 3195, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 270, 3195, 3285, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 540, 3285, 3435, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 1125, 0, 45, 270, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 1530, 180, 270, 540, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 2340, 990, 1125, 1530, r_ab, 2, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 2340, 2, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 75, skbuffer, 2430, 2, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 150, skbuffer, 2520, 2, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 225, skbuffer, 2610, 2, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 300, skbuffer, 2700, 2, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 375, skbuffer, 2790, 2, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 450, skbuffer, 2880, 2, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 525, skbuffer, 2970, 2, 0);

            t4cfunc::bra_transform<2, 1>(sbuffer, 600, skbuffer, 3060, 2, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 2, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDPDS_hpp */
