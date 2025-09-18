#ifndef ElectronRepulsionGeom1010RecFSDS_hpp
#define ElectronRepulsionGeom1010RecFSDS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FS|1/|r-r'||DS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fsds(T& distributor,
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

    CSimdArray<double> cbuffer(1080, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(3915, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(4725, 1);

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

                t2cfunc::reduce(cbuffer, 0, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1, pbuffer, 8, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4, pbuffer, 115, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 7, pbuffer, 127, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 355, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 22, pbuffer, 373, 18, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {8, 11});

                pbuffer.scale(2.0 * a_exp, {115, 118});

                pbuffer.scale(2.0 * a_exp, {127, 136});

                pbuffer.scale(2.0 * a_exp, {355, 361});

                pbuffer.scale(2.0 * a_exp, {373, 391});

                pbuffer.scale(2.0 * a_exp, {715, 725});

                pbuffer.scale(2.0 * a_exp, {735, 765});

                pbuffer.scale(2.0 * a_exp, {1115, 1130});

                pbuffer.scale(2.0 * a_exp, {1130, 1175});

                t2cfunc::reduce(cbuffer, 240, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 241, pbuffer, 8, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 244, pbuffer, 115, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 247, pbuffer, 127, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 256, pbuffer, 355, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 262, pbuffer, 373, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 280, pbuffer, 715, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 290, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 320, pbuffer, 1115, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 335, pbuffer, 1130, 45, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {0, 1});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {8, 11});

                pbuffer.scale(pfactors, 0, 2.0, {29, 35});

                pbuffer.scale(pfactors, 0, 2.0, {65, 75});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {115, 118});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {127, 136});

                pbuffer.scale(pfactors, 0, 2.0, {163, 181});

                pbuffer.scale(pfactors, 0, 2.0, {235, 265});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {355, 361});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {373, 391});

                pbuffer.scale(pfactors, 0, 2.0, {427, 463});

                pbuffer.scale(pfactors, 0, 2.0, {535, 595});

                t2cfunc::reduce(cbuffer, 40, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 41, pbuffer, 8, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 44, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 50, pbuffer, 65, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 60, pbuffer, 115, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 63, pbuffer, 127, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 72, pbuffer, 163, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 235, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 120, pbuffer, 355, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 126, pbuffer, 373, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 144, pbuffer, 427, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 180, pbuffer, 535, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {0, 1});

                pbuffer.scale(2.0 * a_exp, {8, 11});

                pbuffer.scale(2.0 * a_exp, {29, 35});

                pbuffer.scale(2.0 * a_exp, {65, 75});

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

                t2cfunc::reduce(cbuffer, 380, pbuffer, 0, 1, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 381, pbuffer, 8, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 384, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 390, pbuffer, 65, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 400, pbuffer, 115, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 403, pbuffer, 127, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 412, pbuffer, 163, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 430, pbuffer, 235, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 460, pbuffer, 355, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 466, pbuffer, 373, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 484, pbuffer, 427, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 520, pbuffer, 535, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 580, pbuffer, 715, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 590, pbuffer, 735, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 620, pbuffer, 795, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 680, pbuffer, 915, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 780, pbuffer, 1115, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 795, pbuffer, 1130, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 840, pbuffer, 1175, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 930, pbuffer, 1265, 150, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 30, cbuffer, 0, 1, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 177, cbuffer, 4, 7, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 528, cbuffer, 16, 22, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 900, cbuffer, 240, 241, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1047, cbuffer, 244, 247, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1398, cbuffer, 256, 262, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2040, cbuffer, 280, 290, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 3060, cbuffer, 320, 335, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 40, 41, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3, cbuffer, 41, 44, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 12, cbuffer, 44, 50, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 33, cbuffer, 0, 0, 3, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 42, cbuffer, 1, 3, 12, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 69, 30, 33, 42, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 87, cbuffer, 60, 63, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 96, cbuffer, 63, 72, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 123, cbuffer, 72, 90, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 186, cbuffer, 4, 87, 96, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 213, cbuffer, 7, 96, 123, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 294, 177, 186, 213, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 348, cbuffer, 120, 126, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 366, cbuffer, 126, 144, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 420, cbuffer, 144, 180, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 546, cbuffer, 16, 348, 366, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 600, cbuffer, 22, 366, 420, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 762, 528, 546, 600, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 870, cbuffer, 380, 381, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 873, cbuffer, 381, 384, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 882, cbuffer, 384, 390, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 903, cbuffer, 240, 870, 873, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 912, cbuffer, 241, 873, 882, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 939, 900, 903, 912, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 957, cbuffer, 400, 403, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 966, cbuffer, 403, 412, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 993, cbuffer, 412, 430, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1056, cbuffer, 244, 957, 966, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1083, cbuffer, 247, 966, 993, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1164, 1047, 1056, 1083, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1218, cbuffer, 460, 466, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1236, cbuffer, 466, 484, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1290, cbuffer, 484, 520, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1416, cbuffer, 256, 1218, 1236, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1470, cbuffer, 262, 1236, 1290, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1632, 1398, 1416, 1470, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1740, cbuffer, 580, 590, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1770, cbuffer, 590, 620, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1860, cbuffer, 620, 680, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2070, cbuffer, 280, 1740, 1770, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2160, cbuffer, 290, 1770, 1860, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 2430, 2040, 2070, 2160, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2610, cbuffer, 780, 795, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2655, cbuffer, 795, 840, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2790, cbuffer, 840, 930, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 3105, cbuffer, 320, 2610, 2655, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3240, cbuffer, 335, 2655, 2790, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 3645, 3060, 3105, 3240, cfactors, 6, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 0, ckbuffer, 69, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 5, ckbuffer, 75, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 10, ckbuffer, 81, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 60, ckbuffer, 294, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 75, ckbuffer, 312, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 90, ckbuffer, 330, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 240, ckbuffer, 762, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 270, ckbuffer, 798, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 300, ckbuffer, 834, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4200, ckbuffer, 939, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4205, ckbuffer, 945, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4210, ckbuffer, 951, 0, 0);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4215, ckbuffer, 1164, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4230, ckbuffer, 1182, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4245, ckbuffer, 1200, 0, 1);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4260, ckbuffer, 1632, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4290, ckbuffer, 1668, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4320, ckbuffer, 1704, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4350, ckbuffer, 2430, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4400, ckbuffer, 2490, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4450, ckbuffer, 2550, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4500, ckbuffer, 3645, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4575, ckbuffer, 3735, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 4650, ckbuffer, 3825, 0, 4);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 1050, 0, 60, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 1065, 5, 75, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 1080, 10, 90, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1230, 60, 240, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1275, 75, 270, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1320, 90, 300, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 2580, 1050, 1230, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 2610, 1065, 1275, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_dsxx(skbuffer, 2640, 1080, 1320, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 15, 4200, 4215, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 105, 4215, 4260, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 330, 4260, 4350, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 600, 4350, 4500, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 1095, 0, 15, 105, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 1365, 60, 105, 330, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 1770, 240, 330, 600, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 2670, 1050, 1095, 1365, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 2940, 1230, 1365, 1770, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fsxx(skbuffer, 3750, 2580, 2670, 2940, r_ab, 2, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 0, skbuffer, 3750, 2, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 35, skbuffer, 3800, 2, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 70, skbuffer, 3850, 2, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 105, skbuffer, 3900, 2, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 140, skbuffer, 3950, 2, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 175, skbuffer, 4000, 2, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 210, skbuffer, 4050, 2, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 245, skbuffer, 4100, 2, 0);

            t4cfunc::bra_transform<3, 0>(sbuffer, 280, skbuffer, 4150, 2, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 0, 2, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFSDS_hpp */
