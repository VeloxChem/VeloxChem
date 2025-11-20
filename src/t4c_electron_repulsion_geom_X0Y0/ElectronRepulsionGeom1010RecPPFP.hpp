#ifndef ElectronRepulsionGeom1010RecPPFP_hpp
#define ElectronRepulsionGeom1010RecPPFP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXDP.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXFP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(PP|1/|r-r'||FP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_ppfp(T& distributor,
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

    CSimdArray<double> pbuffer(2011, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1628, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(12540, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(4788, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1701, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 294, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 297, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 300, 1, 9, 12, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 309, 2, 12, 15, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 318, 3, 15, 18, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 327, 12, 33, 39, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 345, 15, 39, 45, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 363, 18, 45, 51, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 381, 39, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 411, 45, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 441, 51, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 471, 85, 135, 150, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 516, 95, 150, 165, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 561, 105, 165, 180, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 606, 150, 210, 231, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 669, 165, 231, 252, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 732, 180, 252, 273, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 795, 1, 2, 294, 297, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 801, 9, 12, 294, 300, 309, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 819, 12, 15, 297, 309, 318, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 837, 33, 39, 309, 327, 345, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 873, 39, 45, 318, 345, 363, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 909, 75, 85, 345, 381, 411, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 969, 85, 95, 363, 411, 441, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1029, 135, 150, 411, 471, 516, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1119, 150, 165, 441, 516, 561, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1209, 210, 231, 516, 606, 669, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1335, 231, 252, 561, 669, 732, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1461, 300, 309, 795, 801, 819, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1491, 327, 345, 819, 837, 873, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1551, 381, 411, 873, 909, 969, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1651, 471, 516, 969, 1029, 1119, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 1801, 606, 669, 1119, 1209, 1335, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 300, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 327, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 27, pbuffer, 381, 30, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {300, 309});

                pbuffer.scale(2.0 * a_exp, {327, 345});

                pbuffer.scale(2.0 * a_exp, {381, 411});

                pbuffer.scale(2.0 * a_exp, {801, 819});

                pbuffer.scale(2.0 * a_exp, {837, 873});

                pbuffer.scale(2.0 * a_exp, {909, 969});

                pbuffer.scale(2.0 * a_exp, {1461, 1491});

                pbuffer.scale(2.0 * a_exp, {1491, 1551});

                pbuffer.scale(2.0 * a_exp, {1551, 1651});

                t2cfunc::reduce(cbuffer, 222, pbuffer, 300, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 231, pbuffer, 327, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 249, pbuffer, 381, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 279, pbuffer, 801, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 297, pbuffer, 837, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 333, pbuffer, 909, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 393, pbuffer, 1461, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 423, pbuffer, 1491, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 483, pbuffer, 1551, 100, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {300, 309});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {327, 345});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {381, 411});

                pbuffer.scale(pfactors, 0, 2.0, {471, 516});

                pbuffer.scale(pfactors, 0, 2.0, {606, 669});

                t2cfunc::reduce(cbuffer, 57, pbuffer, 300, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 66, pbuffer, 327, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 84, pbuffer, 381, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 114, pbuffer, 471, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 159, pbuffer, 606, 63, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {300, 309});

                pbuffer.scale(2.0 * a_exp, {327, 345});

                pbuffer.scale(2.0 * a_exp, {381, 411});

                pbuffer.scale(2.0 * a_exp, {471, 516});

                pbuffer.scale(2.0 * a_exp, {606, 669});

                pbuffer.scale(pfactors, 0, 2.0, {801, 819});

                pbuffer.scale(pfactors, 0, 2.0, {837, 873});

                pbuffer.scale(pfactors, 0, 2.0, {909, 969});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1029, 1119});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1209, 1335});

                pbuffer.scale(pfactors, 0, 2.0, {1461, 1491});

                pbuffer.scale(pfactors, 0, 2.0, {1491, 1551});

                pbuffer.scale(pfactors, 0, 2.0, {1551, 1651});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1651, 1801});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1801, 2011});

                t2cfunc::reduce(cbuffer, 583, pbuffer, 300, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 592, pbuffer, 327, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 610, pbuffer, 381, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 640, pbuffer, 471, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 685, pbuffer, 606, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 748, pbuffer, 801, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 766, pbuffer, 837, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 802, pbuffer, 909, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 862, pbuffer, 1029, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 952, pbuffer, 1209, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1078, pbuffer, 1461, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1108, pbuffer, 1491, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1168, pbuffer, 1551, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1268, pbuffer, 1651, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1418, pbuffer, 1801, 210, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 306, cbuffer, 0, 9, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 414, cbuffer, 9, 27, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 900, 306, 414, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 2016, cbuffer, 222, 231, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 2124, cbuffer, 231, 249, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 2610, 2016, 2124, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 4032, cbuffer, 279, 297, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4248, cbuffer, 297, 333, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 5220, 4032, 4248, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 7860, cbuffer, 393, 423, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8220, cbuffer, 423, 483, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdp(ckbuffer, 9840, 7860, 8220, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 0, cbuffer, 57, 66, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 27, cbuffer, 66, 84, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 81, cbuffer, 84, 114, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 171, cbuffer, 114, 159, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 333, cbuffer, 0, 0, 27, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 468, cbuffer, 9, 27, 81, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 630, cbuffer, 27, 81, 171, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 954, 306, 333, 468, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 1116, 414, 468, 630, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 1440, 900, 954, 1116, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1710, cbuffer, 583, 592, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1737, cbuffer, 592, 610, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1791, cbuffer, 610, 640, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 1881, cbuffer, 640, 685, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2043, cbuffer, 222, 1710, 1737, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 2178, cbuffer, 231, 1737, 1791, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2340, cbuffer, 249, 1791, 1881, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 2664, 2016, 2043, 2178, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 2826, 2124, 2178, 2340, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 3150, 2610, 2664, 2826, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3420, cbuffer, 748, 766, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3474, cbuffer, 766, 802, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3582, cbuffer, 802, 862, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 3762, cbuffer, 862, 952, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4086, cbuffer, 279, 3420, 3474, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4356, cbuffer, 297, 3474, 3582, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 4680, cbuffer, 333, 3582, 3762, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 5328, 4032, 4086, 4356, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 5652, 4248, 4356, 4680, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 6300, 5220, 5328, 5652, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 6840, cbuffer, 1078, 1108, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6930, cbuffer, 1108, 1168, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 7110, cbuffer, 1168, 1268, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 7410, cbuffer, 1268, 1418, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 7950, cbuffer, 393, 6840, 6930, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 8400, cbuffer, 423, 6930, 7110, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 8940, cbuffer, 483, 7110, 7410, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdp(ckbuffer, 10020, 7860, 7950, 8400, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 10560, 8220, 8400, 8940, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfp(ckbuffer, 11640, 9840, 10020, 10560, cfactors, 6, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 0, ckbuffer, 1440, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 63, ckbuffer, 1530, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 126, ckbuffer, 1620, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 3591, ckbuffer, 3150, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 3654, ckbuffer, 3240, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 3717, ckbuffer, 3330, 0, 1);

            t4cfunc::ket_transform<3, 1>(skbuffer, 3780, ckbuffer, 6300, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 3906, ckbuffer, 6480, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 4032, ckbuffer, 6660, 0, 2);

            t4cfunc::ket_transform<3, 1>(skbuffer, 4158, ckbuffer, 11640, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 4368, ckbuffer, 11940, 0, 3);

            t4cfunc::ket_transform<3, 1>(skbuffer, 4578, ckbuffer, 12240, 0, 3);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 189, 3591, 3780, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 756, 3780, 4158, r_ab, 3, 1);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 1890, 0, 189, 756, r_ab, 3, 1);

            t4cfunc::bra_transform<1, 1>(sbuffer, 0, skbuffer, 1890, 3, 1);

            t4cfunc::bra_transform<1, 1>(sbuffer, 189, skbuffer, 2079, 3, 1);

            t4cfunc::bra_transform<1, 1>(sbuffer, 378, skbuffer, 2268, 3, 1);

            t4cfunc::bra_transform<1, 1>(sbuffer, 567, skbuffer, 2457, 3, 1);

            t4cfunc::bra_transform<1, 1>(sbuffer, 756, skbuffer, 2646, 3, 1);

            t4cfunc::bra_transform<1, 1>(sbuffer, 945, skbuffer, 2835, 3, 1);

            t4cfunc::bra_transform<1, 1>(sbuffer, 1134, skbuffer, 3024, 3, 1);

            t4cfunc::bra_transform<1, 1>(sbuffer, 1323, skbuffer, 3213, 3, 1);

            t4cfunc::bra_transform<1, 1>(sbuffer, 1512, skbuffer, 3402, 3, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 1, 1, 3, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecPPFP_hpp */
