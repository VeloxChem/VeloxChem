#ifndef ElectronRepulsionGeom1010RecDSFD_hpp
#define ElectronRepulsionGeom1010RecDSFD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXDD.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXFD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSH.hpp"
#include "ElectronRepulsionGeom1010ContrRecDSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPSXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSSXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSI.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSI.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
#include "ElectronRepulsionPrimRecSPSI.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSH.hpp"
#include "ElectronRepulsionPrimRecSSSI.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DS|1/|r-r'||FD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dsfd(T& distributor,
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

    CSimdArray<double> pbuffer(2981, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2664, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(23688, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(11655, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1575, 1);

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

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 155, 37, 43, 85, 95, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 170, 43, 49, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 185, 49, 55, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 200, 55, 61, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 215, 61, 67, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 230, 67, 73, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 245, 85, 95, 155, 170, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 266, 95, 105, 170, 185, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 287, 105, 115, 185, 200, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 308, 115, 125, 200, 215, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 329, 125, 135, 215, 230, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 350, 155, 170, 245, 266, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 378, 170, 185, 266, 287, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 406, 185, 200, 287, 308, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssi(pbuffer, 434, 200, 215, 308, 329, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 462, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 465, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 474, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 483, 13, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 501, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 519, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 537, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 567, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 597, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 627, 95, 155, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 672, 105, 170, 185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 717, 115, 185, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 762, 170, 245, 266, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 825, 185, 266, 287, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 888, 200, 287, 308, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 951, 266, 350, 378, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1035, 287, 378, 406, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsi(pbuffer, 1119, 308, 406, 434, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1203, 13, 16, 462, 465, 474, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1221, 37, 43, 465, 483, 501, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1257, 43, 49, 474, 501, 519, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1293, 85, 95, 501, 537, 567, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1353, 95, 105, 519, 567, 597, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1413, 155, 170, 567, 627, 672, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1503, 170, 185, 597, 672, 717, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1593, 245, 266, 672, 762, 825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1719, 266, 287, 717, 825, 888, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 1845, 350, 378, 825, 951, 1035, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsi(pbuffer, 2013, 378, 406, 888, 1035, 1119, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2181, 483, 501, 1203, 1221, 1257, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2241, 537, 567, 1257, 1293, 1353, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2341, 627, 672, 1353, 1413, 1503, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 2491, 762, 825, 1503, 1593, 1719, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsi(pbuffer, 2701, 951, 1035, 1719, 1845, 2013, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 16, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 31, pbuffer, 483, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 49, pbuffer, 537, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 79, pbuffer, 627, 45, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {37, 43});

                pbuffer.scale(2.0 * a_exp, {85, 95});

                pbuffer.scale(2.0 * a_exp, {155, 170});

                pbuffer.scale(2.0 * a_exp, {483, 501});

                pbuffer.scale(2.0 * a_exp, {537, 567});

                pbuffer.scale(2.0 * a_exp, {627, 672});

                pbuffer.scale(2.0 * a_exp, {1221, 1257});

                pbuffer.scale(2.0 * a_exp, {1293, 1353});

                pbuffer.scale(2.0 * a_exp, {1413, 1503});

                pbuffer.scale(2.0 * a_exp, {2181, 2241});

                pbuffer.scale(2.0 * a_exp, {2241, 2341});

                pbuffer.scale(2.0 * a_exp, {2341, 2491});

                t2cfunc::reduce(cbuffer, 444, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 450, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 460, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 475, pbuffer, 483, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 493, pbuffer, 537, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 523, pbuffer, 627, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 568, pbuffer, 1221, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 604, pbuffer, 1293, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 664, pbuffer, 1413, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 754, pbuffer, 2181, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 814, pbuffer, 2241, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 914, pbuffer, 2341, 150, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {37, 43});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {85, 95});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {155, 170});

                pbuffer.scale(pfactors, 0, 2.0, {245, 266});

                pbuffer.scale(pfactors, 0, 2.0, {350, 378});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {483, 501});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {537, 567});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {627, 672});

                pbuffer.scale(pfactors, 0, 2.0, {762, 825});

                pbuffer.scale(pfactors, 0, 2.0, {951, 1035});

                t2cfunc::reduce(cbuffer, 124, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 130, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 140, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 155, pbuffer, 245, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 176, pbuffer, 350, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 204, pbuffer, 483, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 222, pbuffer, 537, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 252, pbuffer, 627, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 297, pbuffer, 762, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 360, pbuffer, 951, 84, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {37, 43});

                pbuffer.scale(2.0 * a_exp, {85, 95});

                pbuffer.scale(2.0 * a_exp, {155, 170});

                pbuffer.scale(2.0 * a_exp, {245, 266});

                pbuffer.scale(2.0 * a_exp, {350, 378});

                pbuffer.scale(2.0 * a_exp, {483, 501});

                pbuffer.scale(2.0 * a_exp, {537, 567});

                pbuffer.scale(2.0 * a_exp, {627, 672});

                pbuffer.scale(2.0 * a_exp, {762, 825});

                pbuffer.scale(2.0 * a_exp, {951, 1035});

                pbuffer.scale(pfactors, 0, 2.0, {1221, 1257});

                pbuffer.scale(pfactors, 0, 2.0, {1293, 1353});

                pbuffer.scale(pfactors, 0, 2.0, {1413, 1503});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1593, 1719});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {1845, 2013});

                pbuffer.scale(pfactors, 0, 2.0, {2181, 2241});

                pbuffer.scale(pfactors, 0, 2.0, {2241, 2341});

                pbuffer.scale(pfactors, 0, 2.0, {2341, 2491});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2491, 2701});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2701, 2981});

                t2cfunc::reduce(cbuffer, 1064, pbuffer, 37, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1070, pbuffer, 85, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1080, pbuffer, 155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1095, pbuffer, 245, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1116, pbuffer, 350, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1144, pbuffer, 483, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1162, pbuffer, 537, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1192, pbuffer, 627, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1237, pbuffer, 762, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1300, pbuffer, 951, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1384, pbuffer, 1221, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1420, pbuffer, 1293, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1480, pbuffer, 1413, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1570, pbuffer, 1593, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1696, pbuffer, 1845, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1864, pbuffer, 2181, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1924, pbuffer, 2241, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2024, pbuffer, 2341, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2174, pbuffer, 2491, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2384, pbuffer, 2701, 280, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 156, cbuffer, 0, 6, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 228, cbuffer, 6, 16, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 483, 156, 228, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1455, cbuffer, 31, 49, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1671, cbuffer, 49, 79, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 2436, 1455, 1671, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4104, cbuffer, 444, 450, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 4176, cbuffer, 450, 460, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 4431, 4104, 4176, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 5403, cbuffer, 475, 493, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 5619, cbuffer, 493, 523, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 6384, 5403, 5619, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8832, cbuffer, 568, 604, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 9264, cbuffer, 604, 664, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 10794, 8832, 9264, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 15378, cbuffer, 754, 814, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 16098, cbuffer, 814, 914, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 18648, 15378, 16098, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 0, cbuffer, 124, 130, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 18, cbuffer, 130, 140, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 48, cbuffer, 140, 155, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 93, cbuffer, 155, 176, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 174, cbuffer, 0, 0, 18, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 258, cbuffer, 6, 18, 48, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 348, cbuffer, 16, 48, 93, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 519, 156, 174, 258, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 627, 228, 258, 348, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 807, 483, 519, 627, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 987, cbuffer, 204, 222, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1041, cbuffer, 222, 252, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 1131, cbuffer, 252, 297, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 1266, cbuffer, 297, 360, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1509, cbuffer, 31, 987, 1041, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 1761, cbuffer, 49, 1041, 1131, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 2031, cbuffer, 79, 1131, 1266, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 2544, 1455, 1509, 1761, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 2868, 1671, 1761, 2031, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 3408, 2436, 2544, 2868, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3948, cbuffer, 1064, 1070, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3966, cbuffer, 1070, 1080, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 3996, cbuffer, 1080, 1095, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 4041, cbuffer, 1095, 1116, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4122, cbuffer, 444, 3948, 3966, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 4206, cbuffer, 450, 3966, 3996, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 4296, cbuffer, 460, 3996, 4041, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 4467, 4104, 4122, 4206, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 4575, 4176, 4206, 4296, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 4755, 4431, 4467, 4575, cfactors, 6, 0, 0);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4935, cbuffer, 1144, 1162, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 4989, cbuffer, 1162, 1192, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 5079, cbuffer, 1192, 1237, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 5214, cbuffer, 1237, 1300, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 5457, cbuffer, 475, 4935, 4989, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 5709, cbuffer, 493, 4989, 5079, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 5979, cbuffer, 523, 5079, 5214, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 6492, 5403, 5457, 5709, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 6816, 5619, 5709, 5979, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 7356, 6384, 6492, 6816, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 7896, cbuffer, 1384, 1420, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 8004, cbuffer, 1420, 1480, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 8184, cbuffer, 1480, 1570, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 8454, cbuffer, 1570, 1696, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 8940, cbuffer, 568, 7896, 8004, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 9444, cbuffer, 604, 8004, 8184, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 9984, cbuffer, 664, 8184, 8454, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 11010, 8832, 8940, 9444, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 11658, 9264, 9444, 9984, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 12738, 10794, 11010, 11658, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 13818, cbuffer, 1864, 1924, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 13998, cbuffer, 1924, 2024, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 14298, cbuffer, 2024, 2174, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsh(ckbuffer, 14748, cbuffer, 2174, 2384, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 15558, cbuffer, 754, 13818, 13998, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 16398, cbuffer, 814, 13998, 14298, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpg(ckbuffer, 17298, cbuffer, 914, 14298, 14748, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 19008, 15378, 15558, 16398, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdf(ckbuffer, 20088, 16098, 16398, 17298, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxfd(ckbuffer, 21888, 18648, 19008, 20088, cfactors, 6, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 0, ckbuffer, 807, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 35, ckbuffer, 867, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 70, ckbuffer, 927, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 420, ckbuffer, 3408, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 525, ckbuffer, 3588, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 630, ckbuffer, 3768, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 9555, ckbuffer, 4755, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 9590, ckbuffer, 4815, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 9625, ckbuffer, 4875, 0, 0);

            t4cfunc::ket_transform<3, 2>(skbuffer, 9660, ckbuffer, 7356, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 9765, ckbuffer, 7536, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 9870, ckbuffer, 7716, 0, 1);

            t4cfunc::ket_transform<3, 2>(skbuffer, 9975, ckbuffer, 12738, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 10185, ckbuffer, 13098, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 10395, ckbuffer, 13458, 0, 2);

            t4cfunc::ket_transform<3, 2>(skbuffer, 10605, ckbuffer, 21888, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 10955, ckbuffer, 22488, 0, 3);

            t4cfunc::ket_transform<3, 2>(skbuffer, 11305, ckbuffer, 23088, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 3570, 0, 420, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 3675, 35, 525, r_ab, 3, 2);

            erirec::comp_bra_hrr_electron_repulsion_psxx(skbuffer, 3780, 70, 630, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ssxx(skbuffer, 105, 9555, 9660, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 735, 9660, 9975, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 1680, 9975, 10605, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_psxx(skbuffer, 3885, 0, 105, 735, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 4830, 420, 735, 1680, r_ab, 3, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dsxx(skbuffer, 7665, 3570, 3885, 4830, r_ab, 3, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 0, skbuffer, 7665, 3, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 175, skbuffer, 7875, 3, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 350, skbuffer, 8085, 3, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 525, skbuffer, 8295, 3, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 700, skbuffer, 8505, 3, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 875, skbuffer, 8715, 3, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 1050, skbuffer, 8925, 3, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 1225, skbuffer, 9135, 3, 2);

            t4cfunc::bra_transform<2, 0>(sbuffer, 1400, skbuffer, 9345, 3, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 0, 3, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDSFD_hpp */
