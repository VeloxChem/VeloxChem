#ifndef ElectronRepulsionGeom1010RecFPDD_hpp
#define ElectronRepulsionGeom1010RecFPDD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSF.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSG.hpp"
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
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
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

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FP|1/|r-r'||DD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fpdd(T& distributor,
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

    CSimdArray<double> pbuffer(6872, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(5032, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(26862, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(48300, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(4725, 1);

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

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 175, 41, 47, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 190, 47, 53, 105, 115, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 205, 53, 59, 115, 125, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 220, 59, 65, 125, 135, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 235, 65, 71, 135, 145, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 250, 71, 77, 145, 155, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 265, 77, 83, 155, 165, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 280, 95, 105, 175, 190, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 301, 105, 115, 190, 205, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 322, 115, 125, 205, 220, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 343, 125, 135, 220, 235, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 364, 135, 145, 235, 250, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssh(pbuffer, 385, 145, 155, 250, 265, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 406, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 409, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 412, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 415, 2, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 424, 3, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 433, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 442, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 451, 14, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 469, 17, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 487, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 505, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 523, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 541, 47, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 571, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 601, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 631, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 661, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 691, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 736, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 781, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 826, 135, 220, 235, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 871, 145, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 916, 190, 280, 301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 979, 205, 301, 322, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1042, 220, 322, 343, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1105, 235, 343, 364, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1168, 250, 364, 385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1231, 2, 3, 406, 409, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 1237, 3, 4, 409, 412, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1243, 14, 17, 406, 415, 424, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1261, 17, 20, 409, 424, 433, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1279, 20, 23, 412, 433, 442, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1297, 41, 47, 415, 451, 469, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1333, 47, 53, 424, 469, 487, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1369, 53, 59, 433, 487, 505, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1405, 59, 65, 442, 505, 523, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1441, 95, 105, 469, 541, 571, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1501, 105, 115, 487, 571, 601, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1561, 115, 125, 505, 601, 631, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1621, 125, 135, 523, 631, 661, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1681, 175, 190, 571, 691, 736, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1771, 190, 205, 601, 736, 781, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1861, 205, 220, 631, 781, 826, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1951, 220, 235, 661, 826, 871, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2041, 280, 301, 736, 916, 979, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2167, 301, 322, 781, 979, 1042, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2293, 322, 343, 826, 1042, 1105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2419, 343, 364, 871, 1105, 1168, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 2545, 406, 409, 1231, 1237, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2555, 415, 424, 1231, 1243, 1261, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 2585, 424, 433, 1237, 1261, 1279, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2615, 451, 469, 1243, 1297, 1333, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2675, 469, 487, 1261, 1333, 1369, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2735, 487, 505, 1279, 1369, 1405, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2795, 541, 571, 1333, 1441, 1501, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2895, 571, 601, 1369, 1501, 1561, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2995, 601, 631, 1405, 1561, 1621, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3095, 691, 736, 1501, 1681, 1771, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3245, 736, 781, 1561, 1771, 1861, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 3395, 781, 826, 1621, 1861, 1951, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3545, 916, 979, 1771, 2041, 2167, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3755, 979, 1042, 1861, 2167, 2293, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3965, 1042, 1105, 1951, 2293, 2419, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 4175, 1243, 1261, 2545, 2555, 2585, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4220, 1297, 1333, 2555, 2615, 2675, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 4310, 1333, 1369, 2585, 2675, 2735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4400, 1441, 1501, 2675, 2795, 2895, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 4550, 1501, 1561, 2735, 2895, 2995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4700, 1681, 1771, 2895, 3095, 3245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 4925, 1771, 1861, 2995, 3245, 3395, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 5150, 2041, 2167, 3245, 3545, 3755, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 5465, 2167, 2293, 3395, 3755, 3965, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 5780, 2615, 2675, 4175, 4220, 4310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 5906, 2795, 2895, 4310, 4400, 4550, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsg(pbuffer, 6116, 3095, 3245, 4550, 4700, 4925, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsh(pbuffer, 6431, 3545, 3755, 4925, 5150, 5465, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 451, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 541, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 48, pbuffer, 1297, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 84, pbuffer, 1441, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 144, pbuffer, 2615, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 204, pbuffer, 2795, 100, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {451, 469});

                pbuffer.scale(2.0 * a_exp, {541, 571});

                pbuffer.scale(2.0 * a_exp, {1297, 1333});

                pbuffer.scale(2.0 * a_exp, {1441, 1501});

                pbuffer.scale(2.0 * a_exp, {2615, 2675});

                pbuffer.scale(2.0 * a_exp, {2795, 2895});

                pbuffer.scale(2.0 * a_exp, {4220, 4310});

                pbuffer.scale(2.0 * a_exp, {4400, 4550});

                pbuffer.scale(2.0 * a_exp, {5780, 5906});

                pbuffer.scale(2.0 * a_exp, {5906, 6116});

                t2cfunc::reduce(cbuffer, 1292, pbuffer, 451, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1310, pbuffer, 541, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1340, pbuffer, 1297, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1376, pbuffer, 1441, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1436, pbuffer, 2615, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1496, pbuffer, 2795, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1596, pbuffer, 4220, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1686, pbuffer, 4400, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1836, pbuffer, 5780, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1962, pbuffer, 5906, 210, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {451, 469});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {541, 571});

                pbuffer.scale(pfactors, 0, 2.0, {691, 736});

                pbuffer.scale(pfactors, 0, 2.0, {916, 979});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1297, 1333});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1441, 1501});

                pbuffer.scale(pfactors, 0, 2.0, {1681, 1771});

                pbuffer.scale(pfactors, 0, 2.0, {2041, 2167});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2615, 2675});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {2795, 2895});

                pbuffer.scale(pfactors, 0, 2.0, {3095, 3245});

                pbuffer.scale(pfactors, 0, 2.0, {3545, 3755});

                t2cfunc::reduce(cbuffer, 304, pbuffer, 451, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 322, pbuffer, 541, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 352, pbuffer, 691, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 397, pbuffer, 916, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 460, pbuffer, 1297, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 496, pbuffer, 1441, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 556, pbuffer, 1681, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 646, pbuffer, 2041, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 772, pbuffer, 2615, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 832, pbuffer, 2795, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 932, pbuffer, 3095, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1082, pbuffer, 3545, 210, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {451, 469});

                pbuffer.scale(2.0 * a_exp, {541, 571});

                pbuffer.scale(2.0 * a_exp, {691, 736});

                pbuffer.scale(2.0 * a_exp, {916, 979});

                pbuffer.scale(2.0 * a_exp, {1297, 1333});

                pbuffer.scale(2.0 * a_exp, {1441, 1501});

                pbuffer.scale(2.0 * a_exp, {1681, 1771});

                pbuffer.scale(2.0 * a_exp, {2041, 2167});

                pbuffer.scale(2.0 * a_exp, {2615, 2675});

                pbuffer.scale(2.0 * a_exp, {2795, 2895});

                pbuffer.scale(2.0 * a_exp, {3095, 3245});

                pbuffer.scale(2.0 * a_exp, {3545, 3755});

                pbuffer.scale(pfactors, 0, 2.0, {4220, 4310});

                pbuffer.scale(pfactors, 0, 2.0, {4400, 4550});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {4700, 4925});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {5150, 5465});

                pbuffer.scale(pfactors, 0, 2.0, {5780, 5906});

                pbuffer.scale(pfactors, 0, 2.0, {5906, 6116});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {6116, 6431});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {6431, 6872});

                t2cfunc::reduce(cbuffer, 2172, pbuffer, 451, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2190, pbuffer, 541, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2220, pbuffer, 691, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2265, pbuffer, 916, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2328, pbuffer, 1297, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2364, pbuffer, 1441, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2424, pbuffer, 1681, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2514, pbuffer, 2041, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2640, pbuffer, 2615, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2700, pbuffer, 2795, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2800, pbuffer, 3095, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2950, pbuffer, 3545, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3160, pbuffer, 4220, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3250, pbuffer, 4400, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3400, pbuffer, 4700, 225, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3625, pbuffer, 5150, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3940, pbuffer, 5780, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4066, pbuffer, 5906, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4276, pbuffer, 6116, 315, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 4591, pbuffer, 6431, 441, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 279, cbuffer, 0, 18, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1647, cbuffer, 48, 84, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 4197, cbuffer, 144, 204, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 7176, cbuffer, 1292, 1310, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 8544, cbuffer, 1340, 1376, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 11094, cbuffer, 1436, 1496, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 15189, cbuffer, 1596, 1686, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 21192, cbuffer, 1836, 1962, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 0, cbuffer, 304, 322, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 54, cbuffer, 322, 352, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 144, cbuffer, 352, 397, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 333, cbuffer, 0, 0, 54, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 495, cbuffer, 18, 54, 144, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 765, 279, 333, 495, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1089, cbuffer, 460, 496, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 1197, cbuffer, 496, 556, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 1377, cbuffer, 556, 646, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 1755, cbuffer, 48, 1089, 1197, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 2079, cbuffer, 84, 1197, 1377, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 2619, 1647, 1755, 2079, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3267, cbuffer, 772, 832, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 3447, cbuffer, 832, 932, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 3747, cbuffer, 932, 1082, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 4377, cbuffer, 144, 3267, 3447, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 4917, cbuffer, 204, 3447, 3747, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 5817, 4197, 4377, 4917, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6897, cbuffer, 2172, 2190, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 6951, cbuffer, 2190, 2220, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 7041, cbuffer, 2220, 2265, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 7230, cbuffer, 1292, 6897, 6951, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 7392, cbuffer, 1310, 6951, 7041, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 7662, 7176, 7230, 7392, cfactors, 6, 0, 1);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 7986, cbuffer, 2328, 2364, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 8094, cbuffer, 2364, 2424, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 8274, cbuffer, 2424, 2514, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 8652, cbuffer, 1340, 7986, 8094, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 8976, cbuffer, 1376, 8094, 8274, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 9516, 8544, 8652, 8976, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 10164, cbuffer, 2640, 2700, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 10344, cbuffer, 2700, 2800, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 10644, cbuffer, 2800, 2950, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 11274, cbuffer, 1436, 10164, 10344, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 11814, cbuffer, 1496, 10344, 10644, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 12714, 11094, 11274, 11814, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 13794, cbuffer, 3160, 3250, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 14064, cbuffer, 3250, 3400, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 14514, cbuffer, 3400, 3625, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 15459, cbuffer, 1596, 13794, 14064, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 16269, cbuffer, 1686, 14064, 14514, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 17619, 15189, 15459, 16269, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 19239, cbuffer, 3940, 4066, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsf(ckbuffer, 19617, cbuffer, 4066, 4276, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsg(ckbuffer, 20247, cbuffer, 4276, 4591, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpd(ckbuffer, 21570, cbuffer, 1836, 19239, 19617, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpf(ckbuffer, 22704, cbuffer, 1962, 19617, 20247, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxdd(ckbuffer, 24594, 21192, 21570, 22704, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 0, ckbuffer, 765, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 75, ckbuffer, 873, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 150, ckbuffer, 981, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 900, ckbuffer, 2619, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 1050, ckbuffer, 2835, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 1200, ckbuffer, 3051, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 2700, ckbuffer, 5817, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 2950, ckbuffer, 6177, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 3200, ckbuffer, 6537, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 44175, ckbuffer, 7662, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 44250, ckbuffer, 7770, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 44325, ckbuffer, 7878, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 44400, ckbuffer, 9516, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 44550, ckbuffer, 9732, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 44700, ckbuffer, 9948, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 44850, ckbuffer, 12714, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 45100, ckbuffer, 13074, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 45350, ckbuffer, 13434, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 45600, ckbuffer, 17619, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 45975, ckbuffer, 18159, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 46350, ckbuffer, 18699, 0, 4);

            t4cfunc::ket_transform<2, 2>(skbuffer, 46725, ckbuffer, 24594, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 47250, ckbuffer, 25350, 0, 5);

            t4cfunc::ket_transform<2, 2>(skbuffer, 47775, ckbuffer, 26106, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 9075, 0, 900, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 9300, 75, 1050, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 9525, 150, 1200, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 11775, 900, 2700, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 12225, 1050, 2950, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 12675, 1200, 3200, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 23925, 9075, 11775, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 24375, 9300, 12225, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 24825, 9525, 12675, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_spxx(skbuffer, 225, 44175, 44400, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 1350, 44400, 44850, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 3450, 44850, 45600, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 5700, 45600, 46725, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ppxx(skbuffer, 9750, 0, 225, 1350, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 13125, 900, 1350, 3450, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 17175, 2700, 3450, 5700, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dpxx(skbuffer, 25275, 9075, 9750, 13125, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 29325, 11775, 13125, 17175, r_ab, 2, 2);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fpxx(skbuffer, 37425, 23925, 25275, 29325, r_ab, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 0, skbuffer, 37425, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 525, skbuffer, 38175, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1050, skbuffer, 38925, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 1575, skbuffer, 39675, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2100, skbuffer, 40425, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 2625, skbuffer, 41175, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 3150, skbuffer, 41925, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 3675, skbuffer, 42675, 2, 2);

            t4cfunc::bra_transform<3, 1>(sbuffer, 4200, skbuffer, 43425, 2, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 1, 2, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFPDD_hpp */
