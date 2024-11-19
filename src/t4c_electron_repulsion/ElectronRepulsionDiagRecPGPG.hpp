#ifndef ElectronRepulsionDiagRecPGPG_hpp
#define ElectronRepulsionDiagRecPGPG_hpp

#include <cstddef>
#include <utility>
#include <vector>

#include "BoysFunc.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
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
#include "GtoPairBlock.hpp"
#include "SimdArray.hpp"
#include "T2CUtils.hpp"
#include "T4CUtils.hpp"

namespace erirec {  // erirec namespace

/// @brief Computes (PG|1/|r-r'||PG)  integrals for GTOs pair block.
/// @param distributor The pointer to screening data distributor.
/// @param gto_pair_block The GTOs pair block.
/// @param gto_indices The range [gto_first, gto_last) of GTOs on bra and ket sides.
template <class T>
auto
comp_diag_electron_repulsion_pgpg(T& distributor, const CGtoPairBlock& gto_pair_block, const std::pair<size_t, size_t>& gto_indices) -> void
{
    // intialize GTOs pair data

    const auto a_coords = gto_pair_block.bra_coordinates();

    const auto b_coords = gto_pair_block.ket_coordinates();

    const auto a_vec_exps = gto_pair_block.bra_exponents();

    const auto b_vec_exps = gto_pair_block.ket_exponents();

    const auto ab_vec_norms = gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = gto_pair_block.overlap_factors();

    const auto a_indices = gto_pair_block.bra_orbital_indices();

    const auto b_indices = gto_pair_block.ket_orbital_indices();

    const auto ncgtos = gto_pair_block.number_of_contracted_pairs();

    const auto npgtos = gto_pair_block.number_of_primitive_pairs();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> pfactors(29, npgtos);

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(5601, npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1296, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(1620, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(2187, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(729, 1);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, npgtos);

    // allocate aligned array to store max. integral values

    const auto gto_dims = gto_indices.second - gto_indices.first;

    std::vector<double> max_values(gto_dims, 0.0);

    // loop over contracted GTOs on bra and ket sides

    for (auto i = gto_indices.first; i < gto_indices.second; i++)
    {
        // set up indices on ket side

        auto ket_range = std::pair<size_t, size_t>{i, i + 1};

        pfactors.load(a_vec_exps, ket_range, 0, npgtos);

        pfactors.load(b_vec_exps, ket_range, 1, npgtos);

        pfactors.load(ab_vec_ovls, ket_range, 2, npgtos);

        pfactors.load(ab_vec_norms, ket_range, 3, npgtos);

        pfactors.replicate_points(a_coords, ket_range, 4, npgtos);

        pfactors.replicate_points(b_coords, ket_range, 7, npgtos);

        cfactors.replicate_points(a_coords, ket_range, 0, 1);

        cfactors.replicate_points(b_coords, ket_range, 3, 1);

        t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        ckbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // zero integral buffers

        cbuffer.zero();

        ckbuffer.zero();

        skbuffer.zero();

        sbuffer.zero();

        // set up coordinates on bra side

        const auto r_a = a_coords[i];

        const auto r_b = b_coords[i];

        const auto a_xyz = r_a.coordinates();

        const auto b_xyz = r_b.coordinates();

        const auto r_ab = TPoint<double>({a_xyz[0] - b_xyz[0], a_xyz[1] - b_xyz[1], a_xyz[2] - b_xyz[2]});

        for (int j = 0; j < npgtos; j++)
        {
            const auto a_exp = a_vec_exps[j * ncgtos + i];

            const auto b_exp = b_vec_exps[j * ncgtos + i];

            const auto ab_norm = ab_vec_norms[j * ncgtos + i];

            const auto ab_ovl = ab_vec_ovls[j * ncgtos + i];

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

            t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp);

            bf_table.compute(bf_data, 0, 11);

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

            erirec::comp_prim_electron_repulsion_spss(pbuffer, 406, 4, 5, pfactors, 26, r_pb);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 409, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 418, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 427, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 445, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 463, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 481, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 511, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 541, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 571, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 601, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 646, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 691, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 736, 135, 220, 235, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 781, 145, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 826, 190, 280, 301, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 889, 205, 301, 322, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 952, 220, 322, 343, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1015, 235, 343, 364, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1078, 250, 364, 385, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1141, 20, 23, 406, 409, 418, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1159, 53, 59, 409, 427, 445, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1195, 59, 65, 418, 445, 463, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1231, 105, 115, 427, 481, 511, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1291, 115, 125, 445, 511, 541, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1351, 125, 135, 463, 541, 571, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1411, 175, 190, 481, 601, 646, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1501, 190, 205, 511, 646, 691, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1591, 205, 220, 541, 691, 736, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1681, 220, 235, 571, 736, 781, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1771, 280, 301, 646, 826, 889, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1897, 301, 322, 691, 889, 952, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2023, 322, 343, 736, 952, 1015, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2149, 343, 364, 781, 1015, 1078, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2275, 427, 445, 1141, 1159, 1195, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2335, 481, 511, 1159, 1231, 1291, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2435, 511, 541, 1195, 1291, 1351, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2535, 601, 646, 1231, 1411, 1501, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2685, 646, 691, 1291, 1501, 1591, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2835, 691, 736, 1351, 1591, 1681, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 2985, 826, 889, 1501, 1771, 1897, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3195, 889, 952, 1591, 1897, 2023, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3405, 952, 1015, 1681, 2023, 2149, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3615, 1231, 1291, 2275, 2335, 2435, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3765, 1411, 1501, 2335, 2535, 2685, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3990, 1501, 1591, 2435, 2685, 2835, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 4215, 1771, 1897, 2685, 2985, 3195, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 4530, 1897, 2023, 2835, 3195, 3405, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsg(pbuffer, 4845, 2535, 2685, 3615, 3765, 3990, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsh(pbuffer, 5160, 2985, 3195, 3990, 4215, 4530, pfactors, 26, r_pb, a_exp, b_exp);

            t2cfunc::reduce(cbuffer, 0, pbuffer, 3765, 225, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 225, pbuffer, 4215, 315, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 540, pbuffer, 4845, 315, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 855, pbuffer, 5160, 441, ket_width, npgtos);
        }

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 0, cbuffer, 0, 225, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 675, cbuffer, 540, 855, cfactors, 6, 0, 5);

        t4cfunc::ket_transform<1, 4>(skbuffer, 0, ckbuffer, 0, 0, 4);

        t4cfunc::ket_transform<1, 4>(skbuffer, 405, ckbuffer, 675, 0, 5);

        erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 972, 0, 405, r_ab, 1, 4);

        t4cfunc::bra_transform<1, 4>(sbuffer, 0, skbuffer, 972, 1, 4);

        t4cfunc::update_max_values(max_values, sbuffer, i - gto_indices.first);
    }

    distributor.distribute(max_values, gto_indices);
}

template <class T>
auto
comp_diag_electron_repulsion_gpgp(T& distributor, const CGtoPairBlock& gto_pair_block, const std::pair<size_t, size_t>& gto_indices) -> void
{
    // intialize GTOs pair data

    const auto a_coords = gto_pair_block.ket_coordinates();

    const auto b_coords = gto_pair_block.bra_coordinates();

    const auto a_vec_exps = gto_pair_block.ket_exponents();

    const auto b_vec_exps = gto_pair_block.bra_exponents();

    const auto ab_vec_norms = gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = gto_pair_block.overlap_factors();

    const auto a_indices = gto_pair_block.ket_orbital_indices();

    const auto b_indices = gto_pair_block.bra_orbital_indices();

    const auto ncgtos = gto_pair_block.number_of_contracted_pairs();

    const auto npgtos = gto_pair_block.number_of_primitive_pairs();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> pfactors(29, npgtos);

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(5601, npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1296, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(1620, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(2187, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(729, 1);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_data(12, npgtos);

    // allocate aligned array to store max. integral values

    const auto gto_dims = gto_indices.second - gto_indices.first;

    std::vector<double> max_values(gto_dims, 0.0);

    // loop over contracted GTOs on bra and ket sides

    for (auto i = gto_indices.first; i < gto_indices.second; i++)
    {
        // set up indices on ket side

        auto ket_range = std::pair<size_t, size_t>{i, i + 1};

        pfactors.load(a_vec_exps, ket_range, 0, npgtos);

        pfactors.load(b_vec_exps, ket_range, 1, npgtos);

        pfactors.load(ab_vec_ovls, ket_range, 2, npgtos);

        pfactors.load(ab_vec_norms, ket_range, 3, npgtos);

        pfactors.replicate_points(a_coords, ket_range, 4, npgtos);

        pfactors.replicate_points(b_coords, ket_range, 7, npgtos);

        cfactors.replicate_points(a_coords, ket_range, 0, 1);

        cfactors.replicate_points(b_coords, ket_range, 3, 1);

        t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        ckbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // zero integral buffers

        cbuffer.zero();

        ckbuffer.zero();

        skbuffer.zero();

        sbuffer.zero();

        // set up coordinates on bra side

        const auto r_a = a_coords[i];

        const auto r_b = b_coords[i];

        const auto a_xyz = r_a.coordinates();

        const auto b_xyz = r_b.coordinates();

        const auto r_ab = TPoint<double>({a_xyz[0] - b_xyz[0], a_xyz[1] - b_xyz[1], a_xyz[2] - b_xyz[2]});

        for (int j = 0; j < npgtos; j++)
        {
            const auto a_exp = a_vec_exps[j * ncgtos + i];

            const auto b_exp = b_vec_exps[j * ncgtos + i];

            const auto ab_norm = ab_vec_norms[j * ncgtos + i];

            const auto ab_ovl = ab_vec_ovls[j * ncgtos + i];

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

            t4cfunc::comp_boys_args(bf_data, 11, pfactors, 13, a_exp, b_exp);

            bf_table.compute(bf_data, 0, 11);

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

            erirec::comp_prim_electron_repulsion_spss(pbuffer, 406, 4, 5, pfactors, 26, r_pb);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 409, 4, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsp(pbuffer, 418, 5, 23, 26, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 427, 20, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 445, 23, 59, 65, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsd(pbuffer, 463, 26, 65, 71, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 481, 53, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 511, 59, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 541, 65, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsf(pbuffer, 571, 71, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 601, 105, 175, 190, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 646, 115, 190, 205, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 691, 125, 205, 220, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 736, 135, 220, 235, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsg(pbuffer, 781, 145, 235, 250, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 826, 190, 280, 301, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 889, 205, 301, 322, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 952, 220, 322, 343, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1015, 235, 343, 364, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_spsh(pbuffer, 1078, 250, 364, 385, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 1141, 20, 23, 406, 409, 418, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1159, 53, 59, 409, 427, 445, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 1195, 59, 65, 418, 445, 463, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1231, 105, 115, 427, 481, 511, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1291, 115, 125, 445, 511, 541, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1351, 125, 135, 463, 541, 571, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1411, 175, 190, 481, 601, 646, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1501, 190, 205, 511, 646, 691, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1591, 205, 220, 541, 691, 736, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 1681, 220, 235, 571, 736, 781, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1771, 280, 301, 646, 826, 889, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 1897, 301, 322, 691, 889, 952, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2023, 322, 343, 736, 952, 1015, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sdsh(pbuffer, 2149, 343, 364, 781, 1015, 1078, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 2275, 427, 445, 1141, 1159, 1195, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2335, 481, 511, 1159, 1231, 1291, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 2435, 511, 541, 1195, 1291, 1351, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2535, 601, 646, 1231, 1411, 1501, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2685, 646, 691, 1291, 1501, 1591, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 2835, 691, 736, 1351, 1591, 1681, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 2985, 826, 889, 1501, 1771, 1897, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3195, 889, 952, 1591, 1897, 2023, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sfsh(pbuffer, 3405, 952, 1015, 1681, 2023, 2149, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 3615, 1231, 1291, 2275, 2335, 2435, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3765, 1411, 1501, 2335, 2535, 2685, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsg(pbuffer, 3990, 1501, 1591, 2435, 2685, 2835, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 4215, 1771, 1897, 2685, 2985, 3195, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_sgsh(pbuffer, 4530, 1897, 2023, 2835, 3195, 3405, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsg(pbuffer, 4845, 2535, 2685, 3615, 3765, 3990, pfactors, 26, r_pb, a_exp, b_exp);

            erirec::comp_prim_electron_repulsion_shsh(pbuffer, 5160, 2985, 3195, 3990, 4215, 4530, pfactors, 26, r_pb, a_exp, b_exp);

            t2cfunc::reduce(cbuffer, 0, pbuffer, 3765, 225, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 225, pbuffer, 4215, 315, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 540, pbuffer, 4845, 315, ket_width, npgtos);

            t2cfunc::reduce(cbuffer, 855, pbuffer, 5160, 441, ket_width, npgtos);
        }

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 0, cbuffer, 0, 225, cfactors, 6, 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(ckbuffer, 675, cbuffer, 540, 855, cfactors, 6, 0, 5);

        t4cfunc::ket_transform<1, 4>(skbuffer, 0, ckbuffer, 0, 0, 4);

        t4cfunc::ket_transform<1, 4>(skbuffer, 405, ckbuffer, 675, 0, 5);

        erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 972, 0, 405, r_ab, 1, 4);

        t4cfunc::bra_transform<1, 4>(sbuffer, 0, skbuffer, 972, 1, 4);

        t4cfunc::update_max_values(max_values, sbuffer, i - gto_indices.first);
    }

    distributor.distribute(max_values, gto_indices);
}

}  // namespace erirec

#endif /* ElectronRepulsionDiagRecPGPG_hpp */
