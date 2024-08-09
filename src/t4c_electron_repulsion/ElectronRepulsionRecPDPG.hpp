#ifndef ElectronRepulsionRecPDPG_hpp
#define ElectronRepulsionRecPDPG_hpp

#include <array>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
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
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// Computes (PD|1/|r-r'||PG)  integrals for two GTOs pair blocks.
/// - Parameter distributor: the pointer to Fock matrix/matrices distributor.
/// - Parameter bra_gto_pair_block: the GTOs pair block on bra side.
/// - Parameter ket_gto_pair_block: the GTOs pair block on ket side.
/// - Parameter bra_indices: the range [bra_first, bra_last) of GTOs on bra side.
/// - Parameter ket_indices: the range [ket_first, ket_last) of GTOs on ket side.
template <class T>
auto
comp_electron_repulsion_pdpg(T* distributor,
                             const CGtoPairBlock& bra_gto_pair_block,
                             const CGtoPairBlock& ket_gto_pair_block,
                             const std::array<int, 2>& bra_indices,
                             const std::array<int, 2>& ket_indices) -> void
{
    // intialize GTOs pair data on bra side

    const auto a_coords_x = bra_gto_pair_block.bra_coordinates_x();

    const auto a_coords_y = bra_gto_pair_block.bra_coordinates_y();

    const auto a_coords_z = bra_gto_pair_block.bra_coordinates_z();

    const auto b_coords_x = bra_gto_pair_block.ket_coordinates_x();

    const auto b_coords_y = bra_gto_pair_block.ket_coordinates_y();

    const auto b_coords_z = bra_gto_pair_block.ket_coordinates_z();

    const auto a_vec_exps = bra_gto_pair_block.bra_exponents();

    const auto b_vec_exps = bra_gto_pair_block.ket_exponents();

    const auto ab_vec_norms = bra_gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = bra_gto_pair_block.overlap_factors();

    const auto a_indices = bra_gto_pair_block.bra_orbital_indices();

    const auto b_indices = bra_gto_pair_block.ket_orbital_indices();

    const auto bra_ncgtos = bra_gto_pair_block.number_of_contracted_pairs();

    const auto bra_npgtos = bra_gto_pair_block.number_of_primitive_pairs();

    // intialize GTOs data on ket side

    const auto c_coords_x = ket_gto_pair_block.bra_coordinates_x();

    const auto c_coords_y = ket_gto_pair_block.bra_coordinates_y();

    const auto c_coords_z = ket_gto_pair_block.bra_coordinates_z();

    const auto d_coords_x = ket_gto_pair_block.ket_coordinates_x();

    const auto d_coords_y = ket_gto_pair_block.ket_coordinates_y();

    const auto d_coords_z = ket_gto_pair_block.ket_coordinates_z();

    const auto c_vec_exps = ket_gto_pair_block.bra_exponents();

    const auto d_vec_exps = ket_gto_pair_block.ket_exponents();

    const auto cd_vec_norms = ket_gto_pair_block.normalization_factors();

    const auto cd_vec_ovls = ket_gto_pair_block.overlap_factors();

    const auto c_indices = ket_gto_pair_block.bra_orbital_indices();

    const auto d_indices = ket_gto_pair_block.ket_orbital_indices();

    const auto ket_npgtos = ket_gto_pair_block.number_of_primitive_pairs();

    // set up dimensions of bra and ket ranges

    const auto bra_dim = bra_indices[1] - bra_indices[0];

    const auto ket_dim = ket_indices[1] - ket_indices[0];

    // allocate aligned 2D arrays for ket side

    const auto ket_pdim = ket_dim * ket_npgtos;

    CSimdArray<double> c_x(1, ket_pdim);

    CSimdArray<double> c_y(1, ket_pdim);

    CSimdArray<double> c_z(1, ket_pdim);

    CSimdArray<double> d_x(1, ket_pdim);

    CSimdArray<double> d_y(1, ket_pdim);

    CSimdArray<double> d_z(1, ket_pdim);

    CSimdArray<double> c_exps(1, ket_pdim);

    CSimdArray<double> d_exps(1, ket_pdim);

    CSimdArray<double> cd_norms(1, ket_pdim);

    CSimdArray<double> cd_ovls(1, ket_pdim);

     // load GTOs data for ket side

    c_x.replicate(c_coords_x, ket_indices, ket_npgtos);

    c_y.replicate(c_coords_y, ket_indices, ket_npgtos);

    c_z.replicate(c_coords_z, ket_indices, ket_npgtos);

    d_x.replicate(d_coords_x, ket_indices, ket_npgtos);

    d_y.replicate(d_coords_y, ket_indices, ket_npgtos);

    d_z.replicate(d_coords_z, ket_indices, ket_npgtos);

    c_exps.load(c_vec_exps, ket_indices, ket_npgtos);

    d_exps.load(d_vec_exps, ket_indices, ket_npgtos);

    cd_norms.load(cd_vec_norms, ket_indices, ket_npgtos);

    cd_ovls.load(cd_vec_ovls, ket_indices, ket_npgtos);

    // allocate aligned coordinates of Q center

    CSimdArray<double> q_x(1, ket_pdim);

    CSimdArray<double> q_y(1, ket_pdim);

    CSimdArray<double> q_z(1, ket_pdim);

    // allocate aligned coordinates of W center

    CSimdArray<double> w_x(1, ket_pdim);

    CSimdArray<double> w_y(1, ket_pdim);

    CSimdArray<double> w_z(1, ket_pdim);

    // allocate aligned distances R(PQ) = P - Q

    CSimdArray<double> pq_x(1, ket_pdim);

    CSimdArray<double> pq_y(1, ket_pdim);

    CSimdArray<double> pq_z(1, ket_pdim);

    // allocate aligned distances R(QD) = Q - D

    CSimdArray<double> qd_x(1, ket_pdim);

    CSimdArray<double> qd_y(1, ket_pdim);

    CSimdArray<double> qd_z(1, ket_pdim);

    // allocate aligned distances R(WQ) = W - Q

    CSimdArray<double> wq_x(1, ket_pdim);

    CSimdArray<double> wq_y(1, ket_pdim);

    CSimdArray<double> wq_z(1, ket_pdim);

    // allocate aligned distances R(WP) = W - P

    CSimdArray<double> wp_x(1, ket_pdim);

    CSimdArray<double> wp_y(1, ket_pdim);

    CSimdArray<double> wp_z(1, ket_pdim);

    // allocate combined overlap factor

    CSimdArray<double> fss_abcd(1, ket_pdim);

    // allocate and initialize aligned distances R(CD) = C - D

    CSimdArray<double> cd_x(1, ket_dim);

    CSimdArray<double> cd_y(1, ket_dim);

    CSimdArray<double> cd_z(1, ket_dim);

    t4cfunc::comp_distances_cd(cd_x[0], cd_y[0], cd_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], ket_dim);

    // allocate aligned primitive integrals

    CSimdArray<double> prim_buffer_0_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_1_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_2_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_3_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_4_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_5_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_6_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_7_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_8_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsf(30, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsf(30, ket_pdim);

    CSimdArray<double> prim_buffer_0_spsg(45, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsg(45, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsg(45, ket_pdim);

    CSimdArray<double> prim_buffer_0_spsh(63, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsh(63, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsh(63, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsf(60, ket_pdim);

    CSimdArray<double> prim_buffer_0_sdsg(90, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsg(90, ket_pdim);

    CSimdArray<double> prim_buffer_0_sdsh(126, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsh(126, ket_pdim);

    CSimdArray<double> prim_buffer_0_sfsg(150, ket_pdim);

    CSimdArray<double> prim_buffer_0_sfsh(210, ket_pdim);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cart_buffer_0_sdsg(90, ket_dim);

    CSimdArray<double> cart_buffer_0_sdsh(126, ket_dim);

    CSimdArray<double> cart_buffer_0_sfsg(150, ket_dim);

    CSimdArray<double> cart_buffer_0_sfsh(210, ket_dim);

    // allocate aligned contracted integrals

    CSimdArray<double> contr_buffer_0_sdpg(270, ket_dim);

    CSimdArray<double> contr_buffer_0_sfpg(450, ket_dim);

    // allocate aligned half transformed integrals

    CSimdArray<double> ket_spher_buffer_0_sdpg(162, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_sfpg(270, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_pdpg(486, ket_dim);

    // allocate aligned spherical integrals

    CSimdArray<double> spher_buffer_0_pdpg(405, ket_dim);

    // allocate accumulation buffer for integrals

    CSimdArray<double> buffer(bra_dim * 405, ket_dim);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_args(1, ket_pdim);

    CSimdArray<double> bf_values(9, ket_pdim);

    // loop over contracted GTOs on bra side

    for (auto i = bra_indices[0]; i < bra_indices[1]; i++)
    {
        // zero integral buffers

        cart_buffer_0_sdsg.zero();

        cart_buffer_0_sdsh.zero();

        cart_buffer_0_sfsg.zero();

        cart_buffer_0_sfsh.zero();

        ket_spher_buffer_0_sdpg.zero();

        ket_spher_buffer_0_sfpg.zero();

        ket_spher_buffer_0_pdpg.zero();

        spher_buffer_0_pdpg.zero();

        // set up coordinates on bra side

        const auto a_x = a_coords_x[i];

        const auto a_y = a_coords_y[i];

        const auto a_z = a_coords_z[i];

        const auto b_x = b_coords_x[i];

        const auto b_y = b_coords_y[i];

        const auto b_z = b_coords_z[i];

        // set up distances on bra side

        const auto ab_x = a_x - b_x;

        const auto ab_y = a_y - b_y;

        const auto ab_z = a_z - b_z;

        for (int j = 0; j < bra_npgtos; j++)
        {
            const auto a_exp = a_vec_exps[j * bra_ncgtos + i];

            const auto b_exp = b_vec_exps[j * bra_ncgtos + i];

            const auto ab_norm = ab_vec_norms[j * bra_ncgtos + i];

            const auto ab_ovl = ab_vec_ovls[j * bra_ncgtos + i];

            const auto p_x = (a_x * a_exp + b_x * b_exp) / (a_exp + b_exp);

            const auto p_y = (a_y * a_exp + b_y * b_exp) / (a_exp + b_exp);

            const auto p_z = (a_z * a_exp + b_z * b_exp) / (a_exp + b_exp);

            const auto pb_x = p_x - b_x;

            const auto pb_y = p_y - b_y;

            const auto pb_z = p_z - b_z;

            t4cfunc::comp_coordinates_q(q_x[0], q_y[0], q_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], c_exps[0], d_exps[0], ket_pdim);

            t4cfunc::comp_coordinates_w(w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], a_exp, b_exp, c_exps[0], d_exps[0], ket_pdim);

            t4cfunc::comp_distances_pq(pq_x[0], pq_y[0], pq_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], ket_pdim);

            t4cfunc::comp_distances_wq(wq_x[0], wq_y[0], wq_z[0], w_x[0], w_y[0], w_z[0], q_x[0], q_y[0], q_z[0], ket_pdim);

            t4cfunc::comp_distances_qd(qd_x[0], qd_y[0], qd_z[0], q_x[0], q_y[0], q_z[0], d_x[0], d_y[0], d_z[0], ket_pdim);

            t4cfunc::comp_distances_wp(wp_x[0], wp_y[0], wp_z[0], w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, ket_pdim);

            t4cfunc::comp_boys_args(bf_args, pq_x[0], pq_y[0], pq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            bf_table.compute(bf_values, bf_args);

            t4cfunc::comp_ovl_factors(fss_abcd, ab_ovl, cd_ovls[0], ab_norm, cd_norms[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_0_ssss, fss_abcd[0], bf_values[0]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_1_ssss, fss_abcd[0], bf_values[1]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_2_ssss, fss_abcd[0], bf_values[2]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_3_ssss, fss_abcd[0], bf_values[3]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_4_ssss, fss_abcd[0], bf_values[4]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_5_ssss, fss_abcd[0], bf_values[5]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_6_ssss, fss_abcd[0], bf_values[6]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_7_ssss, fss_abcd[0], bf_values[7]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_8_ssss, fss_abcd[0], bf_values[8]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_0_sssp, prim_buffer_0_ssss, prim_buffer_1_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_1_sssp, prim_buffer_1_ssss, prim_buffer_2_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_2_sssp, prim_buffer_2_ssss, prim_buffer_3_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_3_sssp, prim_buffer_3_ssss, prim_buffer_4_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_4_sssp, prim_buffer_4_ssss, prim_buffer_5_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_5_sssp, prim_buffer_5_ssss, prim_buffer_6_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_6_sssp, prim_buffer_6_ssss, prim_buffer_7_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_7_sssp, prim_buffer_7_ssss, prim_buffer_8_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_0_sssd, prim_buffer_0_ssss, prim_buffer_1_ssss, prim_buffer_0_sssp, prim_buffer_1_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_1_sssd, prim_buffer_1_ssss, prim_buffer_2_ssss, prim_buffer_1_sssp, prim_buffer_2_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_2_sssd, prim_buffer_2_ssss, prim_buffer_3_ssss, prim_buffer_2_sssp, prim_buffer_3_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_3_sssd, prim_buffer_3_ssss, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_4_sssd, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_5_sssd, prim_buffer_5_ssss, prim_buffer_6_ssss, prim_buffer_5_sssp, prim_buffer_6_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_6_sssd, prim_buffer_6_ssss, prim_buffer_7_ssss, prim_buffer_6_sssp, prim_buffer_7_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_0_sssf, prim_buffer_0_sssp, prim_buffer_1_sssp, prim_buffer_0_sssd, prim_buffer_1_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_1_sssf, prim_buffer_1_sssp, prim_buffer_2_sssp, prim_buffer_1_sssd, prim_buffer_2_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_2_sssf, prim_buffer_2_sssp, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_3_sssf, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_3_sssd, prim_buffer_4_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_4_sssf, prim_buffer_4_sssp, prim_buffer_5_sssp, prim_buffer_4_sssd, prim_buffer_5_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_5_sssf, prim_buffer_5_sssp, prim_buffer_6_sssp, prim_buffer_5_sssd, prim_buffer_6_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_0_sssg, prim_buffer_0_sssd, prim_buffer_1_sssd, prim_buffer_0_sssf, prim_buffer_1_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_1_sssg, prim_buffer_1_sssd, prim_buffer_2_sssd, prim_buffer_1_sssf, prim_buffer_2_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_2_sssg, prim_buffer_2_sssd, prim_buffer_3_sssd, prim_buffer_2_sssf, prim_buffer_3_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_3_sssg, prim_buffer_3_sssd, prim_buffer_4_sssd, prim_buffer_3_sssf, prim_buffer_4_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_4_sssg, prim_buffer_4_sssd, prim_buffer_5_sssd, prim_buffer_4_sssf, prim_buffer_5_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_0_sssh, prim_buffer_0_sssf, prim_buffer_1_sssf, prim_buffer_0_sssg, prim_buffer_1_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_1_sssh, prim_buffer_1_sssf, prim_buffer_2_sssf, prim_buffer_1_sssg, prim_buffer_2_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_2_sssh, prim_buffer_2_sssf, prim_buffer_3_sssf, prim_buffer_2_sssg, prim_buffer_3_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_3_sssh, prim_buffer_3_sssf, prim_buffer_4_sssf, prim_buffer_3_sssg, prim_buffer_4_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_2_spsd, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_1_spsf, prim_buffer_2_sssd, prim_buffer_1_sssf, prim_buffer_2_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_2_spsf, prim_buffer_3_sssd, prim_buffer_2_sssf, prim_buffer_3_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_0_spsg, prim_buffer_1_sssf, prim_buffer_0_sssg, prim_buffer_1_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_1_spsg, prim_buffer_2_sssf, prim_buffer_1_sssg, prim_buffer_2_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_2_spsg, prim_buffer_3_sssf, prim_buffer_2_sssg, prim_buffer_3_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_0_spsh, prim_buffer_1_sssg, prim_buffer_0_sssh, prim_buffer_1_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_1_spsh, prim_buffer_2_sssg, prim_buffer_1_sssh, prim_buffer_2_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_2_spsh, prim_buffer_3_sssg, prim_buffer_2_sssh, prim_buffer_3_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_1_sdsf, prim_buffer_1_sssf, prim_buffer_2_sssf, prim_buffer_2_spsd, prim_buffer_1_spsf, prim_buffer_2_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_0_sdsg, prim_buffer_0_sssg, prim_buffer_1_sssg, prim_buffer_1_spsf, prim_buffer_0_spsg, prim_buffer_1_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_1_sdsg, prim_buffer_1_sssg, prim_buffer_2_sssg, prim_buffer_2_spsf, prim_buffer_1_spsg, prim_buffer_2_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_0_sdsh, prim_buffer_0_sssh, prim_buffer_1_sssh, prim_buffer_1_spsg, prim_buffer_0_spsh, prim_buffer_1_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_1_sdsh, prim_buffer_1_sssh, prim_buffer_2_sssh, prim_buffer_2_spsg, prim_buffer_1_spsh, prim_buffer_2_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_0_sfsg, prim_buffer_0_spsg, prim_buffer_1_spsg, prim_buffer_1_sdsf, prim_buffer_0_sdsg, prim_buffer_1_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_0_sfsh, prim_buffer_0_spsh, prim_buffer_1_spsh, prim_buffer_1_sdsg, prim_buffer_0_sdsh, prim_buffer_1_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            t2cfunc::reduce(cart_buffer_0_sdsg, prim_buffer_0_sdsg, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sdsh, prim_buffer_0_sdsh, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sfsg, prim_buffer_0_sfsg, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sfsh, prim_buffer_0_sfsh, ket_dim, ket_npgtos);

        }

        erirec::comp_ket_hrr_electron_repulsion_xxpg(contr_buffer_0_sdpg, cart_buffer_0_sdsg, cart_buffer_0_sdsh, cd_x[0], cd_y[0], cd_z[0], 0, 2);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(contr_buffer_0_sfpg, cart_buffer_0_sfsg, cart_buffer_0_sfsh, cd_x[0], cd_y[0], cd_z[0], 0, 3);

        t4cfunc::ket_transform<1, 4>(ket_spher_buffer_0_sdpg, contr_buffer_0_sdpg, 0, 2);

        t4cfunc::ket_transform<1, 4>(ket_spher_buffer_0_sfpg, contr_buffer_0_sfpg, 0, 3);

        erirec::comp_bra_hrr_electron_repulsion_pdxx(ket_spher_buffer_0_pdpg, ket_spher_buffer_0_sdpg, ket_spher_buffer_0_sfpg, ab_x, ab_y, ab_z, 1, 4);

        t4cfunc::bra_transform<1, 2>(spher_buffer_0_pdpg, ket_spher_buffer_0_pdpg, 1, 4);

        t4cfunc::store_values(buffer, spher_buffer_0_pdpg, 405 * (i - bra_indices[0]));
    }

    distributor->distribute(buffer, a_indices, b_indices, c_indices, d_indices, 1, 2, 1, 4, bra_indices, ket_indices);
}

} // erirec namespace

#endif /* ElectronRepulsionRecPDPG_hpp */
