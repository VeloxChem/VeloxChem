#ifndef ElectronRepulsionRecSSGG_hpp
#define ElectronRepulsionRecSSGG_hpp

#include <array>

#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXDH.hpp"
#include "ElectronRepulsionContrRecXXDI.hpp"
#include "ElectronRepulsionContrRecXXFG.hpp"
#include "ElectronRepulsionContrRecXXFH.hpp"
#include "ElectronRepulsionContrRecXXGG.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
#include "ElectronRepulsionContrRecXXPI.hpp"
#include "ElectronRepulsionContrRecXXPK.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSH.hpp"
#include "ElectronRepulsionPrimRecSSSI.hpp"
#include "ElectronRepulsionPrimRecSSSK.hpp"
#include "ElectronRepulsionPrimRecSSSL.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// Computes (SS|1/|r-r'||GG)  integrals for two GTOs pair blocks.
/// - Parameter distributor: the pointer to Fock matrix/matrices distributor.
/// - Parameter bra_gto_pair_block: the GTOs pair block on bra side.
/// - Parameter ket_gto_pair_block: the GTOs pair block on ket side.
/// - Parameter bra_indices: the range [bra_first, bra_last) of GTOs on bra side.
/// - Parameter ket_indices: the range [ket_first, ket_last) of GTOs on ket side.
template <class T>
auto
comp_electron_repulsion_ssgg(T* distributor,
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

    CSimdArray<double> prim_buffer_0_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssk(36, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssk(36, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssl(45, ket_pdim);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cart_buffer_0_sssg(15, ket_dim);

    CSimdArray<double> cart_buffer_0_sssh(21, ket_dim);

    CSimdArray<double> cart_buffer_0_sssi(28, ket_dim);

    CSimdArray<double> cart_buffer_0_sssk(36, ket_dim);

    CSimdArray<double> cart_buffer_0_sssl(45, ket_dim);

    // allocate aligned contracted integrals

    CSimdArray<double> contr_buffer_0_sspg(45, ket_dim);

    CSimdArray<double> contr_buffer_0_ssph(63, ket_dim);

    CSimdArray<double> contr_buffer_0_sspi(84, ket_dim);

    CSimdArray<double> contr_buffer_0_sspk(108, ket_dim);

    CSimdArray<double> contr_buffer_0_ssdg(90, ket_dim);

    CSimdArray<double> contr_buffer_0_ssdh(126, ket_dim);

    CSimdArray<double> contr_buffer_0_ssdi(168, ket_dim);

    CSimdArray<double> contr_buffer_0_ssfg(150, ket_dim);

    CSimdArray<double> contr_buffer_0_ssfh(210, ket_dim);

    CSimdArray<double> contr_buffer_0_ssgg(225, ket_dim);

    // allocate aligned half transformed integrals

    CSimdArray<double> ket_spher_buffer_0_ssgg(81, ket_dim);

    // allocate aligned spherical integrals

    CSimdArray<double> spher_buffer_0_ssgg(81, ket_dim);

    // allocate accumulation buffer for integrals

    CSimdArray<double> buffer(bra_dim * 81, ket_dim);

    // setup Boys fuction data

    const CBoysFunc<8> bf_table;

    CSimdArray<double> bf_args(1, ket_pdim);

    CSimdArray<double> bf_values(9, ket_pdim);

    // loop over contracted GTOs on bra side

    for (auto i = bra_indices[0]; i < bra_indices[1]; i++)
    {
        // zero integral buffers

        cart_buffer_0_sssg.zero();

        cart_buffer_0_sssh.zero();

        cart_buffer_0_sssi.zero();

        cart_buffer_0_sssk.zero();

        cart_buffer_0_sssl.zero();

        ket_spher_buffer_0_ssgg.zero();

        spher_buffer_0_ssgg.zero();

        // set up coordinates on bra side

        const auto a_x = a_coords_x[i];

        const auto a_y = a_coords_y[i];

        const auto a_z = a_coords_z[i];

        const auto b_x = b_coords_x[i];

        const auto b_y = b_coords_y[i];

        const auto b_z = b_coords_z[i];

        for (int j = 0; j < bra_npgtos; j++)
        {
            const auto a_exp = a_vec_exps[j * bra_ncgtos + i];

            const auto b_exp = b_vec_exps[j * bra_ncgtos + i];

            const auto ab_norm = ab_vec_norms[j * bra_ncgtos + i];

            const auto ab_ovl = ab_vec_ovls[j * bra_ncgtos + i];

            const auto p_x = (a_x * a_exp + b_x * b_exp) / (a_exp + b_exp);

            const auto p_y = (a_y * a_exp + b_y * b_exp) / (a_exp + b_exp);

            const auto p_z = (a_z * a_exp + b_z * b_exp) / (a_exp + b_exp);

            t4cfunc::comp_coordinates_q(q_x[0], q_y[0], q_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], c_exps[0], d_exps[0], ket_pdim);

            t4cfunc::comp_coordinates_w(w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], a_exp, b_exp, c_exps[0], d_exps[0], ket_pdim);

            t4cfunc::comp_distances_pq(pq_x[0], pq_y[0], pq_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], ket_pdim);

            t4cfunc::comp_distances_wq(wq_x[0], wq_y[0], wq_z[0], w_x[0], w_y[0], w_z[0], q_x[0], q_y[0], q_z[0], ket_pdim);

            t4cfunc::comp_distances_qd(qd_x[0], qd_y[0], qd_z[0], q_x[0], q_y[0], q_z[0], d_x[0], d_y[0], d_z[0], ket_pdim);

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

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_0_sssi, prim_buffer_0_sssg, prim_buffer_1_sssg, prim_buffer_0_sssh, prim_buffer_1_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_1_sssi, prim_buffer_1_sssg, prim_buffer_2_sssg, prim_buffer_1_sssh, prim_buffer_2_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_2_sssi, prim_buffer_2_sssg, prim_buffer_3_sssg, prim_buffer_2_sssh, prim_buffer_3_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_0_sssk, prim_buffer_0_sssh, prim_buffer_1_sssh, prim_buffer_0_sssi, prim_buffer_1_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_1_sssk, prim_buffer_1_sssh, prim_buffer_2_sssh, prim_buffer_1_sssi, prim_buffer_2_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_0_sssl, prim_buffer_0_sssi, prim_buffer_1_sssi, prim_buffer_0_sssk, prim_buffer_1_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            t2cfunc::reduce(cart_buffer_0_sssg, prim_buffer_0_sssg, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sssh, prim_buffer_0_sssh, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sssi, prim_buffer_0_sssi, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sssk, prim_buffer_0_sssk, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sssl, prim_buffer_0_sssl, ket_dim, ket_npgtos);

        }

        erirec::comp_ket_hrr_electron_repulsion_xxpg(contr_buffer_0_sspg, cart_buffer_0_sssg, cart_buffer_0_sssh, cd_x[0], cd_y[0], cd_z[0], 0, 0);

        erirec::comp_ket_hrr_electron_repulsion_xxph(contr_buffer_0_ssph, cart_buffer_0_sssh, cart_buffer_0_sssi, cd_x[0], cd_y[0], cd_z[0], 0, 0);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(contr_buffer_0_sspi, cart_buffer_0_sssi, cart_buffer_0_sssk, cd_x[0], cd_y[0], cd_z[0], 0, 0);

        erirec::comp_ket_hrr_electron_repulsion_xxpk(contr_buffer_0_sspk, cart_buffer_0_sssk, cart_buffer_0_sssl, cd_x[0], cd_y[0], cd_z[0], 0, 0);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(contr_buffer_0_ssdg, contr_buffer_0_sspg, contr_buffer_0_ssph, cd_x[0], cd_y[0], cd_z[0], 0, 0);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(contr_buffer_0_ssdh, contr_buffer_0_ssph, contr_buffer_0_sspi, cd_x[0], cd_y[0], cd_z[0], 0, 0);

        erirec::comp_ket_hrr_electron_repulsion_xxdi(contr_buffer_0_ssdi, contr_buffer_0_sspi, contr_buffer_0_sspk, cd_x[0], cd_y[0], cd_z[0], 0, 0);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(contr_buffer_0_ssfg, contr_buffer_0_ssdg, contr_buffer_0_ssdh, cd_x[0], cd_y[0], cd_z[0], 0, 0);

        erirec::comp_ket_hrr_electron_repulsion_xxfh(contr_buffer_0_ssfh, contr_buffer_0_ssdh, contr_buffer_0_ssdi, cd_x[0], cd_y[0], cd_z[0], 0, 0);

        erirec::comp_ket_hrr_electron_repulsion_xxgg(contr_buffer_0_ssgg, contr_buffer_0_ssfg, contr_buffer_0_ssfh, cd_x[0], cd_y[0], cd_z[0], 0, 0);

        t4cfunc::ket_transform<4, 4>(ket_spher_buffer_0_ssgg, contr_buffer_0_ssgg, 0, 0);

        t4cfunc::bra_transform<0, 0>(spher_buffer_0_ssgg, ket_spher_buffer_0_ssgg, 4, 4);

        t4cfunc::store_values(buffer, spher_buffer_0_ssgg, 81 * (i - bra_indices[0]));
    }

    distributor->distribute(buffer, a_indices, b_indices, c_indices, d_indices, 0, 0, 4, 4, bra_indices, ket_indices);
}

} // erirec namespace

#endif /* ElectronRepulsionRecSSGG_hpp */
