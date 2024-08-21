#ifndef ElectronRepulsionGeom2000RecPDDD_hpp
#define ElectronRepulsionGeom2000RecPDDD_hpp

#include <array>

#include "ElectronRepulsionPrimRecDDDD.hpp"
#include "ElectronRepulsionPrimRecDDDP.hpp"
#include "ElectronRepulsionPrimRecDDPD.hpp"
#include "ElectronRepulsionPrimRecDPDD.hpp"
#include "ElectronRepulsionPrimRecFDDD.hpp"
#include "ElectronRepulsionPrimRecPDDD.hpp"
#include "ElectronRepulsionPrimRecPDDP.hpp"
#include "ElectronRepulsionPrimRecPDDS.hpp"
#include "ElectronRepulsionPrimRecPDPD.hpp"
#include "ElectronRepulsionPrimRecPDPP.hpp"
#include "ElectronRepulsionPrimRecPDSD.hpp"
#include "ElectronRepulsionPrimRecPPDD.hpp"
#include "ElectronRepulsionPrimRecPPDP.hpp"
#include "ElectronRepulsionPrimRecPPPD.hpp"
#include "ElectronRepulsionPrimRecPSDD.hpp"
#include "ElectronRepulsionPrimRecSDDD.hpp"
#include "ElectronRepulsionPrimRecSDDP.hpp"
#include "ElectronRepulsionPrimRecSDDS.hpp"
#include "ElectronRepulsionPrimRecSDPD.hpp"
#include "ElectronRepulsionPrimRecSDPP.hpp"
#include "ElectronRepulsionPrimRecSDPS.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSPDD.hpp"
#include "ElectronRepulsionPrimRecSPDP.hpp"
#include "ElectronRepulsionPrimRecSPDS.hpp"
#include "ElectronRepulsionPrimRecSPPD.hpp"
#include "ElectronRepulsionPrimRecSPPP.hpp"
#include "ElectronRepulsionPrimRecSPPS.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSDD.hpp"
#include "ElectronRepulsionPrimRecSSDP.hpp"
#include "ElectronRepulsionPrimRecSSDS.hpp"
#include "ElectronRepulsionPrimRecSSPD.hpp"
#include "ElectronRepulsionPrimRecSSPP.hpp"
#include "ElectronRepulsionPrimRecSSPS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "GeomDeriv2000OfScalarForPDDD.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// Computes (PD|1/|r-r'||DD)  integrals for two GTOs pair blocks.
/// - Parameter distributor: the pointer to Fock matrix/matrices distributor.
/// - Parameter bra_gto_pair_block: the GTOs pair block on bra side.
/// - Parameter ket_gto_pair_block: the GTOs pair block on ket side.
/// - Parameter bra_indices: the range [bra_first, bra_last) of GTOs on bra side.
/// - Parameter ket_indices: the range [ket_first, ket_last) of GTOs on ket side.
template <class T>
auto
comp_electron_repulsion_geom2000_pddd(T* distributor,
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

    // allocate aligned distances R(QC) = Q - C

    CSimdArray<double> qc_x(1, ket_pdim);

    CSimdArray<double> qc_y(1, ket_pdim);

    CSimdArray<double> qc_z(1, ket_pdim);

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

    CSimdArray<double> prim_buffer_9_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_8_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_1_ssps(3, ket_pdim);

    CSimdArray<double> prim_buffer_2_ssps(3, ket_pdim);

    CSimdArray<double> prim_buffer_3_ssps(3, ket_pdim);

    CSimdArray<double> prim_buffer_4_ssps(3, ket_pdim);

    CSimdArray<double> prim_buffer_5_ssps(3, ket_pdim);

    CSimdArray<double> prim_buffer_6_ssps(3, ket_pdim);

    CSimdArray<double> prim_buffer_0_sspp(9, ket_pdim);

    CSimdArray<double> prim_buffer_1_sspp(9, ket_pdim);

    CSimdArray<double> prim_buffer_2_sspp(9, ket_pdim);

    CSimdArray<double> prim_buffer_3_sspp(9, ket_pdim);

    CSimdArray<double> prim_buffer_4_sspp(9, ket_pdim);

    CSimdArray<double> prim_buffer_5_sspp(9, ket_pdim);

    CSimdArray<double> prim_buffer_6_sspp(9, ket_pdim);

    CSimdArray<double> prim_buffer_0_sspd(18, ket_pdim);

    CSimdArray<double> prim_buffer_1_sspd(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_sspd(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_sspd(18, ket_pdim);

    CSimdArray<double> prim_buffer_4_sspd(18, ket_pdim);

    CSimdArray<double> prim_buffer_5_sspd(18, ket_pdim);

    CSimdArray<double> prim_buffer_6_sspd(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_ssds(6, ket_pdim);

    CSimdArray<double> prim_buffer_3_ssds(6, ket_pdim);

    CSimdArray<double> prim_buffer_4_ssds(6, ket_pdim);

    CSimdArray<double> prim_buffer_5_ssds(6, ket_pdim);

    CSimdArray<double> prim_buffer_1_ssdp(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_ssdp(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_ssdp(18, ket_pdim);

    CSimdArray<double> prim_buffer_4_ssdp(18, ket_pdim);

    CSimdArray<double> prim_buffer_5_ssdp(18, ket_pdim);

    CSimdArray<double> prim_buffer_0_ssdd(36, ket_pdim);

    CSimdArray<double> prim_buffer_1_ssdd(36, ket_pdim);

    CSimdArray<double> prim_buffer_2_ssdd(36, ket_pdim);

    CSimdArray<double> prim_buffer_3_ssdd(36, ket_pdim);

    CSimdArray<double> prim_buffer_4_ssdd(36, ket_pdim);

    CSimdArray<double> prim_buffer_5_ssdd(36, ket_pdim);

    CSimdArray<double> prim_buffer_4_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_spps(9, ket_pdim);

    CSimdArray<double> prim_buffer_4_spps(9, ket_pdim);

    CSimdArray<double> prim_buffer_2_sppp(27, ket_pdim);

    CSimdArray<double> prim_buffer_3_sppp(27, ket_pdim);

    CSimdArray<double> prim_buffer_4_sppp(27, ket_pdim);

    CSimdArray<double> prim_buffer_1_sppd(54, ket_pdim);

    CSimdArray<double> prim_buffer_2_sppd(54, ket_pdim);

    CSimdArray<double> prim_buffer_3_sppd(54, ket_pdim);

    CSimdArray<double> prim_buffer_4_sppd(54, ket_pdim);

    CSimdArray<double> prim_buffer_2_spds(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_spds(18, ket_pdim);

    CSimdArray<double> prim_buffer_4_spds(18, ket_pdim);

    CSimdArray<double> prim_buffer_1_spdp(54, ket_pdim);

    CSimdArray<double> prim_buffer_2_spdp(54, ket_pdim);

    CSimdArray<double> prim_buffer_3_spdp(54, ket_pdim);

    CSimdArray<double> prim_buffer_4_spdp(54, ket_pdim);

    CSimdArray<double> prim_buffer_0_spdd(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_spdd(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_spdd(108, ket_pdim);

    CSimdArray<double> prim_buffer_3_spdd(108, ket_pdim);

    CSimdArray<double> prim_buffer_4_spdd(108, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdps(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdpp(54, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdpp(54, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdpd(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdpd(108, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdpd(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdds(36, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdds(36, ket_pdim);

    CSimdArray<double> prim_buffer_1_sddp(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_sddp(108, ket_pdim);

    CSimdArray<double> prim_buffer_3_sddp(108, ket_pdim);

    CSimdArray<double> prim_buffer_0_sddd(216, ket_pdim);

    CSimdArray<double> prim_buffer_1_sddd(216, ket_pdim);

    CSimdArray<double> prim_buffer_2_sddd(216, ket_pdim);

    CSimdArray<double> prim_buffer_3_sddd(216, ket_pdim);

    CSimdArray<double> prim_buffer_0_psdd(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_psdd(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_psdd(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_pppd(162, ket_pdim);

    CSimdArray<double> prim_buffer_2_pppd(162, ket_pdim);

    CSimdArray<double> prim_buffer_1_ppdp(162, ket_pdim);

    CSimdArray<double> prim_buffer_2_ppdp(162, ket_pdim);

    CSimdArray<double> prim_buffer_0_ppdd(324, ket_pdim);

    CSimdArray<double> prim_buffer_1_ppdd(324, ket_pdim);

    CSimdArray<double> prim_buffer_2_ppdd(324, ket_pdim);

    CSimdArray<double> prim_buffer_2_pdsd(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_pdpp(162, ket_pdim);

    CSimdArray<double> prim_buffer_1_pdpd(324, ket_pdim);

    CSimdArray<double> prim_buffer_2_pdpd(324, ket_pdim);

    CSimdArray<double> prim_buffer_2_pdds(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_pddp(324, ket_pdim);

    CSimdArray<double> prim_buffer_2_pddp(324, ket_pdim);

    CSimdArray<double> prim_buffer_0_pddd(648, ket_pdim);

    CSimdArray<double> prim_buffer_1_pddd(648, ket_pdim);

    CSimdArray<double> prim_buffer_2_pddd(648, ket_pdim);

    CSimdArray<double> prim_buffer_0_dpdd(648, ket_pdim);

    CSimdArray<double> prim_buffer_1_dpdd(648, ket_pdim);

    CSimdArray<double> prim_buffer_1_ddpd(648, ket_pdim);

    CSimdArray<double> prim_buffer_1_dddp(648, ket_pdim);

    CSimdArray<double> prim_buffer_0_dddd(1296, ket_pdim);

    CSimdArray<double> prim_buffer_1_dddd(1296, ket_pdim);

    CSimdArray<double> prim_buffer_0_fddd(2160, ket_pdim);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cart_buffer_0_geom2000_pddd(3888, ket_dim);

    // allocate aligned spherical integrals

    CSimdArray<double> spher_buffer_0_geom2000_pddd(2250, ket_dim);

    // setup Boys fuction data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_args(1, ket_pdim);

    CSimdArray<double> bf_values(10, ket_pdim);

    // loop over contracted GTOs on bra side

    for (auto i = bra_indices[0]; i < bra_indices[1]; i++)
    {
        cart_buffer_0_geom2000_pddd.zero();

        spher_buffer_0_geom2000_pddd.zero();

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

            const auto pa_x = p_x - a_x;

            const auto pa_y = p_y - a_y;

            const auto pa_z = p_z - a_z;

            const auto pb_x = p_x - b_x;

            const auto pb_y = p_y - b_y;

            const auto pb_z = p_z - b_z;

            t4cfunc::comp_coordinates_q(q_x[0], q_y[0], q_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], c_exps[0], d_exps[0], ket_pdim);

            t4cfunc::comp_coordinates_w(w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], a_exp, b_exp, c_exps[0], d_exps[0], ket_pdim);

            t4cfunc::comp_distances_pq(pq_x[0], pq_y[0], pq_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], ket_pdim);

            t4cfunc::comp_distances_wq(wq_x[0], wq_y[0], wq_z[0], w_x[0], w_y[0], w_z[0], q_x[0], q_y[0], q_z[0], ket_pdim);

            t4cfunc::comp_distances_qc(qc_x[0], qc_y[0], qc_z[0], q_x[0], q_y[0], q_z[0], c_x[0], c_y[0], c_z[0], ket_pdim);

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

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_9_ssss, fss_abcd[0], bf_values[9]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_0_sssp, prim_buffer_0_ssss, prim_buffer_1_ssss, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_1_sssp, prim_buffer_1_ssss, prim_buffer_2_ssss, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_2_sssp, prim_buffer_2_ssss, prim_buffer_3_ssss, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_3_sssp, prim_buffer_3_ssss, prim_buffer_4_ssss, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_4_sssp, prim_buffer_4_ssss, prim_buffer_5_ssss, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_5_sssp, prim_buffer_5_ssss, prim_buffer_6_ssss, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_6_sssp, prim_buffer_6_ssss, prim_buffer_7_ssss, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_7_sssp, prim_buffer_7_ssss, prim_buffer_8_ssss, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_8_sssp, prim_buffer_8_ssss, prim_buffer_9_ssss, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_0_sssd, prim_buffer_0_ssss, prim_buffer_1_ssss, prim_buffer_0_sssp, prim_buffer_1_sssp, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_1_sssd, prim_buffer_1_ssss, prim_buffer_2_ssss, prim_buffer_1_sssp, prim_buffer_2_sssp, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_2_sssd, prim_buffer_2_ssss, prim_buffer_3_ssss, prim_buffer_2_sssp, prim_buffer_3_sssp, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_3_sssd, prim_buffer_3_ssss, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_4_sssd, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_5_sssd, prim_buffer_5_ssss, prim_buffer_6_ssss, prim_buffer_5_sssp, prim_buffer_6_sssp, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_6_sssd, prim_buffer_6_ssss, prim_buffer_7_ssss, prim_buffer_6_sssp, prim_buffer_7_sssp, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_7_sssd, prim_buffer_7_ssss, prim_buffer_8_ssss, prim_buffer_7_sssp, prim_buffer_8_sssp, qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssps(prim_buffer_1_ssps, prim_buffer_1_ssss, prim_buffer_2_ssss, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_ssps(prim_buffer_2_ssps, prim_buffer_2_ssss, prim_buffer_3_ssss, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_ssps(prim_buffer_3_ssps, prim_buffer_3_ssss, prim_buffer_4_ssss, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_ssps(prim_buffer_4_ssps, prim_buffer_4_ssss, prim_buffer_5_ssss, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_ssps(prim_buffer_5_ssps, prim_buffer_5_ssss, prim_buffer_6_ssss, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_ssps(prim_buffer_6_ssps, prim_buffer_6_ssss, prim_buffer_7_ssss, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sspp(prim_buffer_0_sspp, prim_buffer_0_ssss, prim_buffer_1_ssss, prim_buffer_0_sssp, prim_buffer_1_sssp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspp(prim_buffer_1_sspp, prim_buffer_1_ssss, prim_buffer_2_ssss, prim_buffer_1_sssp, prim_buffer_2_sssp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspp(prim_buffer_2_sspp, prim_buffer_2_ssss, prim_buffer_3_ssss, prim_buffer_2_sssp, prim_buffer_3_sssp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspp(prim_buffer_3_sspp, prim_buffer_3_ssss, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspp(prim_buffer_4_sspp, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspp(prim_buffer_5_sspp, prim_buffer_5_ssss, prim_buffer_6_ssss, prim_buffer_5_sssp, prim_buffer_6_sssp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspp(prim_buffer_6_sspp, prim_buffer_6_ssss, prim_buffer_7_ssss, prim_buffer_6_sssp, prim_buffer_7_sssp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspd(prim_buffer_0_sspd, prim_buffer_0_sssp, prim_buffer_1_sssp, prim_buffer_0_sssd, prim_buffer_1_sssd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspd(prim_buffer_1_sspd, prim_buffer_1_sssp, prim_buffer_2_sssp, prim_buffer_1_sssd, prim_buffer_2_sssd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspd(prim_buffer_2_sspd, prim_buffer_2_sssp, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspd(prim_buffer_3_sspd, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_3_sssd, prim_buffer_4_sssd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspd(prim_buffer_4_sspd, prim_buffer_4_sssp, prim_buffer_5_sssp, prim_buffer_4_sssd, prim_buffer_5_sssd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspd(prim_buffer_5_sspd, prim_buffer_5_sssp, prim_buffer_6_sssp, prim_buffer_5_sssd, prim_buffer_6_sssd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sspd(prim_buffer_6_sspd, prim_buffer_6_sssp, prim_buffer_7_sssp, prim_buffer_6_sssd, prim_buffer_7_sssd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssds(prim_buffer_2_ssds, prim_buffer_2_ssss, prim_buffer_3_ssss, prim_buffer_2_ssps, prim_buffer_3_ssps, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssds(prim_buffer_3_ssds, prim_buffer_3_ssss, prim_buffer_4_ssss, prim_buffer_3_ssps, prim_buffer_4_ssps, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssds(prim_buffer_4_ssds, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_ssps, prim_buffer_5_ssps, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssds(prim_buffer_5_ssds, prim_buffer_5_ssss, prim_buffer_6_ssss, prim_buffer_5_ssps, prim_buffer_6_ssps, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssdp(prim_buffer_1_ssdp, prim_buffer_1_sssp, prim_buffer_2_sssp, prim_buffer_1_ssps, prim_buffer_2_ssps, prim_buffer_1_sspp, prim_buffer_2_sspp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssdp(prim_buffer_2_ssdp, prim_buffer_2_sssp, prim_buffer_3_sssp, prim_buffer_2_ssps, prim_buffer_3_ssps, prim_buffer_2_sspp, prim_buffer_3_sspp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssdp(prim_buffer_3_ssdp, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_3_ssps, prim_buffer_4_ssps, prim_buffer_3_sspp, prim_buffer_4_sspp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssdp(prim_buffer_4_ssdp, prim_buffer_4_sssp, prim_buffer_5_sssp, prim_buffer_4_ssps, prim_buffer_5_ssps, prim_buffer_4_sspp, prim_buffer_5_sspp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssdp(prim_buffer_5_ssdp, prim_buffer_5_sssp, prim_buffer_6_sssp, prim_buffer_5_ssps, prim_buffer_6_ssps, prim_buffer_5_sspp, prim_buffer_6_sspp, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssdd(prim_buffer_0_ssdd, prim_buffer_0_sssd, prim_buffer_1_sssd, prim_buffer_0_sspp, prim_buffer_1_sspp, prim_buffer_0_sspd, prim_buffer_1_sspd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssdd(prim_buffer_1_ssdd, prim_buffer_1_sssd, prim_buffer_2_sssd, prim_buffer_1_sspp, prim_buffer_2_sspp, prim_buffer_1_sspd, prim_buffer_2_sspd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssdd(prim_buffer_2_ssdd, prim_buffer_2_sssd, prim_buffer_3_sssd, prim_buffer_2_sspp, prim_buffer_3_sspp, prim_buffer_2_sspd, prim_buffer_3_sspd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssdd(prim_buffer_3_ssdd, prim_buffer_3_sssd, prim_buffer_4_sssd, prim_buffer_3_sspp, prim_buffer_4_sspp, prim_buffer_3_sspd, prim_buffer_4_sspd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssdd(prim_buffer_4_ssdd, prim_buffer_4_sssd, prim_buffer_5_sssd, prim_buffer_4_sspp, prim_buffer_5_sspp, prim_buffer_4_sspd, prim_buffer_5_sspd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssdd(prim_buffer_5_ssdd, prim_buffer_5_sssd, prim_buffer_6_sssd, prim_buffer_5_sspp, prim_buffer_6_sspp, prim_buffer_5_sspd, prim_buffer_6_sspd, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_4_spss, prim_buffer_4_ssss, prim_buffer_5_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_3_spsp, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_4_spsp, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_2_spsd, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_3_spsd, prim_buffer_4_sssp, prim_buffer_3_sssd, prim_buffer_4_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_4_spsd, prim_buffer_5_sssp, prim_buffer_4_sssd, prim_buffer_5_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spps(prim_buffer_3_spps, prim_buffer_4_ssss, prim_buffer_3_ssps, prim_buffer_4_ssps, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spps(prim_buffer_4_spps, prim_buffer_5_ssss, prim_buffer_4_ssps, prim_buffer_5_ssps, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppp(prim_buffer_2_sppp, prim_buffer_3_sssp, prim_buffer_3_ssps, prim_buffer_2_sspp, prim_buffer_3_sspp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppp(prim_buffer_3_sppp, prim_buffer_4_sssp, prim_buffer_4_ssps, prim_buffer_3_sspp, prim_buffer_4_sspp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppp(prim_buffer_4_sppp, prim_buffer_5_sssp, prim_buffer_5_ssps, prim_buffer_4_sspp, prim_buffer_5_sspp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppd(prim_buffer_1_sppd, prim_buffer_2_sssd, prim_buffer_2_sspp, prim_buffer_1_sspd, prim_buffer_2_sspd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppd(prim_buffer_2_sppd, prim_buffer_3_sssd, prim_buffer_3_sspp, prim_buffer_2_sspd, prim_buffer_3_sspd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppd(prim_buffer_3_sppd, prim_buffer_4_sssd, prim_buffer_4_sspp, prim_buffer_3_sspd, prim_buffer_4_sspd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppd(prim_buffer_4_sppd, prim_buffer_5_sssd, prim_buffer_5_sspp, prim_buffer_4_sspd, prim_buffer_5_sspd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spds(prim_buffer_2_spds, prim_buffer_3_ssps, prim_buffer_2_ssds, prim_buffer_3_ssds, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spds(prim_buffer_3_spds, prim_buffer_4_ssps, prim_buffer_3_ssds, prim_buffer_4_ssds, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spds(prim_buffer_4_spds, prim_buffer_5_ssps, prim_buffer_4_ssds, prim_buffer_5_ssds, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spdp(prim_buffer_1_spdp, prim_buffer_2_sspp, prim_buffer_2_ssds, prim_buffer_1_ssdp, prim_buffer_2_ssdp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spdp(prim_buffer_2_spdp, prim_buffer_3_sspp, prim_buffer_3_ssds, prim_buffer_2_ssdp, prim_buffer_3_ssdp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spdp(prim_buffer_3_spdp, prim_buffer_4_sspp, prim_buffer_4_ssds, prim_buffer_3_ssdp, prim_buffer_4_ssdp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spdp(prim_buffer_4_spdp, prim_buffer_5_sspp, prim_buffer_5_ssds, prim_buffer_4_ssdp, prim_buffer_5_ssdp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spdd(prim_buffer_0_spdd, prim_buffer_1_sspd, prim_buffer_1_ssdp, prim_buffer_0_ssdd, prim_buffer_1_ssdd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spdd(prim_buffer_1_spdd, prim_buffer_2_sspd, prim_buffer_2_ssdp, prim_buffer_1_ssdd, prim_buffer_2_ssdd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spdd(prim_buffer_2_spdd, prim_buffer_3_sspd, prim_buffer_3_ssdp, prim_buffer_2_ssdd, prim_buffer_3_ssdd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spdd(prim_buffer_3_spdd, prim_buffer_4_sspd, prim_buffer_4_ssdp, prim_buffer_3_ssdd, prim_buffer_4_ssdd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spdd(prim_buffer_4_spdd, prim_buffer_5_sspd, prim_buffer_5_ssdp, prim_buffer_4_ssdd, prim_buffer_5_ssdd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_3_sdsp, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_4_spss, prim_buffer_3_spsp, prim_buffer_4_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_2_sdsd, prim_buffer_2_sssd, prim_buffer_3_sssd, prim_buffer_3_spsp, prim_buffer_2_spsd, prim_buffer_3_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_3_sdsd, prim_buffer_3_sssd, prim_buffer_4_sssd, prim_buffer_4_spsp, prim_buffer_3_spsd, prim_buffer_4_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdps(prim_buffer_3_sdps, prim_buffer_3_ssps, prim_buffer_4_ssps, prim_buffer_4_spss, prim_buffer_3_spps, prim_buffer_4_spps, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpp(prim_buffer_2_sdpp, prim_buffer_2_sspp, prim_buffer_3_sspp, prim_buffer_3_spsp, prim_buffer_3_spps, prim_buffer_2_sppp, prim_buffer_3_sppp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpp(prim_buffer_3_sdpp, prim_buffer_3_sspp, prim_buffer_4_sspp, prim_buffer_4_spsp, prim_buffer_4_spps, prim_buffer_3_sppp, prim_buffer_4_sppp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpd(prim_buffer_1_sdpd, prim_buffer_1_sspd, prim_buffer_2_sspd, prim_buffer_2_spsd, prim_buffer_2_sppp, prim_buffer_1_sppd, prim_buffer_2_sppd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpd(prim_buffer_2_sdpd, prim_buffer_2_sspd, prim_buffer_3_sspd, prim_buffer_3_spsd, prim_buffer_3_sppp, prim_buffer_2_sppd, prim_buffer_3_sppd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpd(prim_buffer_3_sdpd, prim_buffer_3_sspd, prim_buffer_4_sspd, prim_buffer_4_spsd, prim_buffer_4_sppp, prim_buffer_3_sppd, prim_buffer_4_sppd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdds(prim_buffer_2_sdds, prim_buffer_2_ssds, prim_buffer_3_ssds, prim_buffer_3_spps, prim_buffer_2_spds, prim_buffer_3_spds, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdds(prim_buffer_3_sdds, prim_buffer_3_ssds, prim_buffer_4_ssds, prim_buffer_4_spps, prim_buffer_3_spds, prim_buffer_4_spds, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sddp(prim_buffer_1_sddp, prim_buffer_1_ssdp, prim_buffer_2_ssdp, prim_buffer_2_sppp, prim_buffer_2_spds, prim_buffer_1_spdp, prim_buffer_2_spdp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sddp(prim_buffer_2_sddp, prim_buffer_2_ssdp, prim_buffer_3_ssdp, prim_buffer_3_sppp, prim_buffer_3_spds, prim_buffer_2_spdp, prim_buffer_3_spdp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sddp(prim_buffer_3_sddp, prim_buffer_3_ssdp, prim_buffer_4_ssdp, prim_buffer_4_sppp, prim_buffer_4_spds, prim_buffer_3_spdp, prim_buffer_4_spdp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sddd(prim_buffer_0_sddd, prim_buffer_0_ssdd, prim_buffer_1_ssdd, prim_buffer_1_sppd, prim_buffer_1_spdp, prim_buffer_0_spdd, prim_buffer_1_spdd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sddd(prim_buffer_1_sddd, prim_buffer_1_ssdd, prim_buffer_2_ssdd, prim_buffer_2_sppd, prim_buffer_2_spdp, prim_buffer_1_spdd, prim_buffer_2_spdd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sddd(prim_buffer_2_sddd, prim_buffer_2_ssdd, prim_buffer_3_ssdd, prim_buffer_3_sppd, prim_buffer_3_spdp, prim_buffer_2_spdd, prim_buffer_3_spdd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sddd(prim_buffer_3_sddd, prim_buffer_3_ssdd, prim_buffer_4_ssdd, prim_buffer_4_sppd, prim_buffer_4_spdp, prim_buffer_3_spdd, prim_buffer_4_spdd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_psdd(prim_buffer_0_psdd, prim_buffer_1_sspd, prim_buffer_1_ssdp, prim_buffer_0_ssdd, prim_buffer_1_ssdd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_psdd(prim_buffer_1_psdd, prim_buffer_2_sspd, prim_buffer_2_ssdp, prim_buffer_1_ssdd, prim_buffer_2_ssdd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_psdd(prim_buffer_2_psdd, prim_buffer_3_sspd, prim_buffer_3_ssdp, prim_buffer_2_ssdd, prim_buffer_3_ssdd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pppd(prim_buffer_1_pppd, prim_buffer_1_sspd, prim_buffer_2_sspd, prim_buffer_2_spsd, prim_buffer_2_sppp, prim_buffer_1_sppd, prim_buffer_2_sppd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pppd(prim_buffer_2_pppd, prim_buffer_2_sspd, prim_buffer_3_sspd, prim_buffer_3_spsd, prim_buffer_3_sppp, prim_buffer_2_sppd, prim_buffer_3_sppd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppdp(prim_buffer_1_ppdp, prim_buffer_1_ssdp, prim_buffer_2_ssdp, prim_buffer_2_sppp, prim_buffer_2_spds, prim_buffer_1_spdp, prim_buffer_2_spdp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppdp(prim_buffer_2_ppdp, prim_buffer_2_ssdp, prim_buffer_3_ssdp, prim_buffer_3_sppp, prim_buffer_3_spds, prim_buffer_2_spdp, prim_buffer_3_spdp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppdd(prim_buffer_0_ppdd, prim_buffer_0_ssdd, prim_buffer_1_ssdd, prim_buffer_1_sppd, prim_buffer_1_spdp, prim_buffer_0_spdd, prim_buffer_1_spdd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppdd(prim_buffer_1_ppdd, prim_buffer_1_ssdd, prim_buffer_2_ssdd, prim_buffer_2_sppd, prim_buffer_2_spdp, prim_buffer_1_spdd, prim_buffer_2_spdd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppdd(prim_buffer_2_ppdd, prim_buffer_2_ssdd, prim_buffer_3_ssdd, prim_buffer_3_sppd, prim_buffer_3_spdp, prim_buffer_2_spdd, prim_buffer_3_spdd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdsd(prim_buffer_2_pdsd, prim_buffer_2_spsd, prim_buffer_3_spsd, prim_buffer_3_sdsp, prim_buffer_2_sdsd, prim_buffer_3_sdsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdpp(prim_buffer_2_pdpp, prim_buffer_2_sppp, prim_buffer_3_sppp, prim_buffer_3_sdsp, prim_buffer_3_sdps, prim_buffer_2_sdpp, prim_buffer_3_sdpp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdpd(prim_buffer_1_pdpd, prim_buffer_1_sppd, prim_buffer_2_sppd, prim_buffer_2_sdsd, prim_buffer_2_sdpp, prim_buffer_1_sdpd, prim_buffer_2_sdpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdpd(prim_buffer_2_pdpd, prim_buffer_2_sppd, prim_buffer_3_sppd, prim_buffer_3_sdsd, prim_buffer_3_sdpp, prim_buffer_2_sdpd, prim_buffer_3_sdpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdds(prim_buffer_2_pdds, prim_buffer_2_spds, prim_buffer_3_spds, prim_buffer_3_sdps, prim_buffer_2_sdds, prim_buffer_3_sdds, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pddp(prim_buffer_1_pddp, prim_buffer_1_spdp, prim_buffer_2_spdp, prim_buffer_2_sdpp, prim_buffer_2_sdds, prim_buffer_1_sddp, prim_buffer_2_sddp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pddp(prim_buffer_2_pddp, prim_buffer_2_spdp, prim_buffer_3_spdp, prim_buffer_3_sdpp, prim_buffer_3_sdds, prim_buffer_2_sddp, prim_buffer_3_sddp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pddd(prim_buffer_0_pddd, prim_buffer_0_spdd, prim_buffer_1_spdd, prim_buffer_1_sdpd, prim_buffer_1_sddp, prim_buffer_0_sddd, prim_buffer_1_sddd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pddd(prim_buffer_1_pddd, prim_buffer_1_spdd, prim_buffer_2_spdd, prim_buffer_2_sdpd, prim_buffer_2_sddp, prim_buffer_1_sddd, prim_buffer_2_sddd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pddd(prim_buffer_2_pddd, prim_buffer_2_spdd, prim_buffer_3_spdd, prim_buffer_3_sdpd, prim_buffer_3_sddp, prim_buffer_2_sddd, prim_buffer_3_sddd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dpdd(prim_buffer_0_dpdd, prim_buffer_0_spdd, prim_buffer_1_spdd, prim_buffer_0_psdd, prim_buffer_1_psdd, prim_buffer_1_pppd, prim_buffer_1_ppdp, prim_buffer_0_ppdd, prim_buffer_1_ppdd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dpdd(prim_buffer_1_dpdd, prim_buffer_1_spdd, prim_buffer_2_spdd, prim_buffer_1_psdd, prim_buffer_2_psdd, prim_buffer_2_pppd, prim_buffer_2_ppdp, prim_buffer_1_ppdd, prim_buffer_2_ppdd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ddpd(prim_buffer_1_ddpd, prim_buffer_1_sdpd, prim_buffer_2_sdpd, prim_buffer_1_pppd, prim_buffer_2_pppd, prim_buffer_2_pdsd, prim_buffer_2_pdpp, prim_buffer_1_pdpd, prim_buffer_2_pdpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dddp(prim_buffer_1_dddp, prim_buffer_1_sddp, prim_buffer_2_sddp, prim_buffer_1_ppdp, prim_buffer_2_ppdp, prim_buffer_2_pdpp, prim_buffer_2_pdds, prim_buffer_1_pddp, prim_buffer_2_pddp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dddd(prim_buffer_0_dddd, prim_buffer_0_sddd, prim_buffer_1_sddd, prim_buffer_0_ppdd, prim_buffer_1_ppdd, prim_buffer_1_pdpd, prim_buffer_1_pddp, prim_buffer_0_pddd, prim_buffer_1_pddd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dddd(prim_buffer_1_dddd, prim_buffer_1_sddd, prim_buffer_2_sddd, prim_buffer_1_ppdd, prim_buffer_2_ppdd, prim_buffer_2_pdpd, prim_buffer_2_pddp, prim_buffer_1_pddd, prim_buffer_2_pddd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_fddd(prim_buffer_0_fddd, prim_buffer_0_pddd, prim_buffer_1_pddd, prim_buffer_0_dpdd, prim_buffer_1_dpdd, prim_buffer_1_ddpd, prim_buffer_1_dddp, prim_buffer_0_dddd, prim_buffer_1_dddd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            t4c_geom::comp_geom2000_pddd_0(prim_buffer_0_geom2000_pddd, prim_buffer_0_pddd, prim_buffer_0_fddd, a_exp );

            t2cfunc::reduce(cart_buffer_0_geom2000_pddd, prim_buffer_0_geom2000_pddd, ket_dim, ket_npgtos);
        }

        t4cfunc::full_transform<1, 2, 2, 2>(spher_buffer_0_geom2000_pddd, cart_spher_buffer_0_geom2000_pddd);

        distributor->distribute(spher_buffer_0_geom2000_pddd, a_indices, b_indices, c_indices, d_indices, 1, 2, 2, 2, 2, 0, 0, 0, i, ket_indices);
    }
}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecPDDD_hpp */