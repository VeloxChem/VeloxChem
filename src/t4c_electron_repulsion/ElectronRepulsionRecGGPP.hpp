#ifndef ElectronRepulsionRecGGPP_hpp
#define ElectronRepulsionRecGGPP_hpp

#include <array>

#include "ElectronRepulsionContrRecDGXX.hpp"
#include "ElectronRepulsionContrRecDHXX.hpp"
#include "ElectronRepulsionContrRecDIXX.hpp"
#include "ElectronRepulsionContrRecFGXX.hpp"
#include "ElectronRepulsionContrRecFHXX.hpp"
#include "ElectronRepulsionContrRecGGXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecPIXX.hpp"
#include "ElectronRepulsionContrRecPKXX.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
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
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSP.hpp"
#include "ElectronRepulsionPrimRecSKSS.hpp"
#include "ElectronRepulsionPrimRecSLSD.hpp"
#include "ElectronRepulsionPrimRecSLSP.hpp"
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
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// Computes (GG|1/|r-r'||PP)  integrals for two GTOs pair blocks.
/// - Parameter distributor: the pointer to Fock matrix/matrices distributor.
/// - Parameter bra_gto_pair_block: the GTOs pair block on bra side.
/// - Parameter ket_gto_pair_block: the GTOs pair block on ket side.
/// - Parameter bra_indices: the range [bra_first, bra_last) of GTOs on bra side.
/// - Parameter ket_indices: the range [ket_first, ket_last) of GTOs on ket side.
template <class T>
auto
comp_electron_repulsion_ggpp(T* distributor,
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

    CSimdArray<double> prim_buffer_9_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_10_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_8_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_9_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_8_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_1_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_2_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_3_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_4_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_5_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_6_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_7_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_0_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_6_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_7_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_0_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_6_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_7_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdss(6, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdss(6, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdss(6, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdss(6, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdss(6, ket_pdim);

    CSimdArray<double> prim_buffer_6_sdss(6, ket_pdim);

    CSimdArray<double> prim_buffer_0_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_6_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_0_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_6_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_1_sfss(10, ket_pdim);

    CSimdArray<double> prim_buffer_2_sfss(10, ket_pdim);

    CSimdArray<double> prim_buffer_3_sfss(10, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfss(10, ket_pdim);

    CSimdArray<double> prim_buffer_5_sfss(10, ket_pdim);

    CSimdArray<double> prim_buffer_0_sfsp(30, ket_pdim);

    CSimdArray<double> prim_buffer_1_sfsp(30, ket_pdim);

    CSimdArray<double> prim_buffer_2_sfsp(30, ket_pdim);

    CSimdArray<double> prim_buffer_3_sfsp(30, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfsp(30, ket_pdim);

    CSimdArray<double> prim_buffer_5_sfsp(30, ket_pdim);

    CSimdArray<double> prim_buffer_0_sfsd(60, ket_pdim);

    CSimdArray<double> prim_buffer_1_sfsd(60, ket_pdim);

    CSimdArray<double> prim_buffer_2_sfsd(60, ket_pdim);

    CSimdArray<double> prim_buffer_3_sfsd(60, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfsd(60, ket_pdim);

    CSimdArray<double> prim_buffer_5_sfsd(60, ket_pdim);

    CSimdArray<double> prim_buffer_1_sgss(15, ket_pdim);

    CSimdArray<double> prim_buffer_2_sgss(15, ket_pdim);

    CSimdArray<double> prim_buffer_3_sgss(15, ket_pdim);

    CSimdArray<double> prim_buffer_4_sgss(15, ket_pdim);

    CSimdArray<double> prim_buffer_0_sgsp(45, ket_pdim);

    CSimdArray<double> prim_buffer_1_sgsp(45, ket_pdim);

    CSimdArray<double> prim_buffer_2_sgsp(45, ket_pdim);

    CSimdArray<double> prim_buffer_3_sgsp(45, ket_pdim);

    CSimdArray<double> prim_buffer_4_sgsp(45, ket_pdim);

    CSimdArray<double> prim_buffer_0_sgsd(90, ket_pdim);

    CSimdArray<double> prim_buffer_1_sgsd(90, ket_pdim);

    CSimdArray<double> prim_buffer_2_sgsd(90, ket_pdim);

    CSimdArray<double> prim_buffer_3_sgsd(90, ket_pdim);

    CSimdArray<double> prim_buffer_4_sgsd(90, ket_pdim);

    CSimdArray<double> prim_buffer_1_shss(21, ket_pdim);

    CSimdArray<double> prim_buffer_2_shss(21, ket_pdim);

    CSimdArray<double> prim_buffer_3_shss(21, ket_pdim);

    CSimdArray<double> prim_buffer_0_shsp(63, ket_pdim);

    CSimdArray<double> prim_buffer_1_shsp(63, ket_pdim);

    CSimdArray<double> prim_buffer_2_shsp(63, ket_pdim);

    CSimdArray<double> prim_buffer_3_shsp(63, ket_pdim);

    CSimdArray<double> prim_buffer_0_shsd(126, ket_pdim);

    CSimdArray<double> prim_buffer_1_shsd(126, ket_pdim);

    CSimdArray<double> prim_buffer_2_shsd(126, ket_pdim);

    CSimdArray<double> prim_buffer_3_shsd(126, ket_pdim);

    CSimdArray<double> prim_buffer_1_siss(28, ket_pdim);

    CSimdArray<double> prim_buffer_2_siss(28, ket_pdim);

    CSimdArray<double> prim_buffer_0_sisp(84, ket_pdim);

    CSimdArray<double> prim_buffer_1_sisp(84, ket_pdim);

    CSimdArray<double> prim_buffer_2_sisp(84, ket_pdim);

    CSimdArray<double> prim_buffer_0_sisd(168, ket_pdim);

    CSimdArray<double> prim_buffer_1_sisd(168, ket_pdim);

    CSimdArray<double> prim_buffer_2_sisd(168, ket_pdim);

    CSimdArray<double> prim_buffer_1_skss(36, ket_pdim);

    CSimdArray<double> prim_buffer_0_sksp(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_sksp(108, ket_pdim);

    CSimdArray<double> prim_buffer_0_sksd(216, ket_pdim);

    CSimdArray<double> prim_buffer_1_sksd(216, ket_pdim);

    CSimdArray<double> prim_buffer_0_slsp(135, ket_pdim);

    CSimdArray<double> prim_buffer_0_slsd(270, ket_pdim);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cart_buffer_0_sgsp(45, ket_dim);

    CSimdArray<double> cart_buffer_0_sgsd(90, ket_dim);

    CSimdArray<double> cart_buffer_0_shsp(63, ket_dim);

    CSimdArray<double> cart_buffer_0_shsd(126, ket_dim);

    CSimdArray<double> cart_buffer_0_sisp(84, ket_dim);

    CSimdArray<double> cart_buffer_0_sisd(168, ket_dim);

    CSimdArray<double> cart_buffer_0_sksp(108, ket_dim);

    CSimdArray<double> cart_buffer_0_sksd(216, ket_dim);

    CSimdArray<double> cart_buffer_0_slsp(135, ket_dim);

    CSimdArray<double> cart_buffer_0_slsd(270, ket_dim);

    // allocate aligned contracted integrals

    CSimdArray<double> contr_buffer_0_sgpp(135, ket_dim);

    CSimdArray<double> contr_buffer_0_shpp(189, ket_dim);

    CSimdArray<double> contr_buffer_0_sipp(252, ket_dim);

    CSimdArray<double> contr_buffer_0_skpp(324, ket_dim);

    CSimdArray<double> contr_buffer_0_slpp(405, ket_dim);

    // allocate aligned half transformed integrals

    CSimdArray<double> ket_spher_buffer_0_sgpp(135, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_shpp(189, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_sipp(252, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_skpp(324, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_slpp(405, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_pgpp(405, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_phpp(567, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_pipp(756, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_pkpp(972, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_dgpp(810, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_dhpp(1134, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_dipp(1512, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_fgpp(1350, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_fhpp(1890, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_ggpp(2025, ket_dim);

    // allocate aligned spherical integrals

    CSimdArray<double> spher_buffer_0_ggpp(729, ket_dim);

    // allocate accumulation buffer for integrals

    CSimdArray<double> buffer(bra_dim * 729, ket_dim);

    // setup Boys fuction data

    const CBoysFunc<10> bf_table;

    CSimdArray<double> bf_args(1, ket_pdim);

    CSimdArray<double> bf_values(11, ket_pdim);

    // loop over contracted GTOs on bra side

    for (auto i = bra_indices[0]; i < bra_indices[1]; i++)
    {
        // zero integral buffers

        cart_buffer_0_sgsp.zero();

        cart_buffer_0_sgsd.zero();

        cart_buffer_0_shsp.zero();

        cart_buffer_0_shsd.zero();

        cart_buffer_0_sisp.zero();

        cart_buffer_0_sisd.zero();

        cart_buffer_0_sksp.zero();

        cart_buffer_0_sksd.zero();

        cart_buffer_0_slsp.zero();

        cart_buffer_0_slsd.zero();

        ket_spher_buffer_0_sgpp.zero();

        ket_spher_buffer_0_shpp.zero();

        ket_spher_buffer_0_sipp.zero();

        ket_spher_buffer_0_skpp.zero();

        ket_spher_buffer_0_slpp.zero();

        ket_spher_buffer_0_pgpp.zero();

        ket_spher_buffer_0_phpp.zero();

        ket_spher_buffer_0_pipp.zero();

        ket_spher_buffer_0_pkpp.zero();

        ket_spher_buffer_0_dgpp.zero();

        ket_spher_buffer_0_dhpp.zero();

        ket_spher_buffer_0_dipp.zero();

        ket_spher_buffer_0_fgpp.zero();

        ket_spher_buffer_0_fhpp.zero();

        ket_spher_buffer_0_ggpp.zero();

        spher_buffer_0_ggpp.zero();

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

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_9_ssss, fss_abcd[0], bf_values[9]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_10_ssss, fss_abcd[0], bf_values[10]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_0_sssp, prim_buffer_0_ssss, prim_buffer_1_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_1_sssp, prim_buffer_1_ssss, prim_buffer_2_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_2_sssp, prim_buffer_2_ssss, prim_buffer_3_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_3_sssp, prim_buffer_3_ssss, prim_buffer_4_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_4_sssp, prim_buffer_4_ssss, prim_buffer_5_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_5_sssp, prim_buffer_5_ssss, prim_buffer_6_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_6_sssp, prim_buffer_6_ssss, prim_buffer_7_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_7_sssp, prim_buffer_7_ssss, prim_buffer_8_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_8_sssp, prim_buffer_8_ssss, prim_buffer_9_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_9_sssp, prim_buffer_9_ssss, prim_buffer_10_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_0_sssd, prim_buffer_0_ssss, prim_buffer_1_ssss, prim_buffer_0_sssp, prim_buffer_1_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_1_sssd, prim_buffer_1_ssss, prim_buffer_2_ssss, prim_buffer_1_sssp, prim_buffer_2_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_2_sssd, prim_buffer_2_ssss, prim_buffer_3_ssss, prim_buffer_2_sssp, prim_buffer_3_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_3_sssd, prim_buffer_3_ssss, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_4_sssd, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_5_sssd, prim_buffer_5_ssss, prim_buffer_6_ssss, prim_buffer_5_sssp, prim_buffer_6_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_6_sssd, prim_buffer_6_ssss, prim_buffer_7_ssss, prim_buffer_6_sssp, prim_buffer_7_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_7_sssd, prim_buffer_7_ssss, prim_buffer_8_ssss, prim_buffer_7_sssp, prim_buffer_8_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_8_sssd, prim_buffer_8_ssss, prim_buffer_9_ssss, prim_buffer_8_sssp, prim_buffer_9_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_1_spss, prim_buffer_1_ssss, prim_buffer_2_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_2_spss, prim_buffer_2_ssss, prim_buffer_3_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_3_spss, prim_buffer_3_ssss, prim_buffer_4_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_4_spss, prim_buffer_4_ssss, prim_buffer_5_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_5_spss, prim_buffer_5_ssss, prim_buffer_6_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_6_spss, prim_buffer_6_ssss, prim_buffer_7_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_7_spss, prim_buffer_7_ssss, prim_buffer_8_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_0_spsp, prim_buffer_1_ssss, prim_buffer_0_sssp, prim_buffer_1_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_1_spsp, prim_buffer_2_ssss, prim_buffer_1_sssp, prim_buffer_2_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_2_spsp, prim_buffer_3_ssss, prim_buffer_2_sssp, prim_buffer_3_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_3_spsp, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_4_spsp, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_5_spsp, prim_buffer_6_ssss, prim_buffer_5_sssp, prim_buffer_6_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_6_spsp, prim_buffer_7_ssss, prim_buffer_6_sssp, prim_buffer_7_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_7_spsp, prim_buffer_8_ssss, prim_buffer_7_sssp, prim_buffer_8_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_0_spsd, prim_buffer_1_sssp, prim_buffer_0_sssd, prim_buffer_1_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_1_spsd, prim_buffer_2_sssp, prim_buffer_1_sssd, prim_buffer_2_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_2_spsd, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_3_spsd, prim_buffer_4_sssp, prim_buffer_3_sssd, prim_buffer_4_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_4_spsd, prim_buffer_5_sssp, prim_buffer_4_sssd, prim_buffer_5_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_5_spsd, prim_buffer_6_sssp, prim_buffer_5_sssd, prim_buffer_6_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_6_spsd, prim_buffer_7_sssp, prim_buffer_6_sssd, prim_buffer_7_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_7_spsd, prim_buffer_8_sssp, prim_buffer_7_sssd, prim_buffer_8_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_1_sdss, prim_buffer_1_ssss, prim_buffer_2_ssss, prim_buffer_1_spss, prim_buffer_2_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_2_sdss, prim_buffer_2_ssss, prim_buffer_3_ssss, prim_buffer_2_spss, prim_buffer_3_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_3_sdss, prim_buffer_3_ssss, prim_buffer_4_ssss, prim_buffer_3_spss, prim_buffer_4_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_4_sdss, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_spss, prim_buffer_5_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_5_sdss, prim_buffer_5_ssss, prim_buffer_6_ssss, prim_buffer_5_spss, prim_buffer_6_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_6_sdss, prim_buffer_6_ssss, prim_buffer_7_ssss, prim_buffer_6_spss, prim_buffer_7_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_0_sdsp, prim_buffer_0_sssp, prim_buffer_1_sssp, prim_buffer_1_spss, prim_buffer_0_spsp, prim_buffer_1_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_1_sdsp, prim_buffer_1_sssp, prim_buffer_2_sssp, prim_buffer_2_spss, prim_buffer_1_spsp, prim_buffer_2_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_2_sdsp, prim_buffer_2_sssp, prim_buffer_3_sssp, prim_buffer_3_spss, prim_buffer_2_spsp, prim_buffer_3_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_3_sdsp, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_4_spss, prim_buffer_3_spsp, prim_buffer_4_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_4_sdsp, prim_buffer_4_sssp, prim_buffer_5_sssp, prim_buffer_5_spss, prim_buffer_4_spsp, prim_buffer_5_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_5_sdsp, prim_buffer_5_sssp, prim_buffer_6_sssp, prim_buffer_6_spss, prim_buffer_5_spsp, prim_buffer_6_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_6_sdsp, prim_buffer_6_sssp, prim_buffer_7_sssp, prim_buffer_7_spss, prim_buffer_6_spsp, prim_buffer_7_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_0_sdsd, prim_buffer_0_sssd, prim_buffer_1_sssd, prim_buffer_1_spsp, prim_buffer_0_spsd, prim_buffer_1_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_1_sdsd, prim_buffer_1_sssd, prim_buffer_2_sssd, prim_buffer_2_spsp, prim_buffer_1_spsd, prim_buffer_2_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_2_sdsd, prim_buffer_2_sssd, prim_buffer_3_sssd, prim_buffer_3_spsp, prim_buffer_2_spsd, prim_buffer_3_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_3_sdsd, prim_buffer_3_sssd, prim_buffer_4_sssd, prim_buffer_4_spsp, prim_buffer_3_spsd, prim_buffer_4_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_4_sdsd, prim_buffer_4_sssd, prim_buffer_5_sssd, prim_buffer_5_spsp, prim_buffer_4_spsd, prim_buffer_5_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_5_sdsd, prim_buffer_5_sssd, prim_buffer_6_sssd, prim_buffer_6_spsp, prim_buffer_5_spsd, prim_buffer_6_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_6_sdsd, prim_buffer_6_sssd, prim_buffer_7_sssd, prim_buffer_7_spsp, prim_buffer_6_spsd, prim_buffer_7_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfss(prim_buffer_1_sfss, prim_buffer_1_spss, prim_buffer_2_spss, prim_buffer_1_sdss, prim_buffer_2_sdss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfss(prim_buffer_2_sfss, prim_buffer_2_spss, prim_buffer_3_spss, prim_buffer_2_sdss, prim_buffer_3_sdss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfss(prim_buffer_3_sfss, prim_buffer_3_spss, prim_buffer_4_spss, prim_buffer_3_sdss, prim_buffer_4_sdss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfss(prim_buffer_4_sfss, prim_buffer_4_spss, prim_buffer_5_spss, prim_buffer_4_sdss, prim_buffer_5_sdss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfss(prim_buffer_5_sfss, prim_buffer_5_spss, prim_buffer_6_spss, prim_buffer_5_sdss, prim_buffer_6_sdss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_0_sfsp, prim_buffer_0_spsp, prim_buffer_1_spsp, prim_buffer_1_sdss, prim_buffer_0_sdsp, prim_buffer_1_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_1_sfsp, prim_buffer_1_spsp, prim_buffer_2_spsp, prim_buffer_2_sdss, prim_buffer_1_sdsp, prim_buffer_2_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_2_sfsp, prim_buffer_2_spsp, prim_buffer_3_spsp, prim_buffer_3_sdss, prim_buffer_2_sdsp, prim_buffer_3_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_3_sfsp, prim_buffer_3_spsp, prim_buffer_4_spsp, prim_buffer_4_sdss, prim_buffer_3_sdsp, prim_buffer_4_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_4_sfsp, prim_buffer_4_spsp, prim_buffer_5_spsp, prim_buffer_5_sdss, prim_buffer_4_sdsp, prim_buffer_5_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_5_sfsp, prim_buffer_5_spsp, prim_buffer_6_spsp, prim_buffer_6_sdss, prim_buffer_5_sdsp, prim_buffer_6_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_0_sfsd, prim_buffer_0_spsd, prim_buffer_1_spsd, prim_buffer_1_sdsp, prim_buffer_0_sdsd, prim_buffer_1_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_1_sfsd, prim_buffer_1_spsd, prim_buffer_2_spsd, prim_buffer_2_sdsp, prim_buffer_1_sdsd, prim_buffer_2_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_2_sfsd, prim_buffer_2_spsd, prim_buffer_3_spsd, prim_buffer_3_sdsp, prim_buffer_2_sdsd, prim_buffer_3_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_3_sfsd, prim_buffer_3_spsd, prim_buffer_4_spsd, prim_buffer_4_sdsp, prim_buffer_3_sdsd, prim_buffer_4_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_4_sfsd, prim_buffer_4_spsd, prim_buffer_5_spsd, prim_buffer_5_sdsp, prim_buffer_4_sdsd, prim_buffer_5_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_5_sfsd, prim_buffer_5_spsd, prim_buffer_6_spsd, prim_buffer_6_sdsp, prim_buffer_5_sdsd, prim_buffer_6_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgss(prim_buffer_1_sgss, prim_buffer_1_sdss, prim_buffer_2_sdss, prim_buffer_1_sfss, prim_buffer_2_sfss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgss(prim_buffer_2_sgss, prim_buffer_2_sdss, prim_buffer_3_sdss, prim_buffer_2_sfss, prim_buffer_3_sfss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgss(prim_buffer_3_sgss, prim_buffer_3_sdss, prim_buffer_4_sdss, prim_buffer_3_sfss, prim_buffer_4_sfss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgss(prim_buffer_4_sgss, prim_buffer_4_sdss, prim_buffer_5_sdss, prim_buffer_4_sfss, prim_buffer_5_sfss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsp(prim_buffer_0_sgsp, prim_buffer_0_sdsp, prim_buffer_1_sdsp, prim_buffer_1_sfss, prim_buffer_0_sfsp, prim_buffer_1_sfsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsp(prim_buffer_1_sgsp, prim_buffer_1_sdsp, prim_buffer_2_sdsp, prim_buffer_2_sfss, prim_buffer_1_sfsp, prim_buffer_2_sfsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsp(prim_buffer_2_sgsp, prim_buffer_2_sdsp, prim_buffer_3_sdsp, prim_buffer_3_sfss, prim_buffer_2_sfsp, prim_buffer_3_sfsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsp(prim_buffer_3_sgsp, prim_buffer_3_sdsp, prim_buffer_4_sdsp, prim_buffer_4_sfss, prim_buffer_3_sfsp, prim_buffer_4_sfsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsp(prim_buffer_4_sgsp, prim_buffer_4_sdsp, prim_buffer_5_sdsp, prim_buffer_5_sfss, prim_buffer_4_sfsp, prim_buffer_5_sfsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsd(prim_buffer_0_sgsd, prim_buffer_0_sdsd, prim_buffer_1_sdsd, prim_buffer_1_sfsp, prim_buffer_0_sfsd, prim_buffer_1_sfsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsd(prim_buffer_1_sgsd, prim_buffer_1_sdsd, prim_buffer_2_sdsd, prim_buffer_2_sfsp, prim_buffer_1_sfsd, prim_buffer_2_sfsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsd(prim_buffer_2_sgsd, prim_buffer_2_sdsd, prim_buffer_3_sdsd, prim_buffer_3_sfsp, prim_buffer_2_sfsd, prim_buffer_3_sfsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsd(prim_buffer_3_sgsd, prim_buffer_3_sdsd, prim_buffer_4_sdsd, prim_buffer_4_sfsp, prim_buffer_3_sfsd, prim_buffer_4_sfsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsd(prim_buffer_4_sgsd, prim_buffer_4_sdsd, prim_buffer_5_sdsd, prim_buffer_5_sfsp, prim_buffer_4_sfsd, prim_buffer_5_sfsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shss(prim_buffer_1_shss, prim_buffer_1_sfss, prim_buffer_2_sfss, prim_buffer_1_sgss, prim_buffer_2_sgss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shss(prim_buffer_2_shss, prim_buffer_2_sfss, prim_buffer_3_sfss, prim_buffer_2_sgss, prim_buffer_3_sgss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shss(prim_buffer_3_shss, prim_buffer_3_sfss, prim_buffer_4_sfss, prim_buffer_3_sgss, prim_buffer_4_sgss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsp(prim_buffer_0_shsp, prim_buffer_0_sfsp, prim_buffer_1_sfsp, prim_buffer_1_sgss, prim_buffer_0_sgsp, prim_buffer_1_sgsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsp(prim_buffer_1_shsp, prim_buffer_1_sfsp, prim_buffer_2_sfsp, prim_buffer_2_sgss, prim_buffer_1_sgsp, prim_buffer_2_sgsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsp(prim_buffer_2_shsp, prim_buffer_2_sfsp, prim_buffer_3_sfsp, prim_buffer_3_sgss, prim_buffer_2_sgsp, prim_buffer_3_sgsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsp(prim_buffer_3_shsp, prim_buffer_3_sfsp, prim_buffer_4_sfsp, prim_buffer_4_sgss, prim_buffer_3_sgsp, prim_buffer_4_sgsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsd(prim_buffer_0_shsd, prim_buffer_0_sfsd, prim_buffer_1_sfsd, prim_buffer_1_sgsp, prim_buffer_0_sgsd, prim_buffer_1_sgsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsd(prim_buffer_1_shsd, prim_buffer_1_sfsd, prim_buffer_2_sfsd, prim_buffer_2_sgsp, prim_buffer_1_sgsd, prim_buffer_2_sgsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsd(prim_buffer_2_shsd, prim_buffer_2_sfsd, prim_buffer_3_sfsd, prim_buffer_3_sgsp, prim_buffer_2_sgsd, prim_buffer_3_sgsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsd(prim_buffer_3_shsd, prim_buffer_3_sfsd, prim_buffer_4_sfsd, prim_buffer_4_sgsp, prim_buffer_3_sgsd, prim_buffer_4_sgsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_siss(prim_buffer_1_siss, prim_buffer_1_sgss, prim_buffer_2_sgss, prim_buffer_1_shss, prim_buffer_2_shss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_siss(prim_buffer_2_siss, prim_buffer_2_sgss, prim_buffer_3_sgss, prim_buffer_2_shss, prim_buffer_3_shss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisp(prim_buffer_0_sisp, prim_buffer_0_sgsp, prim_buffer_1_sgsp, prim_buffer_1_shss, prim_buffer_0_shsp, prim_buffer_1_shsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisp(prim_buffer_1_sisp, prim_buffer_1_sgsp, prim_buffer_2_sgsp, prim_buffer_2_shss, prim_buffer_1_shsp, prim_buffer_2_shsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisp(prim_buffer_2_sisp, prim_buffer_2_sgsp, prim_buffer_3_sgsp, prim_buffer_3_shss, prim_buffer_2_shsp, prim_buffer_3_shsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisd(prim_buffer_0_sisd, prim_buffer_0_sgsd, prim_buffer_1_sgsd, prim_buffer_1_shsp, prim_buffer_0_shsd, prim_buffer_1_shsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisd(prim_buffer_1_sisd, prim_buffer_1_sgsd, prim_buffer_2_sgsd, prim_buffer_2_shsp, prim_buffer_1_shsd, prim_buffer_2_shsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisd(prim_buffer_2_sisd, prim_buffer_2_sgsd, prim_buffer_3_sgsd, prim_buffer_3_shsp, prim_buffer_2_shsd, prim_buffer_3_shsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_skss(prim_buffer_1_skss, prim_buffer_1_shss, prim_buffer_2_shss, prim_buffer_1_siss, prim_buffer_2_siss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksp(prim_buffer_0_sksp, prim_buffer_0_shsp, prim_buffer_1_shsp, prim_buffer_1_siss, prim_buffer_0_sisp, prim_buffer_1_sisp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksp(prim_buffer_1_sksp, prim_buffer_1_shsp, prim_buffer_2_shsp, prim_buffer_2_siss, prim_buffer_1_sisp, prim_buffer_2_sisp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksd(prim_buffer_0_sksd, prim_buffer_0_shsd, prim_buffer_1_shsd, prim_buffer_1_sisp, prim_buffer_0_sisd, prim_buffer_1_sisd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksd(prim_buffer_1_sksd, prim_buffer_1_shsd, prim_buffer_2_shsd, prim_buffer_2_sisp, prim_buffer_1_sisd, prim_buffer_2_sisd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_slsp(prim_buffer_0_slsp, prim_buffer_0_sisp, prim_buffer_1_sisp, prim_buffer_1_skss, prim_buffer_0_sksp, prim_buffer_1_sksp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_slsd(prim_buffer_0_slsd, prim_buffer_0_sisd, prim_buffer_1_sisd, prim_buffer_1_sksp, prim_buffer_0_sksd, prim_buffer_1_sksd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            t2cfunc::reduce(cart_buffer_0_sgsp, prim_buffer_0_sgsp, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sgsd, prim_buffer_0_sgsd, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_shsp, prim_buffer_0_shsp, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_shsd, prim_buffer_0_shsd, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sisp, prim_buffer_0_sisp, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sisd, prim_buffer_0_sisd, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sksp, prim_buffer_0_sksp, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sksd, prim_buffer_0_sksd, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_slsp, prim_buffer_0_slsp, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_slsd, prim_buffer_0_slsd, ket_dim, ket_npgtos);

        }

        erirec::comp_ket_hrr_electron_repulsion_xxpp(contr_buffer_0_sgpp, cart_buffer_0_sgsp, cart_buffer_0_sgsd, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxpp(contr_buffer_0_shpp, cart_buffer_0_shsp, cart_buffer_0_shsd, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxpp(contr_buffer_0_sipp, cart_buffer_0_sisp, cart_buffer_0_sisd, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxpp(contr_buffer_0_skpp, cart_buffer_0_sksp, cart_buffer_0_sksd, cd_x[0], cd_y[0], cd_z[0], 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxpp(contr_buffer_0_slpp, cart_buffer_0_slsp, cart_buffer_0_slsd, cd_x[0], cd_y[0], cd_z[0], 0, 8);

        t4cfunc::ket_transform<1, 1>(ket_spher_buffer_0_sgpp, contr_buffer_0_sgpp, 0, 4);

        t4cfunc::ket_transform<1, 1>(ket_spher_buffer_0_shpp, contr_buffer_0_shpp, 0, 5);

        t4cfunc::ket_transform<1, 1>(ket_spher_buffer_0_sipp, contr_buffer_0_sipp, 0, 6);

        t4cfunc::ket_transform<1, 1>(ket_spher_buffer_0_skpp, contr_buffer_0_skpp, 0, 7);

        t4cfunc::ket_transform<1, 1>(ket_spher_buffer_0_slpp, contr_buffer_0_slpp, 0, 8);

        erirec::comp_bra_hrr_electron_repulsion_pgxx(ket_spher_buffer_0_pgpp, ket_spher_buffer_0_sgpp, ket_spher_buffer_0_shpp, ab_x, ab_y, ab_z, 1, 1);

        erirec::comp_bra_hrr_electron_repulsion_phxx(ket_spher_buffer_0_phpp, ket_spher_buffer_0_shpp, ket_spher_buffer_0_sipp, ab_x, ab_y, ab_z, 1, 1);

        erirec::comp_bra_hrr_electron_repulsion_pixx(ket_spher_buffer_0_pipp, ket_spher_buffer_0_sipp, ket_spher_buffer_0_skpp, ab_x, ab_y, ab_z, 1, 1);

        erirec::comp_bra_hrr_electron_repulsion_pkxx(ket_spher_buffer_0_pkpp, ket_spher_buffer_0_skpp, ket_spher_buffer_0_slpp, ab_x, ab_y, ab_z, 1, 1);

        erirec::comp_bra_hrr_electron_repulsion_dgxx(ket_spher_buffer_0_dgpp, ket_spher_buffer_0_pgpp, ket_spher_buffer_0_phpp, ab_x, ab_y, ab_z, 1, 1);

        erirec::comp_bra_hrr_electron_repulsion_dhxx(ket_spher_buffer_0_dhpp, ket_spher_buffer_0_phpp, ket_spher_buffer_0_pipp, ab_x, ab_y, ab_z, 1, 1);

        erirec::comp_bra_hrr_electron_repulsion_dixx(ket_spher_buffer_0_dipp, ket_spher_buffer_0_pipp, ket_spher_buffer_0_pkpp, ab_x, ab_y, ab_z, 1, 1);

        erirec::comp_bra_hrr_electron_repulsion_fgxx(ket_spher_buffer_0_fgpp, ket_spher_buffer_0_dgpp, ket_spher_buffer_0_dhpp, ab_x, ab_y, ab_z, 1, 1);

        erirec::comp_bra_hrr_electron_repulsion_fhxx(ket_spher_buffer_0_fhpp, ket_spher_buffer_0_dhpp, ket_spher_buffer_0_dipp, ab_x, ab_y, ab_z, 1, 1);

        erirec::comp_bra_hrr_electron_repulsion_ggxx(ket_spher_buffer_0_ggpp, ket_spher_buffer_0_fgpp, ket_spher_buffer_0_fhpp, ab_x, ab_y, ab_z, 1, 1);

        t4cfunc::bra_transform<4, 4>(spher_buffer_0_ggpp, ket_spher_buffer_0_ggpp, 1, 1);

        t4cfunc::store_values(buffer, spher_buffer_0_ggpp, 729 * (i - bra_indices[0]));
    }

    distributor->distribute(buffer, a_indices, b_indices, c_indices, d_indices, 4, 4, 1, 1, bra_indices, ket_indices);
}

} // erirec namespace

#endif /* ElectronRepulsionRecGGPP_hpp */
