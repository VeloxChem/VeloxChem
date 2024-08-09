#ifndef ElectronRepulsionDiagRecDGDG_hpp
#define ElectronRepulsionDiagRecDGDG_hpp

#include <vector>
#include <array>

#include "ElectronRepulsionContrRecDGXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSI.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSI.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSI.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSI.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
#include "ElectronRepulsionPrimRecSISI.hpp"
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
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// Computes (DG|1/|r-r'||DG)  integrals for GTOs pair block.
/// - Parameter distributor: the pointer to screening data distributor.
/// - Parameter gto_pair_block: the GTOs pair block.
/// - Parameter go_indices: the range [gto_first, gto_last) of GTOs on bra and ket sides.
template <class T>
auto
comp_diag_electron_repulsion_dgdg(T* distributor,
                                  const CGtoPairBlock& gto_pair_block,
                                  const std::array<int, 2>& gto_indices) -> void
{
    // intialize GTOs pair data

    const auto a_coords_x = gto_pair_block.bra_coordinates_x();

    const auto a_coords_y = gto_pair_block.bra_coordinates_y();

    const auto a_coords_z = gto_pair_block.bra_coordinates_z();

    const auto b_coords_x = gto_pair_block.ket_coordinates_x();

    const auto b_coords_y = gto_pair_block.ket_coordinates_y();

    const auto b_coords_z = gto_pair_block.ket_coordinates_z();

    const auto a_vec_exps = gto_pair_block.bra_exponents();

    const auto b_vec_exps = gto_pair_block.ket_exponents();

    const auto ab_vec_norms = gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = gto_pair_block.overlap_factors();

    const auto a_indices = gto_pair_block.bra_orbital_indices();

    const auto b_indices = gto_pair_block.ket_orbital_indices();

    const auto ncgtos = gto_pair_block.number_of_contracted_pairs();

    const auto npgtos = gto_pair_block.number_of_primitive_pairs();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> c_x(1, npgtos);

    CSimdArray<double> c_y(1, npgtos);

    CSimdArray<double> c_z(1, npgtos);

    CSimdArray<double> d_x(1, npgtos);

    CSimdArray<double> d_y(1, npgtos);

    CSimdArray<double> d_z(1, npgtos);

    CSimdArray<double> c_exps(1, npgtos);

    CSimdArray<double> d_exps(1, npgtos);

    CSimdArray<double> cd_norms(1, npgtos);

    CSimdArray<double> cd_ovls(1, npgtos);

    // allocate aligned coordinates of Q center

    CSimdArray<double> q_x(1, npgtos);

    CSimdArray<double> q_y(1, npgtos);

    CSimdArray<double> q_z(1, npgtos);

    // allocate aligned coordinates of W center

    CSimdArray<double> w_x(1, npgtos);

    CSimdArray<double> w_y(1, npgtos);

    CSimdArray<double> w_z(1, npgtos);

    // allocate aligned distances R(PQ) = P - Q

    CSimdArray<double> pq_x(1, npgtos);

    CSimdArray<double> pq_y(1, npgtos);

    CSimdArray<double> pq_z(1, npgtos);

    // allocate aligned distances R(QD) = Q - D

    CSimdArray<double> qd_x(1, npgtos);

    CSimdArray<double> qd_y(1, npgtos);

    CSimdArray<double> qd_z(1, npgtos);

    // allocate aligned distances R(WQ) = W - Q

    CSimdArray<double> wq_x(1, npgtos);

    CSimdArray<double> wq_y(1, npgtos);

    CSimdArray<double> wq_z(1, npgtos);

    // allocate aligned distances R(WP) = W - P

    CSimdArray<double> wp_x(1, npgtos);

    CSimdArray<double> wp_y(1, npgtos);

    CSimdArray<double> wp_z(1, npgtos);

    // allocate combined overlap factor

    CSimdArray<double> fss_abcd(1, npgtos);

    // allocate and initialize aligned distances R(CD) = C - D

    CSimdArray<double> cd_x(1, 1);

    CSimdArray<double> cd_y(1, 1);

    CSimdArray<double> cd_z(1, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> prim_buffer_0_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_1_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_2_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_3_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_4_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_5_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_6_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_7_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_8_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_9_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_10_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_11_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_12_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_0_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_1_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_2_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_3_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_4_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_5_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_6_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_7_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_8_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_9_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_10_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_11_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_0_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_1_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_2_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_3_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_4_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_5_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_6_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_7_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_8_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_9_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_10_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_0_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_1_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_2_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_3_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_4_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_5_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_6_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_7_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_8_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_9_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_0_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_1_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_2_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_3_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_4_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_5_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_6_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_7_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_8_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_0_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_1_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_2_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_3_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_4_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_5_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_6_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_7_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_0_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_1_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_2_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_3_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_4_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_5_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_6_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_4_spss(3, npgtos);

    CSimdArray<double> prim_buffer_5_spss(3, npgtos);

    CSimdArray<double> prim_buffer_3_spsp(9, npgtos);

    CSimdArray<double> prim_buffer_4_spsp(9, npgtos);

    CSimdArray<double> prim_buffer_5_spsp(9, npgtos);

    CSimdArray<double> prim_buffer_2_spsd(18, npgtos);

    CSimdArray<double> prim_buffer_3_spsd(18, npgtos);

    CSimdArray<double> prim_buffer_4_spsd(18, npgtos);

    CSimdArray<double> prim_buffer_5_spsd(18, npgtos);

    CSimdArray<double> prim_buffer_1_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_2_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_3_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_4_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_5_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_0_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_1_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_2_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_3_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_4_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_5_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_0_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_1_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_2_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_3_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_4_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_5_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_0_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_1_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_2_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_3_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_4_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_5_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_4_sdss(6, npgtos);

    CSimdArray<double> prim_buffer_3_sdsp(18, npgtos);

    CSimdArray<double> prim_buffer_4_sdsp(18, npgtos);

    CSimdArray<double> prim_buffer_2_sdsd(36, npgtos);

    CSimdArray<double> prim_buffer_3_sdsd(36, npgtos);

    CSimdArray<double> prim_buffer_4_sdsd(36, npgtos);

    CSimdArray<double> prim_buffer_1_sdsf(60, npgtos);

    CSimdArray<double> prim_buffer_2_sdsf(60, npgtos);

    CSimdArray<double> prim_buffer_3_sdsf(60, npgtos);

    CSimdArray<double> prim_buffer_4_sdsf(60, npgtos);

    CSimdArray<double> prim_buffer_0_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_1_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_2_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_3_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_4_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_0_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_1_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_2_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_3_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_4_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_0_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_1_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_2_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_3_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_4_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_3_sfsp(30, npgtos);

    CSimdArray<double> prim_buffer_2_sfsd(60, npgtos);

    CSimdArray<double> prim_buffer_3_sfsd(60, npgtos);

    CSimdArray<double> prim_buffer_1_sfsf(100, npgtos);

    CSimdArray<double> prim_buffer_2_sfsf(100, npgtos);

    CSimdArray<double> prim_buffer_3_sfsf(100, npgtos);

    CSimdArray<double> prim_buffer_0_sfsg(150, npgtos);

    CSimdArray<double> prim_buffer_1_sfsg(150, npgtos);

    CSimdArray<double> prim_buffer_2_sfsg(150, npgtos);

    CSimdArray<double> prim_buffer_3_sfsg(150, npgtos);

    CSimdArray<double> prim_buffer_0_sfsh(210, npgtos);

    CSimdArray<double> prim_buffer_1_sfsh(210, npgtos);

    CSimdArray<double> prim_buffer_2_sfsh(210, npgtos);

    CSimdArray<double> prim_buffer_3_sfsh(210, npgtos);

    CSimdArray<double> prim_buffer_0_sfsi(280, npgtos);

    CSimdArray<double> prim_buffer_1_sfsi(280, npgtos);

    CSimdArray<double> prim_buffer_2_sfsi(280, npgtos);

    CSimdArray<double> prim_buffer_3_sfsi(280, npgtos);

    CSimdArray<double> prim_buffer_2_sgsd(90, npgtos);

    CSimdArray<double> prim_buffer_1_sgsf(150, npgtos);

    CSimdArray<double> prim_buffer_2_sgsf(150, npgtos);

    CSimdArray<double> prim_buffer_0_sgsg(225, npgtos);

    CSimdArray<double> prim_buffer_1_sgsg(225, npgtos);

    CSimdArray<double> prim_buffer_2_sgsg(225, npgtos);

    CSimdArray<double> prim_buffer_0_sgsh(315, npgtos);

    CSimdArray<double> prim_buffer_1_sgsh(315, npgtos);

    CSimdArray<double> prim_buffer_2_sgsh(315, npgtos);

    CSimdArray<double> prim_buffer_0_sgsi(420, npgtos);

    CSimdArray<double> prim_buffer_1_sgsi(420, npgtos);

    CSimdArray<double> prim_buffer_2_sgsi(420, npgtos);

    CSimdArray<double> prim_buffer_1_shsf(210, npgtos);

    CSimdArray<double> prim_buffer_0_shsg(315, npgtos);

    CSimdArray<double> prim_buffer_1_shsg(315, npgtos);

    CSimdArray<double> prim_buffer_0_shsh(441, npgtos);

    CSimdArray<double> prim_buffer_1_shsh(441, npgtos);

    CSimdArray<double> prim_buffer_0_shsi(588, npgtos);

    CSimdArray<double> prim_buffer_1_shsi(588, npgtos);

    CSimdArray<double> prim_buffer_0_sisg(420, npgtos);

    CSimdArray<double> prim_buffer_0_sish(588, npgtos);

    CSimdArray<double> prim_buffer_0_sisi(784, npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cart_buffer_0_sgsg(225, 1);

    CSimdArray<double> cart_buffer_0_sgsh(315, 1);

    CSimdArray<double> cart_buffer_0_sgsi(420, 1);

    CSimdArray<double> cart_buffer_0_shsg(315, 1);

    CSimdArray<double> cart_buffer_0_shsh(441, 1);

    CSimdArray<double> cart_buffer_0_shsi(588, 1);

    CSimdArray<double> cart_buffer_0_sisg(420, 1);

    CSimdArray<double> cart_buffer_0_sish(588, 1);

    CSimdArray<double> cart_buffer_0_sisi(784, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> contr_buffer_0_sgpg(675, 1);

    CSimdArray<double> contr_buffer_0_sgph(945, 1);

    CSimdArray<double> contr_buffer_0_sgdg(1350, 1);

    CSimdArray<double> contr_buffer_0_shpg(945, 1);

    CSimdArray<double> contr_buffer_0_shph(1323, 1);

    CSimdArray<double> contr_buffer_0_shdg(1890, 1);

    CSimdArray<double> contr_buffer_0_sipg(1260, 1);

    CSimdArray<double> contr_buffer_0_siph(1764, 1);

    CSimdArray<double> contr_buffer_0_sidg(2520, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> ket_spher_buffer_0_sgdg(675, 1);

    CSimdArray<double> ket_spher_buffer_0_shdg(945, 1);

    CSimdArray<double> ket_spher_buffer_0_sidg(1260, 1);

    CSimdArray<double> ket_spher_buffer_0_pgdg(2025, 1);

    CSimdArray<double> ket_spher_buffer_0_phdg(2835, 1);

    CSimdArray<double> ket_spher_buffer_0_dgdg(4050, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> spher_buffer_0_dgdg(2025, 1);

    // setup Boys fuction data

    const CBoysFunc<12> bf_table;

    CSimdArray<double> bf_args(1, npgtos);

    CSimdArray<double> bf_values(13, npgtos);

    // allocate aligned array to store max. integral values

    const auto gto_dim = gto_indices[1] - gto_indices[0];

    std::vector<double> max_values(gto_dim, 0.0);

    // loop over contracted GTOs on bra and ket sides

    for (auto i = gto_indices[0]; i < gto_indices[1]; i++)
    {
        // set up indices on ket side

        std::array<int, 2> ket_indices({i, i + 1});

        // zero integral buffers

        cart_buffer_0_sgsg.zero();

        cart_buffer_0_sgsh.zero();

        cart_buffer_0_sgsi.zero();

        cart_buffer_0_shsg.zero();

        cart_buffer_0_shsh.zero();

        cart_buffer_0_shsi.zero();

        cart_buffer_0_sisg.zero();

        cart_buffer_0_sish.zero();

        cart_buffer_0_sisi.zero();

        ket_spher_buffer_0_sgdg.zero();

        ket_spher_buffer_0_shdg.zero();

        ket_spher_buffer_0_sidg.zero();

        ket_spher_buffer_0_pgdg.zero();

        ket_spher_buffer_0_phdg.zero();

        ket_spher_buffer_0_dgdg.zero();

        spher_buffer_0_dgdg.zero();

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

         // load GTOs data for ket side

        c_x.replicate(a_coords_x, ket_indices, npgtos);

        c_y.replicate(a_coords_y, ket_indices, npgtos);

        c_z.replicate(a_coords_z, ket_indices, npgtos);

        d_x.replicate(b_coords_x, ket_indices, npgtos);

        d_y.replicate(b_coords_y, ket_indices, npgtos);

        d_z.replicate(b_coords_z, ket_indices, npgtos);

        c_exps.load(a_vec_exps, ket_indices, npgtos);

        d_exps.load(b_vec_exps, ket_indices, npgtos);

        cd_norms.load(ab_vec_norms, ket_indices, npgtos);

        cd_ovls.load(ab_vec_ovls, ket_indices, npgtos);

        // set up distances on bra side

        t4cfunc::comp_distances_cd(cd_x[0], cd_y[0], cd_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], 1);

        for (int j = 0; j < npgtos; j++)
        {
            const auto a_exp = a_vec_exps[j * ncgtos + i];

            const auto b_exp = b_vec_exps[j * ncgtos + i];

            const auto ab_norm = ab_vec_norms[j * ncgtos + i];

            const auto ab_ovl = ab_vec_ovls[j * ncgtos + i];

            const auto p_x = (a_x * a_exp + b_x * b_exp) / (a_exp + b_exp);

            const auto p_y = (a_y * a_exp + b_y * b_exp) / (a_exp + b_exp);

            const auto p_z = (a_z * a_exp + b_z * b_exp) / (a_exp + b_exp);

            const auto pb_x = p_x - b_x;

            const auto pb_y = p_y - b_y;

            const auto pb_z = p_z - b_z;

            t4cfunc::comp_coordinates_q(q_x[0], q_y[0], q_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], c_exps[0], d_exps[0], npgtos);

            t4cfunc::comp_coordinates_w(w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], a_exp, b_exp, c_exps[0], d_exps[0], npgtos);

            t4cfunc::comp_distances_pq(pq_x[0], pq_y[0], pq_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], npgtos);

            t4cfunc::comp_distances_wq(wq_x[0], wq_y[0], wq_z[0], w_x[0], w_y[0], w_z[0], q_x[0], q_y[0], q_z[0], npgtos);

            t4cfunc::comp_distances_qd(qd_x[0], qd_y[0], qd_z[0], q_x[0], q_y[0], q_z[0], d_x[0], d_y[0], d_z[0], npgtos);

            t4cfunc::comp_distances_wp(wp_x[0], wp_y[0], wp_z[0], w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, npgtos);

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

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_11_ssss, fss_abcd[0], bf_values[11]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_12_ssss, fss_abcd[0], bf_values[12]);

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

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_10_sssp, prim_buffer_10_ssss, prim_buffer_11_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_11_sssp, prim_buffer_11_ssss, prim_buffer_12_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_0_sssd, prim_buffer_0_ssss, prim_buffer_1_ssss, prim_buffer_0_sssp, prim_buffer_1_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_1_sssd, prim_buffer_1_ssss, prim_buffer_2_ssss, prim_buffer_1_sssp, prim_buffer_2_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_2_sssd, prim_buffer_2_ssss, prim_buffer_3_ssss, prim_buffer_2_sssp, prim_buffer_3_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_3_sssd, prim_buffer_3_ssss, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_4_sssd, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_5_sssd, prim_buffer_5_ssss, prim_buffer_6_ssss, prim_buffer_5_sssp, prim_buffer_6_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_6_sssd, prim_buffer_6_ssss, prim_buffer_7_ssss, prim_buffer_6_sssp, prim_buffer_7_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_7_sssd, prim_buffer_7_ssss, prim_buffer_8_ssss, prim_buffer_7_sssp, prim_buffer_8_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_8_sssd, prim_buffer_8_ssss, prim_buffer_9_ssss, prim_buffer_8_sssp, prim_buffer_9_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_9_sssd, prim_buffer_9_ssss, prim_buffer_10_ssss, prim_buffer_9_sssp, prim_buffer_10_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_10_sssd, prim_buffer_10_ssss, prim_buffer_11_ssss, prim_buffer_10_sssp, prim_buffer_11_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_0_sssf, prim_buffer_0_sssp, prim_buffer_1_sssp, prim_buffer_0_sssd, prim_buffer_1_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_1_sssf, prim_buffer_1_sssp, prim_buffer_2_sssp, prim_buffer_1_sssd, prim_buffer_2_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_2_sssf, prim_buffer_2_sssp, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_3_sssf, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_3_sssd, prim_buffer_4_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_4_sssf, prim_buffer_4_sssp, prim_buffer_5_sssp, prim_buffer_4_sssd, prim_buffer_5_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_5_sssf, prim_buffer_5_sssp, prim_buffer_6_sssp, prim_buffer_5_sssd, prim_buffer_6_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_6_sssf, prim_buffer_6_sssp, prim_buffer_7_sssp, prim_buffer_6_sssd, prim_buffer_7_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_7_sssf, prim_buffer_7_sssp, prim_buffer_8_sssp, prim_buffer_7_sssd, prim_buffer_8_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_8_sssf, prim_buffer_8_sssp, prim_buffer_9_sssp, prim_buffer_8_sssd, prim_buffer_9_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_9_sssf, prim_buffer_9_sssp, prim_buffer_10_sssp, prim_buffer_9_sssd, prim_buffer_10_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_0_sssg, prim_buffer_0_sssd, prim_buffer_1_sssd, prim_buffer_0_sssf, prim_buffer_1_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_1_sssg, prim_buffer_1_sssd, prim_buffer_2_sssd, prim_buffer_1_sssf, prim_buffer_2_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_2_sssg, prim_buffer_2_sssd, prim_buffer_3_sssd, prim_buffer_2_sssf, prim_buffer_3_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_3_sssg, prim_buffer_3_sssd, prim_buffer_4_sssd, prim_buffer_3_sssf, prim_buffer_4_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_4_sssg, prim_buffer_4_sssd, prim_buffer_5_sssd, prim_buffer_4_sssf, prim_buffer_5_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_5_sssg, prim_buffer_5_sssd, prim_buffer_6_sssd, prim_buffer_5_sssf, prim_buffer_6_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_6_sssg, prim_buffer_6_sssd, prim_buffer_7_sssd, prim_buffer_6_sssf, prim_buffer_7_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_7_sssg, prim_buffer_7_sssd, prim_buffer_8_sssd, prim_buffer_7_sssf, prim_buffer_8_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_8_sssg, prim_buffer_8_sssd, prim_buffer_9_sssd, prim_buffer_8_sssf, prim_buffer_9_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_0_sssh, prim_buffer_0_sssf, prim_buffer_1_sssf, prim_buffer_0_sssg, prim_buffer_1_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_1_sssh, prim_buffer_1_sssf, prim_buffer_2_sssf, prim_buffer_1_sssg, prim_buffer_2_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_2_sssh, prim_buffer_2_sssf, prim_buffer_3_sssf, prim_buffer_2_sssg, prim_buffer_3_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_3_sssh, prim_buffer_3_sssf, prim_buffer_4_sssf, prim_buffer_3_sssg, prim_buffer_4_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_4_sssh, prim_buffer_4_sssf, prim_buffer_5_sssf, prim_buffer_4_sssg, prim_buffer_5_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_5_sssh, prim_buffer_5_sssf, prim_buffer_6_sssf, prim_buffer_5_sssg, prim_buffer_6_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_6_sssh, prim_buffer_6_sssf, prim_buffer_7_sssf, prim_buffer_6_sssg, prim_buffer_7_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_7_sssh, prim_buffer_7_sssf, prim_buffer_8_sssf, prim_buffer_7_sssg, prim_buffer_8_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_0_sssi, prim_buffer_0_sssg, prim_buffer_1_sssg, prim_buffer_0_sssh, prim_buffer_1_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_1_sssi, prim_buffer_1_sssg, prim_buffer_2_sssg, prim_buffer_1_sssh, prim_buffer_2_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_2_sssi, prim_buffer_2_sssg, prim_buffer_3_sssg, prim_buffer_2_sssh, prim_buffer_3_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_3_sssi, prim_buffer_3_sssg, prim_buffer_4_sssg, prim_buffer_3_sssh, prim_buffer_4_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_4_sssi, prim_buffer_4_sssg, prim_buffer_5_sssg, prim_buffer_4_sssh, prim_buffer_5_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_5_sssi, prim_buffer_5_sssg, prim_buffer_6_sssg, prim_buffer_5_sssh, prim_buffer_6_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_6_sssi, prim_buffer_6_sssg, prim_buffer_7_sssg, prim_buffer_6_sssh, prim_buffer_7_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_4_spss, prim_buffer_4_ssss, prim_buffer_5_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_5_spss, prim_buffer_5_ssss, prim_buffer_6_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_3_spsp, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_4_spsp, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_5_spsp, prim_buffer_6_ssss, prim_buffer_5_sssp, prim_buffer_6_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_2_spsd, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_3_spsd, prim_buffer_4_sssp, prim_buffer_3_sssd, prim_buffer_4_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_4_spsd, prim_buffer_5_sssp, prim_buffer_4_sssd, prim_buffer_5_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_5_spsd, prim_buffer_6_sssp, prim_buffer_5_sssd, prim_buffer_6_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_1_spsf, prim_buffer_2_sssd, prim_buffer_1_sssf, prim_buffer_2_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_2_spsf, prim_buffer_3_sssd, prim_buffer_2_sssf, prim_buffer_3_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_3_spsf, prim_buffer_4_sssd, prim_buffer_3_sssf, prim_buffer_4_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_4_spsf, prim_buffer_5_sssd, prim_buffer_4_sssf, prim_buffer_5_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_5_spsf, prim_buffer_6_sssd, prim_buffer_5_sssf, prim_buffer_6_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_0_spsg, prim_buffer_1_sssf, prim_buffer_0_sssg, prim_buffer_1_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_1_spsg, prim_buffer_2_sssf, prim_buffer_1_sssg, prim_buffer_2_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_2_spsg, prim_buffer_3_sssf, prim_buffer_2_sssg, prim_buffer_3_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_3_spsg, prim_buffer_4_sssf, prim_buffer_3_sssg, prim_buffer_4_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_4_spsg, prim_buffer_5_sssf, prim_buffer_4_sssg, prim_buffer_5_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_5_spsg, prim_buffer_6_sssf, prim_buffer_5_sssg, prim_buffer_6_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_0_spsh, prim_buffer_1_sssg, prim_buffer_0_sssh, prim_buffer_1_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_1_spsh, prim_buffer_2_sssg, prim_buffer_1_sssh, prim_buffer_2_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_2_spsh, prim_buffer_3_sssg, prim_buffer_2_sssh, prim_buffer_3_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_3_spsh, prim_buffer_4_sssg, prim_buffer_3_sssh, prim_buffer_4_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_4_spsh, prim_buffer_5_sssg, prim_buffer_4_sssh, prim_buffer_5_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_5_spsh, prim_buffer_6_sssg, prim_buffer_5_sssh, prim_buffer_6_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_0_spsi, prim_buffer_1_sssh, prim_buffer_0_sssi, prim_buffer_1_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_1_spsi, prim_buffer_2_sssh, prim_buffer_1_sssi, prim_buffer_2_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_2_spsi, prim_buffer_3_sssh, prim_buffer_2_sssi, prim_buffer_3_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_3_spsi, prim_buffer_4_sssh, prim_buffer_3_sssi, prim_buffer_4_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_4_spsi, prim_buffer_5_sssh, prim_buffer_4_sssi, prim_buffer_5_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_5_spsi, prim_buffer_6_sssh, prim_buffer_5_sssi, prim_buffer_6_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_4_sdss, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_spss, prim_buffer_5_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_3_sdsp, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_4_spss, prim_buffer_3_spsp, prim_buffer_4_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_4_sdsp, prim_buffer_4_sssp, prim_buffer_5_sssp, prim_buffer_5_spss, prim_buffer_4_spsp, prim_buffer_5_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_2_sdsd, prim_buffer_2_sssd, prim_buffer_3_sssd, prim_buffer_3_spsp, prim_buffer_2_spsd, prim_buffer_3_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_3_sdsd, prim_buffer_3_sssd, prim_buffer_4_sssd, prim_buffer_4_spsp, prim_buffer_3_spsd, prim_buffer_4_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_4_sdsd, prim_buffer_4_sssd, prim_buffer_5_sssd, prim_buffer_5_spsp, prim_buffer_4_spsd, prim_buffer_5_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_1_sdsf, prim_buffer_1_sssf, prim_buffer_2_sssf, prim_buffer_2_spsd, prim_buffer_1_spsf, prim_buffer_2_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_2_sdsf, prim_buffer_2_sssf, prim_buffer_3_sssf, prim_buffer_3_spsd, prim_buffer_2_spsf, prim_buffer_3_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_3_sdsf, prim_buffer_3_sssf, prim_buffer_4_sssf, prim_buffer_4_spsd, prim_buffer_3_spsf, prim_buffer_4_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_4_sdsf, prim_buffer_4_sssf, prim_buffer_5_sssf, prim_buffer_5_spsd, prim_buffer_4_spsf, prim_buffer_5_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_0_sdsg, prim_buffer_0_sssg, prim_buffer_1_sssg, prim_buffer_1_spsf, prim_buffer_0_spsg, prim_buffer_1_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_1_sdsg, prim_buffer_1_sssg, prim_buffer_2_sssg, prim_buffer_2_spsf, prim_buffer_1_spsg, prim_buffer_2_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_2_sdsg, prim_buffer_2_sssg, prim_buffer_3_sssg, prim_buffer_3_spsf, prim_buffer_2_spsg, prim_buffer_3_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_3_sdsg, prim_buffer_3_sssg, prim_buffer_4_sssg, prim_buffer_4_spsf, prim_buffer_3_spsg, prim_buffer_4_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_4_sdsg, prim_buffer_4_sssg, prim_buffer_5_sssg, prim_buffer_5_spsf, prim_buffer_4_spsg, prim_buffer_5_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_0_sdsh, prim_buffer_0_sssh, prim_buffer_1_sssh, prim_buffer_1_spsg, prim_buffer_0_spsh, prim_buffer_1_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_1_sdsh, prim_buffer_1_sssh, prim_buffer_2_sssh, prim_buffer_2_spsg, prim_buffer_1_spsh, prim_buffer_2_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_2_sdsh, prim_buffer_2_sssh, prim_buffer_3_sssh, prim_buffer_3_spsg, prim_buffer_2_spsh, prim_buffer_3_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_3_sdsh, prim_buffer_3_sssh, prim_buffer_4_sssh, prim_buffer_4_spsg, prim_buffer_3_spsh, prim_buffer_4_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_4_sdsh, prim_buffer_4_sssh, prim_buffer_5_sssh, prim_buffer_5_spsg, prim_buffer_4_spsh, prim_buffer_5_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_0_sdsi, prim_buffer_0_sssi, prim_buffer_1_sssi, prim_buffer_1_spsh, prim_buffer_0_spsi, prim_buffer_1_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_1_sdsi, prim_buffer_1_sssi, prim_buffer_2_sssi, prim_buffer_2_spsh, prim_buffer_1_spsi, prim_buffer_2_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_2_sdsi, prim_buffer_2_sssi, prim_buffer_3_sssi, prim_buffer_3_spsh, prim_buffer_2_spsi, prim_buffer_3_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_3_sdsi, prim_buffer_3_sssi, prim_buffer_4_sssi, prim_buffer_4_spsh, prim_buffer_3_spsi, prim_buffer_4_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_4_sdsi, prim_buffer_4_sssi, prim_buffer_5_sssi, prim_buffer_5_spsh, prim_buffer_4_spsi, prim_buffer_5_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_3_sfsp, prim_buffer_3_spsp, prim_buffer_4_spsp, prim_buffer_4_sdss, prim_buffer_3_sdsp, prim_buffer_4_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_2_sfsd, prim_buffer_2_spsd, prim_buffer_3_spsd, prim_buffer_3_sdsp, prim_buffer_2_sdsd, prim_buffer_3_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_3_sfsd, prim_buffer_3_spsd, prim_buffer_4_spsd, prim_buffer_4_sdsp, prim_buffer_3_sdsd, prim_buffer_4_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_1_sfsf, prim_buffer_1_spsf, prim_buffer_2_spsf, prim_buffer_2_sdsd, prim_buffer_1_sdsf, prim_buffer_2_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_2_sfsf, prim_buffer_2_spsf, prim_buffer_3_spsf, prim_buffer_3_sdsd, prim_buffer_2_sdsf, prim_buffer_3_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_3_sfsf, prim_buffer_3_spsf, prim_buffer_4_spsf, prim_buffer_4_sdsd, prim_buffer_3_sdsf, prim_buffer_4_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_0_sfsg, prim_buffer_0_spsg, prim_buffer_1_spsg, prim_buffer_1_sdsf, prim_buffer_0_sdsg, prim_buffer_1_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_1_sfsg, prim_buffer_1_spsg, prim_buffer_2_spsg, prim_buffer_2_sdsf, prim_buffer_1_sdsg, prim_buffer_2_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_2_sfsg, prim_buffer_2_spsg, prim_buffer_3_spsg, prim_buffer_3_sdsf, prim_buffer_2_sdsg, prim_buffer_3_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_3_sfsg, prim_buffer_3_spsg, prim_buffer_4_spsg, prim_buffer_4_sdsf, prim_buffer_3_sdsg, prim_buffer_4_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_0_sfsh, prim_buffer_0_spsh, prim_buffer_1_spsh, prim_buffer_1_sdsg, prim_buffer_0_sdsh, prim_buffer_1_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_1_sfsh, prim_buffer_1_spsh, prim_buffer_2_spsh, prim_buffer_2_sdsg, prim_buffer_1_sdsh, prim_buffer_2_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_2_sfsh, prim_buffer_2_spsh, prim_buffer_3_spsh, prim_buffer_3_sdsg, prim_buffer_2_sdsh, prim_buffer_3_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_3_sfsh, prim_buffer_3_spsh, prim_buffer_4_spsh, prim_buffer_4_sdsg, prim_buffer_3_sdsh, prim_buffer_4_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_0_sfsi, prim_buffer_0_spsi, prim_buffer_1_spsi, prim_buffer_1_sdsh, prim_buffer_0_sdsi, prim_buffer_1_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_1_sfsi, prim_buffer_1_spsi, prim_buffer_2_spsi, prim_buffer_2_sdsh, prim_buffer_1_sdsi, prim_buffer_2_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_2_sfsi, prim_buffer_2_spsi, prim_buffer_3_spsi, prim_buffer_3_sdsh, prim_buffer_2_sdsi, prim_buffer_3_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_3_sfsi, prim_buffer_3_spsi, prim_buffer_4_spsi, prim_buffer_4_sdsh, prim_buffer_3_sdsi, prim_buffer_4_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsd(prim_buffer_2_sgsd, prim_buffer_2_sdsd, prim_buffer_3_sdsd, prim_buffer_3_sfsp, prim_buffer_2_sfsd, prim_buffer_3_sfsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsf(prim_buffer_1_sgsf, prim_buffer_1_sdsf, prim_buffer_2_sdsf, prim_buffer_2_sfsd, prim_buffer_1_sfsf, prim_buffer_2_sfsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsf(prim_buffer_2_sgsf, prim_buffer_2_sdsf, prim_buffer_3_sdsf, prim_buffer_3_sfsd, prim_buffer_2_sfsf, prim_buffer_3_sfsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_0_sgsg, prim_buffer_0_sdsg, prim_buffer_1_sdsg, prim_buffer_1_sfsf, prim_buffer_0_sfsg, prim_buffer_1_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_1_sgsg, prim_buffer_1_sdsg, prim_buffer_2_sdsg, prim_buffer_2_sfsf, prim_buffer_1_sfsg, prim_buffer_2_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_2_sgsg, prim_buffer_2_sdsg, prim_buffer_3_sdsg, prim_buffer_3_sfsf, prim_buffer_2_sfsg, prim_buffer_3_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_0_sgsh, prim_buffer_0_sdsh, prim_buffer_1_sdsh, prim_buffer_1_sfsg, prim_buffer_0_sfsh, prim_buffer_1_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_1_sgsh, prim_buffer_1_sdsh, prim_buffer_2_sdsh, prim_buffer_2_sfsg, prim_buffer_1_sfsh, prim_buffer_2_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_2_sgsh, prim_buffer_2_sdsh, prim_buffer_3_sdsh, prim_buffer_3_sfsg, prim_buffer_2_sfsh, prim_buffer_3_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_0_sgsi, prim_buffer_0_sdsi, prim_buffer_1_sdsi, prim_buffer_1_sfsh, prim_buffer_0_sfsi, prim_buffer_1_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_1_sgsi, prim_buffer_1_sdsi, prim_buffer_2_sdsi, prim_buffer_2_sfsh, prim_buffer_1_sfsi, prim_buffer_2_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_2_sgsi, prim_buffer_2_sdsi, prim_buffer_3_sdsi, prim_buffer_3_sfsh, prim_buffer_2_sfsi, prim_buffer_3_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsf(prim_buffer_1_shsf, prim_buffer_1_sfsf, prim_buffer_2_sfsf, prim_buffer_2_sgsd, prim_buffer_1_sgsf, prim_buffer_2_sgsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsg(prim_buffer_0_shsg, prim_buffer_0_sfsg, prim_buffer_1_sfsg, prim_buffer_1_sgsf, prim_buffer_0_sgsg, prim_buffer_1_sgsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsg(prim_buffer_1_shsg, prim_buffer_1_sfsg, prim_buffer_2_sfsg, prim_buffer_2_sgsf, prim_buffer_1_sgsg, prim_buffer_2_sgsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsh(prim_buffer_0_shsh, prim_buffer_0_sfsh, prim_buffer_1_sfsh, prim_buffer_1_sgsg, prim_buffer_0_sgsh, prim_buffer_1_sgsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsh(prim_buffer_1_shsh, prim_buffer_1_sfsh, prim_buffer_2_sfsh, prim_buffer_2_sgsg, prim_buffer_1_sgsh, prim_buffer_2_sgsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsi(prim_buffer_0_shsi, prim_buffer_0_sfsi, prim_buffer_1_sfsi, prim_buffer_1_sgsh, prim_buffer_0_sgsi, prim_buffer_1_sgsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsi(prim_buffer_1_shsi, prim_buffer_1_sfsi, prim_buffer_2_sfsi, prim_buffer_2_sgsh, prim_buffer_1_sgsi, prim_buffer_2_sgsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisg(prim_buffer_0_sisg, prim_buffer_0_sgsg, prim_buffer_1_sgsg, prim_buffer_1_shsf, prim_buffer_0_shsg, prim_buffer_1_shsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sish(prim_buffer_0_sish, prim_buffer_0_sgsh, prim_buffer_1_sgsh, prim_buffer_1_shsg, prim_buffer_0_shsh, prim_buffer_1_shsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisi(prim_buffer_0_sisi, prim_buffer_0_sgsi, prim_buffer_1_sgsi, prim_buffer_1_shsh, prim_buffer_0_shsi, prim_buffer_1_shsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            t2cfunc::reduce(cart_buffer_0_sgsg, prim_buffer_0_sgsg, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sgsh, prim_buffer_0_sgsh, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sgsi, prim_buffer_0_sgsi, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_shsg, prim_buffer_0_shsg, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_shsh, prim_buffer_0_shsh, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_shsi, prim_buffer_0_shsi, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sisg, prim_buffer_0_sisg, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sish, prim_buffer_0_sish, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sisi, prim_buffer_0_sisi, 1, npgtos);

        }

        erirec::comp_ket_hrr_electron_repulsion_xxpg(contr_buffer_0_sgpg, cart_buffer_0_sgsg, cart_buffer_0_sgsh, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxph(contr_buffer_0_sgph, cart_buffer_0_sgsh, cart_buffer_0_sgsi, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(contr_buffer_0_sgdg, contr_buffer_0_sgpg, contr_buffer_0_sgph, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(contr_buffer_0_shpg, cart_buffer_0_shsg, cart_buffer_0_shsh, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxph(contr_buffer_0_shph, cart_buffer_0_shsh, cart_buffer_0_shsi, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(contr_buffer_0_shdg, contr_buffer_0_shpg, contr_buffer_0_shph, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(contr_buffer_0_sipg, cart_buffer_0_sisg, cart_buffer_0_sish, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxph(contr_buffer_0_siph, cart_buffer_0_sish, cart_buffer_0_sisi, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(contr_buffer_0_sidg, contr_buffer_0_sipg, contr_buffer_0_siph, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        t4cfunc::ket_transform<2, 4>(ket_spher_buffer_0_sgdg, contr_buffer_0_sgdg, 0, 4);

        t4cfunc::ket_transform<2, 4>(ket_spher_buffer_0_shdg, contr_buffer_0_shdg, 0, 5);

        t4cfunc::ket_transform<2, 4>(ket_spher_buffer_0_sidg, contr_buffer_0_sidg, 0, 6);

        erirec::comp_bra_hrr_electron_repulsion_pgxx(ket_spher_buffer_0_pgdg, ket_spher_buffer_0_sgdg, ket_spher_buffer_0_shdg, ab_x, ab_y, ab_z, 2, 4);

        erirec::comp_bra_hrr_electron_repulsion_phxx(ket_spher_buffer_0_phdg, ket_spher_buffer_0_shdg, ket_spher_buffer_0_sidg, ab_x, ab_y, ab_z, 2, 4);

        erirec::comp_bra_hrr_electron_repulsion_dgxx(ket_spher_buffer_0_dgdg, ket_spher_buffer_0_pgdg, ket_spher_buffer_0_phdg, ab_x, ab_y, ab_z, 2, 4);

        t4cfunc::bra_transform<2, 4>(spher_buffer_0_dgdg, ket_spher_buffer_0_dgdg, 2, 4);

        t4cfunc::update_max_values(max_values, spher_buffer_0_dgdg, i - gto_indices[0]); 
    }

    distributor->distribute(max_values, gto_indices);
}

} // erirec namespace

#endif /* ElectronRepulsionDiagRecDGDG_hpp */
