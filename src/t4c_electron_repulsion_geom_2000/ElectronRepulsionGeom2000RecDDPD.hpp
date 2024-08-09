#ifndef ElectronRepulsionGeom2000RecDDPD_hpp
#define ElectronRepulsionGeom2000RecDDPD_hpp

#include <array>

#include "ElectronRepulsionPrimRecDDPD.hpp"
#include "ElectronRepulsionPrimRecDDPP.hpp"
#include "ElectronRepulsionPrimRecDDPS.hpp"
#include "ElectronRepulsionPrimRecDDSD.hpp"
#include "ElectronRepulsionPrimRecDDSP.hpp"
#include "ElectronRepulsionPrimRecDPPD.hpp"
#include "ElectronRepulsionPrimRecDPPP.hpp"
#include "ElectronRepulsionPrimRecDPSD.hpp"
#include "ElectronRepulsionPrimRecDSPD.hpp"
#include "ElectronRepulsionPrimRecFDPD.hpp"
#include "ElectronRepulsionPrimRecFDPP.hpp"
#include "ElectronRepulsionPrimRecFDSD.hpp"
#include "ElectronRepulsionPrimRecFPPD.hpp"
#include "ElectronRepulsionPrimRecGDPD.hpp"
#include "ElectronRepulsionPrimRecPDPD.hpp"
#include "ElectronRepulsionPrimRecPDPP.hpp"
#include "ElectronRepulsionPrimRecPDPS.hpp"
#include "ElectronRepulsionPrimRecPDSD.hpp"
#include "ElectronRepulsionPrimRecPDSP.hpp"
#include "ElectronRepulsionPrimRecPDSS.hpp"
#include "ElectronRepulsionPrimRecPPPD.hpp"
#include "ElectronRepulsionPrimRecPPPP.hpp"
#include "ElectronRepulsionPrimRecPPPS.hpp"
#include "ElectronRepulsionPrimRecPPSD.hpp"
#include "ElectronRepulsionPrimRecPPSP.hpp"
#include "ElectronRepulsionPrimRecPSPD.hpp"
#include "ElectronRepulsionPrimRecPSPP.hpp"
#include "ElectronRepulsionPrimRecPSSD.hpp"
#include "ElectronRepulsionPrimRecSDPD.hpp"
#include "ElectronRepulsionPrimRecSDPP.hpp"
#include "ElectronRepulsionPrimRecSDPS.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSPPD.hpp"
#include "ElectronRepulsionPrimRecSPPP.hpp"
#include "ElectronRepulsionPrimRecSPPS.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSPD.hpp"
#include "ElectronRepulsionPrimRecSSPP.hpp"
#include "ElectronRepulsionPrimRecSSPS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "GeomDeriv2000OfScalarForDDPD.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// Computes (DD|1/|r-r'||PD)  integrals for two GTOs pair blocks.
/// - Parameter distributor: the pointer to Fock matrix/matrices distributor.
/// - Parameter bra_gto_pair_block: the GTOs pair block on bra side.
/// - Parameter ket_gto_pair_block: the GTOs pair block on ket side.
/// - Parameter bra_indices: the range [bra_first, bra_last) of GTOs on bra side.
/// - Parameter ket_indices: the range [ket_first, ket_last) of GTOs on ket side.
template <class T>
auto
comp_electron_repulsion_geom2000_ddpd(T* distributor,
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

    CSimdArray<double> prim_buffer_2_ssps(3, ket_pdim);

    CSimdArray<double> prim_buffer_3_ssps(3, ket_pdim);

    CSimdArray<double> prim_buffer_4_ssps(3, ket_pdim);

    CSimdArray<double> prim_buffer_5_ssps(3, ket_pdim);

    CSimdArray<double> prim_buffer_6_ssps(3, ket_pdim);

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

    CSimdArray<double> prim_buffer_3_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_4_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_5_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_spps(9, ket_pdim);

    CSimdArray<double> prim_buffer_3_spps(9, ket_pdim);

    CSimdArray<double> prim_buffer_4_spps(9, ket_pdim);

    CSimdArray<double> prim_buffer_5_spps(9, ket_pdim);

    CSimdArray<double> prim_buffer_1_sppp(27, ket_pdim);

    CSimdArray<double> prim_buffer_2_sppp(27, ket_pdim);

    CSimdArray<double> prim_buffer_3_sppp(27, ket_pdim);

    CSimdArray<double> prim_buffer_4_sppp(27, ket_pdim);

    CSimdArray<double> prim_buffer_5_sppp(27, ket_pdim);

    CSimdArray<double> prim_buffer_0_sppd(54, ket_pdim);

    CSimdArray<double> prim_buffer_1_sppd(54, ket_pdim);

    CSimdArray<double> prim_buffer_2_sppd(54, ket_pdim);

    CSimdArray<double> prim_buffer_3_sppd(54, ket_pdim);

    CSimdArray<double> prim_buffer_4_sppd(54, ket_pdim);

    CSimdArray<double> prim_buffer_5_sppd(54, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdss(6, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdss(6, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdps(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdps(18, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdps(18, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdpp(54, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdpp(54, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdpp(54, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdpp(54, ket_pdim);

    CSimdArray<double> prim_buffer_0_sdpd(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdpd(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdpd(108, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdpd(108, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdpd(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_pssd(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_pssd(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_pssd(18, ket_pdim);

    CSimdArray<double> prim_buffer_1_pspp(27, ket_pdim);

    CSimdArray<double> prim_buffer_2_pspp(27, ket_pdim);

    CSimdArray<double> prim_buffer_3_pspp(27, ket_pdim);

    CSimdArray<double> prim_buffer_0_pspd(54, ket_pdim);

    CSimdArray<double> prim_buffer_1_pspd(54, ket_pdim);

    CSimdArray<double> prim_buffer_2_pspd(54, ket_pdim);

    CSimdArray<double> prim_buffer_3_pspd(54, ket_pdim);

    CSimdArray<double> prim_buffer_2_ppsp(27, ket_pdim);

    CSimdArray<double> prim_buffer_3_ppsp(27, ket_pdim);

    CSimdArray<double> prim_buffer_1_ppsd(54, ket_pdim);

    CSimdArray<double> prim_buffer_2_ppsd(54, ket_pdim);

    CSimdArray<double> prim_buffer_3_ppsd(54, ket_pdim);

    CSimdArray<double> prim_buffer_2_ppps(27, ket_pdim);

    CSimdArray<double> prim_buffer_3_ppps(27, ket_pdim);

    CSimdArray<double> prim_buffer_1_pppp(81, ket_pdim);

    CSimdArray<double> prim_buffer_2_pppp(81, ket_pdim);

    CSimdArray<double> prim_buffer_3_pppp(81, ket_pdim);

    CSimdArray<double> prim_buffer_0_pppd(162, ket_pdim);

    CSimdArray<double> prim_buffer_1_pppd(162, ket_pdim);

    CSimdArray<double> prim_buffer_2_pppd(162, ket_pdim);

    CSimdArray<double> prim_buffer_3_pppd(162, ket_pdim);

    CSimdArray<double> prim_buffer_3_pdss(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_pdsp(54, ket_pdim);

    CSimdArray<double> prim_buffer_3_pdsp(54, ket_pdim);

    CSimdArray<double> prim_buffer_1_pdsd(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_pdsd(108, ket_pdim);

    CSimdArray<double> prim_buffer_3_pdsd(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_pdps(54, ket_pdim);

    CSimdArray<double> prim_buffer_3_pdps(54, ket_pdim);

    CSimdArray<double> prim_buffer_1_pdpp(162, ket_pdim);

    CSimdArray<double> prim_buffer_2_pdpp(162, ket_pdim);

    CSimdArray<double> prim_buffer_3_pdpp(162, ket_pdim);

    CSimdArray<double> prim_buffer_0_pdpd(324, ket_pdim);

    CSimdArray<double> prim_buffer_1_pdpd(324, ket_pdim);

    CSimdArray<double> prim_buffer_2_pdpd(324, ket_pdim);

    CSimdArray<double> prim_buffer_3_pdpd(324, ket_pdim);

    CSimdArray<double> prim_buffer_0_dspd(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_dspd(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_dspd(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_dpsd(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_dpsd(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_dppp(162, ket_pdim);

    CSimdArray<double> prim_buffer_2_dppp(162, ket_pdim);

    CSimdArray<double> prim_buffer_0_dppd(324, ket_pdim);

    CSimdArray<double> prim_buffer_1_dppd(324, ket_pdim);

    CSimdArray<double> prim_buffer_2_dppd(324, ket_pdim);

    CSimdArray<double> prim_buffer_2_ddsp(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_ddsd(216, ket_pdim);

    CSimdArray<double> prim_buffer_2_ddsd(216, ket_pdim);

    CSimdArray<double> prim_buffer_2_ddps(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_ddpp(324, ket_pdim);

    CSimdArray<double> prim_buffer_2_ddpp(324, ket_pdim);

    CSimdArray<double> prim_buffer_0_ddpd(648, ket_pdim);

    CSimdArray<double> prim_buffer_1_ddpd(648, ket_pdim);

    CSimdArray<double> prim_buffer_2_ddpd(648, ket_pdim);

    CSimdArray<double> prim_buffer_0_fppd(540, ket_pdim);

    CSimdArray<double> prim_buffer_1_fppd(540, ket_pdim);

    CSimdArray<double> prim_buffer_1_fdsd(360, ket_pdim);

    CSimdArray<double> prim_buffer_1_fdpp(540, ket_pdim);

    CSimdArray<double> prim_buffer_0_fdpd(1080, ket_pdim);

    CSimdArray<double> prim_buffer_1_fdpd(1080, ket_pdim);

    CSimdArray<double> prim_buffer_0_gdpd(1620, ket_pdim);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cart_buffer_0_geom2000_ddpd(3888, ket_dim);

    // allocate aligned spherical integrals

    CSimdArray<double> spher_buffer_0_geom2000_ddpd(2250, ket_dim);

    // setup Boys fuction data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_args(1, ket_pdim);

    CSimdArray<double> bf_values(10, ket_pdim);

    // loop over contracted GTOs on bra side

    for (auto i = bra_indices[0]; i < bra_indices[1]; i++)
    {
        cart_buffer_0_geom2000_ddpd.zero();

        spher_buffer_0_geom2000_ddpd.zero();

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

            erirec::comp_prim_electron_repulsion_ssps(prim_buffer_2_ssps, prim_buffer_2_ssss, prim_buffer_3_ssss, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_ssps(prim_buffer_3_ssps, prim_buffer_3_ssss, prim_buffer_4_ssss, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_ssps(prim_buffer_4_ssps, prim_buffer_4_ssss, prim_buffer_5_ssss, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_ssps(prim_buffer_5_ssps, prim_buffer_5_ssss, prim_buffer_6_ssss, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_ssps(prim_buffer_6_ssps, prim_buffer_6_ssss, prim_buffer_7_ssss, qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]);

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

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_3_spss, prim_buffer_3_ssss, prim_buffer_4_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_4_spss, prim_buffer_4_ssss, prim_buffer_5_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_5_spss, prim_buffer_5_ssss, prim_buffer_6_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_2_spsp, prim_buffer_3_ssss, prim_buffer_2_sssp, prim_buffer_3_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_3_spsp, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_4_spsp, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_5_spsp, prim_buffer_6_ssss, prim_buffer_5_sssp, prim_buffer_6_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_1_spsd, prim_buffer_2_sssp, prim_buffer_1_sssd, prim_buffer_2_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_2_spsd, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_3_spsd, prim_buffer_4_sssp, prim_buffer_3_sssd, prim_buffer_4_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_4_spsd, prim_buffer_5_sssp, prim_buffer_4_sssd, prim_buffer_5_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_5_spsd, prim_buffer_6_sssp, prim_buffer_5_sssd, prim_buffer_6_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spps(prim_buffer_2_spps, prim_buffer_3_ssss, prim_buffer_2_ssps, prim_buffer_3_ssps, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spps(prim_buffer_3_spps, prim_buffer_4_ssss, prim_buffer_3_ssps, prim_buffer_4_ssps, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spps(prim_buffer_4_spps, prim_buffer_5_ssss, prim_buffer_4_ssps, prim_buffer_5_ssps, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spps(prim_buffer_5_spps, prim_buffer_6_ssss, prim_buffer_5_ssps, prim_buffer_6_ssps, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppp(prim_buffer_1_sppp, prim_buffer_2_sssp, prim_buffer_2_ssps, prim_buffer_1_sspp, prim_buffer_2_sspp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppp(prim_buffer_2_sppp, prim_buffer_3_sssp, prim_buffer_3_ssps, prim_buffer_2_sspp, prim_buffer_3_sspp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppp(prim_buffer_3_sppp, prim_buffer_4_sssp, prim_buffer_4_ssps, prim_buffer_3_sspp, prim_buffer_4_sspp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppp(prim_buffer_4_sppp, prim_buffer_5_sssp, prim_buffer_5_ssps, prim_buffer_4_sspp, prim_buffer_5_sspp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppp(prim_buffer_5_sppp, prim_buffer_6_sssp, prim_buffer_6_ssps, prim_buffer_5_sspp, prim_buffer_6_sspp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppd(prim_buffer_0_sppd, prim_buffer_1_sssd, prim_buffer_1_sspp, prim_buffer_0_sspd, prim_buffer_1_sspd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppd(prim_buffer_1_sppd, prim_buffer_2_sssd, prim_buffer_2_sspp, prim_buffer_1_sspd, prim_buffer_2_sspd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppd(prim_buffer_2_sppd, prim_buffer_3_sssd, prim_buffer_3_sspp, prim_buffer_2_sspd, prim_buffer_3_sspd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppd(prim_buffer_3_sppd, prim_buffer_4_sssd, prim_buffer_4_sspp, prim_buffer_3_sspd, prim_buffer_4_sspd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppd(prim_buffer_4_sppd, prim_buffer_5_sssd, prim_buffer_5_sspp, prim_buffer_4_sspd, prim_buffer_5_sspd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sppd(prim_buffer_5_sppd, prim_buffer_6_sssd, prim_buffer_6_sspp, prim_buffer_5_sspd, prim_buffer_6_sspd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_3_sdss, prim_buffer_3_ssss, prim_buffer_4_ssss, prim_buffer_3_spss, prim_buffer_4_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_4_sdss, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_spss, prim_buffer_5_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_2_sdsp, prim_buffer_2_sssp, prim_buffer_3_sssp, prim_buffer_3_spss, prim_buffer_2_spsp, prim_buffer_3_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_3_sdsp, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_4_spss, prim_buffer_3_spsp, prim_buffer_4_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_4_sdsp, prim_buffer_4_sssp, prim_buffer_5_sssp, prim_buffer_5_spss, prim_buffer_4_spsp, prim_buffer_5_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_1_sdsd, prim_buffer_1_sssd, prim_buffer_2_sssd, prim_buffer_2_spsp, prim_buffer_1_spsd, prim_buffer_2_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_2_sdsd, prim_buffer_2_sssd, prim_buffer_3_sssd, prim_buffer_3_spsp, prim_buffer_2_spsd, prim_buffer_3_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_3_sdsd, prim_buffer_3_sssd, prim_buffer_4_sssd, prim_buffer_4_spsp, prim_buffer_3_spsd, prim_buffer_4_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_4_sdsd, prim_buffer_4_sssd, prim_buffer_5_sssd, prim_buffer_5_spsp, prim_buffer_4_spsd, prim_buffer_5_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdps(prim_buffer_2_sdps, prim_buffer_2_ssps, prim_buffer_3_ssps, prim_buffer_3_spss, prim_buffer_2_spps, prim_buffer_3_spps, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdps(prim_buffer_3_sdps, prim_buffer_3_ssps, prim_buffer_4_ssps, prim_buffer_4_spss, prim_buffer_3_spps, prim_buffer_4_spps, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdps(prim_buffer_4_sdps, prim_buffer_4_ssps, prim_buffer_5_ssps, prim_buffer_5_spss, prim_buffer_4_spps, prim_buffer_5_spps, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpp(prim_buffer_1_sdpp, prim_buffer_1_sspp, prim_buffer_2_sspp, prim_buffer_2_spsp, prim_buffer_2_spps, prim_buffer_1_sppp, prim_buffer_2_sppp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpp(prim_buffer_2_sdpp, prim_buffer_2_sspp, prim_buffer_3_sspp, prim_buffer_3_spsp, prim_buffer_3_spps, prim_buffer_2_sppp, prim_buffer_3_sppp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpp(prim_buffer_3_sdpp, prim_buffer_3_sspp, prim_buffer_4_sspp, prim_buffer_4_spsp, prim_buffer_4_spps, prim_buffer_3_sppp, prim_buffer_4_sppp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpp(prim_buffer_4_sdpp, prim_buffer_4_sspp, prim_buffer_5_sspp, prim_buffer_5_spsp, prim_buffer_5_spps, prim_buffer_4_sppp, prim_buffer_5_sppp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpd(prim_buffer_0_sdpd, prim_buffer_0_sspd, prim_buffer_1_sspd, prim_buffer_1_spsd, prim_buffer_1_sppp, prim_buffer_0_sppd, prim_buffer_1_sppd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpd(prim_buffer_1_sdpd, prim_buffer_1_sspd, prim_buffer_2_sspd, prim_buffer_2_spsd, prim_buffer_2_sppp, prim_buffer_1_sppd, prim_buffer_2_sppd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpd(prim_buffer_2_sdpd, prim_buffer_2_sspd, prim_buffer_3_sspd, prim_buffer_3_spsd, prim_buffer_3_sppp, prim_buffer_2_sppd, prim_buffer_3_sppd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpd(prim_buffer_3_sdpd, prim_buffer_3_sspd, prim_buffer_4_sspd, prim_buffer_4_spsd, prim_buffer_4_sppp, prim_buffer_3_sppd, prim_buffer_4_sppd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdpd(prim_buffer_4_sdpd, prim_buffer_4_sspd, prim_buffer_5_sspd, prim_buffer_5_spsd, prim_buffer_5_sppp, prim_buffer_4_sppd, prim_buffer_5_sppd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pssd(prim_buffer_1_pssd, prim_buffer_2_sssp, prim_buffer_1_sssd, prim_buffer_2_sssd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pssd(prim_buffer_2_pssd, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pssd(prim_buffer_3_pssd, prim_buffer_4_sssp, prim_buffer_3_sssd, prim_buffer_4_sssd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pspp(prim_buffer_1_pspp, prim_buffer_2_sssp, prim_buffer_2_ssps, prim_buffer_1_sspp, prim_buffer_2_sspp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pspp(prim_buffer_2_pspp, prim_buffer_3_sssp, prim_buffer_3_ssps, prim_buffer_2_sspp, prim_buffer_3_sspp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pspp(prim_buffer_3_pspp, prim_buffer_4_sssp, prim_buffer_4_ssps, prim_buffer_3_sspp, prim_buffer_4_sspp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pspd(prim_buffer_0_pspd, prim_buffer_1_sssd, prim_buffer_1_sspp, prim_buffer_0_sspd, prim_buffer_1_sspd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pspd(prim_buffer_1_pspd, prim_buffer_2_sssd, prim_buffer_2_sspp, prim_buffer_1_sspd, prim_buffer_2_sspd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pspd(prim_buffer_2_pspd, prim_buffer_3_sssd, prim_buffer_3_sspp, prim_buffer_2_sspd, prim_buffer_3_sspd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pspd(prim_buffer_3_pspd, prim_buffer_4_sssd, prim_buffer_4_sspp, prim_buffer_3_sspd, prim_buffer_4_sspd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppsp(prim_buffer_2_ppsp, prim_buffer_2_sssp, prim_buffer_3_sssp, prim_buffer_3_spss, prim_buffer_2_spsp, prim_buffer_3_spsp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppsp(prim_buffer_3_ppsp, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_4_spss, prim_buffer_3_spsp, prim_buffer_4_spsp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppsd(prim_buffer_1_ppsd, prim_buffer_1_sssd, prim_buffer_2_sssd, prim_buffer_2_spsp, prim_buffer_1_spsd, prim_buffer_2_spsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppsd(prim_buffer_2_ppsd, prim_buffer_2_sssd, prim_buffer_3_sssd, prim_buffer_3_spsp, prim_buffer_2_spsd, prim_buffer_3_spsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppsd(prim_buffer_3_ppsd, prim_buffer_3_sssd, prim_buffer_4_sssd, prim_buffer_4_spsp, prim_buffer_3_spsd, prim_buffer_4_spsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppps(prim_buffer_2_ppps, prim_buffer_2_ssps, prim_buffer_3_ssps, prim_buffer_3_spss, prim_buffer_2_spps, prim_buffer_3_spps, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ppps(prim_buffer_3_ppps, prim_buffer_3_ssps, prim_buffer_4_ssps, prim_buffer_4_spss, prim_buffer_3_spps, prim_buffer_4_spps, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pppp(prim_buffer_1_pppp, prim_buffer_1_sspp, prim_buffer_2_sspp, prim_buffer_2_spsp, prim_buffer_2_spps, prim_buffer_1_sppp, prim_buffer_2_sppp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pppp(prim_buffer_2_pppp, prim_buffer_2_sspp, prim_buffer_3_sspp, prim_buffer_3_spsp, prim_buffer_3_spps, prim_buffer_2_sppp, prim_buffer_3_sppp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pppp(prim_buffer_3_pppp, prim_buffer_3_sspp, prim_buffer_4_sspp, prim_buffer_4_spsp, prim_buffer_4_spps, prim_buffer_3_sppp, prim_buffer_4_sppp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pppd(prim_buffer_0_pppd, prim_buffer_0_sspd, prim_buffer_1_sspd, prim_buffer_1_spsd, prim_buffer_1_sppp, prim_buffer_0_sppd, prim_buffer_1_sppd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pppd(prim_buffer_1_pppd, prim_buffer_1_sspd, prim_buffer_2_sspd, prim_buffer_2_spsd, prim_buffer_2_sppp, prim_buffer_1_sppd, prim_buffer_2_sppd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pppd(prim_buffer_2_pppd, prim_buffer_2_sspd, prim_buffer_3_sspd, prim_buffer_3_spsd, prim_buffer_3_sppp, prim_buffer_2_sppd, prim_buffer_3_sppd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pppd(prim_buffer_3_pppd, prim_buffer_3_sspd, prim_buffer_4_sspd, prim_buffer_4_spsd, prim_buffer_4_sppp, prim_buffer_3_sppd, prim_buffer_4_sppd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdss(prim_buffer_3_pdss, prim_buffer_3_spss, prim_buffer_4_spss, prim_buffer_3_sdss, prim_buffer_4_sdss, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdsp(prim_buffer_2_pdsp, prim_buffer_2_spsp, prim_buffer_3_spsp, prim_buffer_3_sdss, prim_buffer_2_sdsp, prim_buffer_3_sdsp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdsp(prim_buffer_3_pdsp, prim_buffer_3_spsp, prim_buffer_4_spsp, prim_buffer_4_sdss, prim_buffer_3_sdsp, prim_buffer_4_sdsp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdsd(prim_buffer_1_pdsd, prim_buffer_1_spsd, prim_buffer_2_spsd, prim_buffer_2_sdsp, prim_buffer_1_sdsd, prim_buffer_2_sdsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdsd(prim_buffer_2_pdsd, prim_buffer_2_spsd, prim_buffer_3_spsd, prim_buffer_3_sdsp, prim_buffer_2_sdsd, prim_buffer_3_sdsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdsd(prim_buffer_3_pdsd, prim_buffer_3_spsd, prim_buffer_4_spsd, prim_buffer_4_sdsp, prim_buffer_3_sdsd, prim_buffer_4_sdsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdps(prim_buffer_2_pdps, prim_buffer_2_spps, prim_buffer_3_spps, prim_buffer_3_sdss, prim_buffer_2_sdps, prim_buffer_3_sdps, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdps(prim_buffer_3_pdps, prim_buffer_3_spps, prim_buffer_4_spps, prim_buffer_4_sdss, prim_buffer_3_sdps, prim_buffer_4_sdps, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdpp(prim_buffer_1_pdpp, prim_buffer_1_sppp, prim_buffer_2_sppp, prim_buffer_2_sdsp, prim_buffer_2_sdps, prim_buffer_1_sdpp, prim_buffer_2_sdpp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdpp(prim_buffer_2_pdpp, prim_buffer_2_sppp, prim_buffer_3_sppp, prim_buffer_3_sdsp, prim_buffer_3_sdps, prim_buffer_2_sdpp, prim_buffer_3_sdpp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdpp(prim_buffer_3_pdpp, prim_buffer_3_sppp, prim_buffer_4_sppp, prim_buffer_4_sdsp, prim_buffer_4_sdps, prim_buffer_3_sdpp, prim_buffer_4_sdpp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdpd(prim_buffer_0_pdpd, prim_buffer_0_sppd, prim_buffer_1_sppd, prim_buffer_1_sdsd, prim_buffer_1_sdpp, prim_buffer_0_sdpd, prim_buffer_1_sdpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdpd(prim_buffer_1_pdpd, prim_buffer_1_sppd, prim_buffer_2_sppd, prim_buffer_2_sdsd, prim_buffer_2_sdpp, prim_buffer_1_sdpd, prim_buffer_2_sdpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdpd(prim_buffer_2_pdpd, prim_buffer_2_sppd, prim_buffer_3_sppd, prim_buffer_3_sdsd, prim_buffer_3_sdpp, prim_buffer_2_sdpd, prim_buffer_3_sdpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_pdpd(prim_buffer_3_pdpd, prim_buffer_3_sppd, prim_buffer_4_sppd, prim_buffer_4_sdsd, prim_buffer_4_sdpp, prim_buffer_3_sdpd, prim_buffer_4_sdpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dspd(prim_buffer_0_dspd, prim_buffer_0_sspd, prim_buffer_1_sspd, prim_buffer_1_pssd, prim_buffer_1_pspp, prim_buffer_0_pspd, prim_buffer_1_pspd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dspd(prim_buffer_1_dspd, prim_buffer_1_sspd, prim_buffer_2_sspd, prim_buffer_2_pssd, prim_buffer_2_pspp, prim_buffer_1_pspd, prim_buffer_2_pspd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dspd(prim_buffer_2_dspd, prim_buffer_2_sspd, prim_buffer_3_sspd, prim_buffer_3_pssd, prim_buffer_3_pspp, prim_buffer_2_pspd, prim_buffer_3_pspd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dpsd(prim_buffer_1_dpsd, prim_buffer_1_spsd, prim_buffer_2_spsd, prim_buffer_1_pssd, prim_buffer_2_pssd, prim_buffer_2_ppsp, prim_buffer_1_ppsd, prim_buffer_2_ppsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dpsd(prim_buffer_2_dpsd, prim_buffer_2_spsd, prim_buffer_3_spsd, prim_buffer_2_pssd, prim_buffer_3_pssd, prim_buffer_3_ppsp, prim_buffer_2_ppsd, prim_buffer_3_ppsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dppp(prim_buffer_1_dppp, prim_buffer_1_sppp, prim_buffer_2_sppp, prim_buffer_1_pspp, prim_buffer_2_pspp, prim_buffer_2_ppsp, prim_buffer_2_ppps, prim_buffer_1_pppp, prim_buffer_2_pppp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dppp(prim_buffer_2_dppp, prim_buffer_2_sppp, prim_buffer_3_sppp, prim_buffer_2_pspp, prim_buffer_3_pspp, prim_buffer_3_ppsp, prim_buffer_3_ppps, prim_buffer_2_pppp, prim_buffer_3_pppp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dppd(prim_buffer_0_dppd, prim_buffer_0_sppd, prim_buffer_1_sppd, prim_buffer_0_pspd, prim_buffer_1_pspd, prim_buffer_1_ppsd, prim_buffer_1_pppp, prim_buffer_0_pppd, prim_buffer_1_pppd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dppd(prim_buffer_1_dppd, prim_buffer_1_sppd, prim_buffer_2_sppd, prim_buffer_1_pspd, prim_buffer_2_pspd, prim_buffer_2_ppsd, prim_buffer_2_pppp, prim_buffer_1_pppd, prim_buffer_2_pppd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_dppd(prim_buffer_2_dppd, prim_buffer_2_sppd, prim_buffer_3_sppd, prim_buffer_2_pspd, prim_buffer_3_pspd, prim_buffer_3_ppsd, prim_buffer_3_pppp, prim_buffer_2_pppd, prim_buffer_3_pppd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ddsp(prim_buffer_2_ddsp, prim_buffer_2_sdsp, prim_buffer_3_sdsp, prim_buffer_2_ppsp, prim_buffer_3_ppsp, prim_buffer_3_pdss, prim_buffer_2_pdsp, prim_buffer_3_pdsp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ddsd(prim_buffer_1_ddsd, prim_buffer_1_sdsd, prim_buffer_2_sdsd, prim_buffer_1_ppsd, prim_buffer_2_ppsd, prim_buffer_2_pdsp, prim_buffer_1_pdsd, prim_buffer_2_pdsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ddsd(prim_buffer_2_ddsd, prim_buffer_2_sdsd, prim_buffer_3_sdsd, prim_buffer_2_ppsd, prim_buffer_3_ppsd, prim_buffer_3_pdsp, prim_buffer_2_pdsd, prim_buffer_3_pdsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ddps(prim_buffer_2_ddps, prim_buffer_2_sdps, prim_buffer_3_sdps, prim_buffer_2_ppps, prim_buffer_3_ppps, prim_buffer_3_pdss, prim_buffer_2_pdps, prim_buffer_3_pdps, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ddpp(prim_buffer_1_ddpp, prim_buffer_1_sdpp, prim_buffer_2_sdpp, prim_buffer_1_pppp, prim_buffer_2_pppp, prim_buffer_2_pdsp, prim_buffer_2_pdps, prim_buffer_1_pdpp, prim_buffer_2_pdpp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ddpp(prim_buffer_2_ddpp, prim_buffer_2_sdpp, prim_buffer_3_sdpp, prim_buffer_2_pppp, prim_buffer_3_pppp, prim_buffer_3_pdsp, prim_buffer_3_pdps, prim_buffer_2_pdpp, prim_buffer_3_pdpp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ddpd(prim_buffer_0_ddpd, prim_buffer_0_sdpd, prim_buffer_1_sdpd, prim_buffer_0_pppd, prim_buffer_1_pppd, prim_buffer_1_pdsd, prim_buffer_1_pdpp, prim_buffer_0_pdpd, prim_buffer_1_pdpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ddpd(prim_buffer_1_ddpd, prim_buffer_1_sdpd, prim_buffer_2_sdpd, prim_buffer_1_pppd, prim_buffer_2_pppd, prim_buffer_2_pdsd, prim_buffer_2_pdpp, prim_buffer_1_pdpd, prim_buffer_2_pdpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ddpd(prim_buffer_2_ddpd, prim_buffer_2_sdpd, prim_buffer_3_sdpd, prim_buffer_2_pppd, prim_buffer_3_pppd, prim_buffer_3_pdsd, prim_buffer_3_pdpp, prim_buffer_2_pdpd, prim_buffer_3_pdpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_fppd(prim_buffer_0_fppd, prim_buffer_0_pppd, prim_buffer_1_pppd, prim_buffer_0_dspd, prim_buffer_1_dspd, prim_buffer_1_dpsd, prim_buffer_1_dppp, prim_buffer_0_dppd, prim_buffer_1_dppd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_fppd(prim_buffer_1_fppd, prim_buffer_1_pppd, prim_buffer_2_pppd, prim_buffer_1_dspd, prim_buffer_2_dspd, prim_buffer_2_dpsd, prim_buffer_2_dppp, prim_buffer_1_dppd, prim_buffer_2_dppd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_fdsd(prim_buffer_1_fdsd, prim_buffer_1_pdsd, prim_buffer_2_pdsd, prim_buffer_1_dpsd, prim_buffer_2_dpsd, prim_buffer_2_ddsp, prim_buffer_1_ddsd, prim_buffer_2_ddsd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_fdpp(prim_buffer_1_fdpp, prim_buffer_1_pdpp, prim_buffer_2_pdpp, prim_buffer_1_dppp, prim_buffer_2_dppp, prim_buffer_2_ddsp, prim_buffer_2_ddps, prim_buffer_1_ddpp, prim_buffer_2_ddpp, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_fdpd(prim_buffer_0_fdpd, prim_buffer_0_pdpd, prim_buffer_1_pdpd, prim_buffer_0_dppd, prim_buffer_1_dppd, prim_buffer_1_ddsd, prim_buffer_1_ddpp, prim_buffer_0_ddpd, prim_buffer_1_ddpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_fdpd(prim_buffer_1_fdpd, prim_buffer_1_pdpd, prim_buffer_2_pdpd, prim_buffer_1_dppd, prim_buffer_2_dppd, prim_buffer_2_ddsd, prim_buffer_2_ddpp, prim_buffer_1_ddpd, prim_buffer_2_ddpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_gdpd(prim_buffer_0_gdpd, prim_buffer_0_ddpd, prim_buffer_1_ddpd, prim_buffer_0_fppd, prim_buffer_1_fppd, prim_buffer_1_fdsd, prim_buffer_1_fdpp, prim_buffer_0_fdpd, prim_buffer_1_fdpd, pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            t4c_geom::comp_geom2000_ddpd_0(prim_buffer_0_geom2000_ddpd, prim_buffer_0_sdpd, prim_buffer_0_ddpd, prim_buffer_0_gdpd, a_exp );

            t2cfunc::reduce(cart_buffer_0_geom2000_ddpd, prim_buffer_0_geom2000_ddpd, ket_dim, ket_npgtos);
        }

        t4cfunc::full_transform<2, 2, 1, 2>(spher_buffer_0_geom2000_ddpd, cart_spher_buffer_0_geom2000_ddpd);

        distributor->distribute(spher_buffer_0_geom2000_ddpd, a_indices, b_indices, c_indices, d_indices, 2, 2, 1, 2, 2, 0, 0, 0, i, ket_indices);
    }
}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecDDPD_hpp */
