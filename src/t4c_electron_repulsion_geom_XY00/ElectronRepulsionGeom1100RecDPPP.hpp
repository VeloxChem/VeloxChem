#ifndef ElectronRepulsionGeom1100RecDPPP_hpp
#define ElectronRepulsionGeom1100RecDPPP_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom0100ContrRecSPXX.hpp"
#include "ElectronRepulsionGeom1000ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecDPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecPPXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1100ContrRecSPXX.hpp"
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
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dB^(1)(DP|1/|r-r'||PP)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1100_dppp(T& distributor,
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

    CSimdArray<double> pbuffer(1210, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(828, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(1989, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(7902, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1215, 1);

    // setup Boys fuction data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

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

                t4cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 8);

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 6, pfactors, 16, bf_data, 6);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 7, pfactors, 16, bf_data, 7);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 8, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 29, 0, 1, 8, 11, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 35, 1, 2, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 41, 2, 3, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 47, 3, 4, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 4, 5, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 5, 6, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 65, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 68, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 71, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 74, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 77, 1, 8, 11, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 86, 2, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 95, 3, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 104, 4, 17, 20, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 113, 5, 20, 23, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 122, 11, 29, 35, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 140, 14, 35, 41, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 158, 17, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 176, 20, 47, 53, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 194, 23, 53, 59, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 212, 1, 2, 65, 68, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 218, 2, 3, 68, 71, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 224, 3, 4, 71, 74, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 230, 8, 11, 65, 77, 86, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 248, 11, 14, 68, 86, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 266, 14, 17, 71, 95, 104, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 284, 17, 20, 74, 104, 113, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 302, 29, 35, 86, 122, 140, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 338, 35, 41, 95, 140, 158, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 374, 41, 47, 104, 158, 176, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 410, 47, 53, 113, 176, 194, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 446, 65, 68, 212, 218, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 456, 68, 71, 218, 224, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 466, 77, 86, 212, 230, 248, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 496, 86, 95, 218, 248, 266, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 526, 95, 104, 224, 266, 284, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 556, 122, 140, 248, 302, 338, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 616, 140, 158, 266, 338, 374, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 676, 158, 176, 284, 374, 410, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 736, 212, 218, 446, 456, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 751, 230, 248, 446, 466, 496, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 796, 248, 266, 456, 496, 526, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 841, 302, 338, 496, 556, 616, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 931, 338, 374, 526, 616, 676, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1021, 466, 496, 736, 751, 796, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 1084, 556, 616, 796, 841, 931, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 8, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 3, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 9, pbuffer, 77, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 122, 18, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {230, 248});

                pbuffer.scale(2.0 * b_exp, {302, 338});

                pbuffer.scale(2.0 * b_exp, {466, 496});

                pbuffer.scale(2.0 * b_exp, {556, 616});

                t2cfunc::reduce(cbuffer, 36, pbuffer, 230, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 54, pbuffer, 302, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 90, pbuffer, 466, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 120, pbuffer, 556, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {8, 11});

                pbuffer.scale(2.0 * a_exp, {29, 35});

                pbuffer.scale(2.0 * a_exp, {77, 86});

                pbuffer.scale(2.0 * a_exp, {122, 140});

                pbuffer.scale(a_exp / b_exp, {230, 248});

                pbuffer.scale(a_exp / b_exp, {302, 338});

                pbuffer.scale(a_exp / b_exp, {466, 496});

                pbuffer.scale(a_exp / b_exp, {556, 616});

                t2cfunc::reduce(cbuffer, 180, pbuffer, 8, 3, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 183, pbuffer, 29, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 189, pbuffer, 77, 9, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 198, pbuffer, 122, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 216, pbuffer, 230, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 234, pbuffer, 302, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 270, pbuffer, 466, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 300, pbuffer, 556, 60, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * b_exp, {230, 248});

                pbuffer.scale(2.0 * b_exp, {302, 338});

                pbuffer.scale(2.0 * b_exp, {466, 496});

                pbuffer.scale(2.0 * b_exp, {556, 616});

                pbuffer.scale(4.0 * a_exp * b_exp, {751, 796});

                pbuffer.scale(4.0 * a_exp * b_exp, {841, 931});

                pbuffer.scale(4.0 * a_exp * b_exp, {1021, 1084});

                pbuffer.scale(4.0 * a_exp * b_exp, {1084, 1210});

                t2cfunc::reduce(cbuffer, 360, pbuffer, 230, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 378, pbuffer, 302, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 414, pbuffer, 466, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 444, pbuffer, 556, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 504, pbuffer, 751, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 549, pbuffer, 841, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 639, pbuffer, 1021, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 702, pbuffer, 1084, 126, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 0, cbuffer, 0, 3, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 9, cbuffer, 9, 18, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 279, cbuffer, 36, 54, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 333, cbuffer, 90, 120, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 423, cbuffer, 180, 183, cfactors, 6, 0, 0);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 432, cbuffer, 189, 198, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 540, cbuffer, 216, 234, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 756, cbuffer, 270, 300, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1521, cbuffer, 360, 378, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1575, cbuffer, 414, 444, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1665, cbuffer, 504, 549, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxpp(ckbuffer, 1800, cbuffer, 639, 702, cfactors, 6, 0, 5);

            t4cfunc::ket_transform<1, 1>(skbuffer, 0, ckbuffer, 0, 0, 0);

            t4cfunc::ket_transform<1, 1>(skbuffer, 9, ckbuffer, 9, 0, 1);

            t4cfunc::ket_transform<1, 1>(skbuffer, 5949, ckbuffer, 279, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 6003, ckbuffer, 333, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 6093, ckbuffer, 423, 0, 0);

            t4cfunc::ket_transform<1, 1>(skbuffer, 6102, ckbuffer, 432, 0, 1);

            t4cfunc::ket_transform<1, 1>(skbuffer, 6210, ckbuffer, 540, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 6426, ckbuffer, 756, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 7434, ckbuffer, 1521, 0, 2);

            t4cfunc::ket_transform<1, 1>(skbuffer, 7488, ckbuffer, 1575, 0, 3);

            t4cfunc::ket_transform<1, 1>(skbuffer, 7578, ckbuffer, 1665, 0, 4);

            t4cfunc::ket_transform<1, 1>(skbuffer, 7713, ckbuffer, 1800, 0, 5);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 7191, 6102, 6210, r_ab, 1, 1);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 7272, 6210, 6426, r_ab, 1, 1);

            erirec::comp_bra_geom10_hrr_electron_repulsion_ppxx(skbuffer, 2061, 9, 7191, 7272, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 36, 0, 5949, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 360, 9, 6003, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_ppxx(skbuffer, 1818, 9, 36, 360, r_ab, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_spxx(skbuffer, 6129, 6093, 7434, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sdxx(skbuffer, 6264, 6102, 7488, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sfxx(skbuffer, 6516, 6210, 7578, 1, 1);

            erirec::comp_bra_geom01_hrr_electron_repulsion_sgxx(skbuffer, 6786, 6426, 7713, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_spxx(skbuffer, 117, 6102, 6129, 6264, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sdxx(skbuffer, 522, 6210, 6264, 6516, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_sfxx(skbuffer, 1008, 6426, 6516, 6786, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_ppxx(skbuffer, 2304, 36, 7191, 117, 522, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_pdxx(skbuffer, 3033, 360, 7272, 522, 1008, r_ab, 1, 1);

            erirec::comp_bra_geom11_hrr_electron_repulsion_dpxx(skbuffer, 4491, 1818, 2061, 2304, 3033, r_ab, 1, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 0, skbuffer, 4491, 1, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 135, skbuffer, 4653, 1, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 270, skbuffer, 4815, 1, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 405, skbuffer, 4977, 1, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 540, skbuffer, 5139, 1, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 675, skbuffer, 5301, 1, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 810, skbuffer, 5463, 1, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 945, skbuffer, 5625, 1, 1);

            t4cfunc::bra_transform<2, 1>(sbuffer, 1080, skbuffer, 5787, 1, 1);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 1, 1, 1, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1100RecDPPP_hpp */