#ifndef ElectronRepulsionGeom1010RecFDDS_hpp
#define ElectronRepulsionGeom1010RecFDDS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
#include "ElectronRepulsionGeom1010ContrRecDDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSDXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FD|1/|r-r'||DS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_fdds(T& distributor,
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

    CSimdArray<double> pbuffer(4215, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2664, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(9657, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(16470, 1);

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

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 155, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 158, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 161, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 164, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 167, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 170, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 173, 1, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 182, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 191, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 200, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 209, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 218, 6, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 227, 13, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 245, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 263, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 281, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 299, 25, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 317, 28, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 335, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 365, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 395, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 425, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 455, 67, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 485, 73, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 515, 0, 1, 155, 158, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 521, 1, 2, 158, 161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 527, 2, 3, 161, 164, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 533, 3, 4, 164, 167, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 539, 4, 5, 167, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 545, 10, 13, 158, 173, 182, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 563, 13, 16, 161, 182, 191, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 581, 16, 19, 164, 191, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 599, 19, 22, 167, 200, 209, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 617, 22, 25, 170, 209, 218, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 635, 37, 43, 182, 227, 245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 671, 43, 49, 191, 245, 263, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 707, 49, 55, 200, 263, 281, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 743, 55, 61, 209, 281, 299, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 779, 61, 67, 218, 299, 317, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 815, 85, 95, 245, 335, 365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 875, 95, 105, 263, 365, 395, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 935, 105, 115, 281, 395, 425, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 995, 115, 125, 299, 425, 455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1055, 125, 135, 317, 455, 485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1115, 155, 158, 515, 521, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1125, 158, 161, 521, 527, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1135, 161, 164, 527, 533, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1145, 164, 167, 533, 539, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1155, 173, 182, 521, 545, 563, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1185, 182, 191, 527, 563, 581, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1215, 191, 200, 533, 581, 599, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1245, 200, 209, 539, 599, 617, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1275, 227, 245, 563, 635, 671, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1335, 245, 263, 581, 671, 707, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1395, 263, 281, 599, 707, 743, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1455, 281, 299, 617, 743, 779, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1515, 335, 365, 671, 815, 875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1615, 365, 395, 707, 875, 935, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1715, 395, 425, 743, 935, 995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1815, 425, 455, 779, 995, 1055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1915, 515, 521, 1115, 1125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1930, 521, 527, 1125, 1135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1945, 527, 533, 1135, 1145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1960, 545, 563, 1125, 1155, 1185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2005, 563, 581, 1135, 1185, 1215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2050, 581, 599, 1145, 1215, 1245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2095, 635, 671, 1185, 1275, 1335, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2185, 671, 707, 1215, 1335, 1395, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2275, 707, 743, 1245, 1395, 1455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2365, 815, 875, 1335, 1515, 1615, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2515, 875, 935, 1395, 1615, 1715, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2665, 935, 995, 1455, 1715, 1815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 2815, 1115, 1125, 1915, 1930, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 2836, 1125, 1135, 1930, 1945, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2857, 1155, 1185, 1930, 1960, 2005, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2920, 1185, 1215, 1945, 2005, 2050, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2983, 1275, 1335, 2005, 2095, 2185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3109, 1335, 1395, 2050, 2185, 2275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3235, 1515, 1615, 2185, 2365, 2515, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3445, 1615, 1715, 2275, 2515, 2665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 3655, 1915, 1930, 2815, 2836, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 3683, 1960, 2005, 2836, 2857, 2920, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 3767, 2095, 2185, 2920, 2983, 3109, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 3935, 2365, 2515, 3109, 3235, 3445, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 515, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 6, pbuffer, 545, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 24, pbuffer, 1115, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 34, pbuffer, 1155, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 64, pbuffer, 1915, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 79, pbuffer, 1960, 45, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {515, 521});

                pbuffer.scale(2.0 * a_exp, {545, 563});

                pbuffer.scale(2.0 * a_exp, {1115, 1125});

                pbuffer.scale(2.0 * a_exp, {1155, 1185});

                pbuffer.scale(2.0 * a_exp, {1915, 1930});

                pbuffer.scale(2.0 * a_exp, {1960, 2005});

                pbuffer.scale(2.0 * a_exp, {2815, 2836});

                pbuffer.scale(2.0 * a_exp, {2857, 2920});

                pbuffer.scale(2.0 * a_exp, {3655, 3683});

                pbuffer.scale(2.0 * a_exp, {3683, 3767});

                t2cfunc::reduce(cbuffer, 744, pbuffer, 515, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 750, pbuffer, 545, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 768, pbuffer, 1115, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 778, pbuffer, 1155, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 808, pbuffer, 1915, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 823, pbuffer, 1960, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 868, pbuffer, 2815, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 889, pbuffer, 2857, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 952, pbuffer, 3655, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 980, pbuffer, 3683, 84, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {515, 521});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {545, 563});

                pbuffer.scale(pfactors, 0, 2.0, {635, 671});

                pbuffer.scale(pfactors, 0, 2.0, {815, 875});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1115, 1125});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1155, 1185});

                pbuffer.scale(pfactors, 0, 2.0, {1275, 1335});

                pbuffer.scale(pfactors, 0, 2.0, {1515, 1615});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1915, 1930});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1960, 2005});

                pbuffer.scale(pfactors, 0, 2.0, {2095, 2185});

                pbuffer.scale(pfactors, 0, 2.0, {2365, 2515});

                t2cfunc::reduce(cbuffer, 124, pbuffer, 515, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 130, pbuffer, 545, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 148, pbuffer, 635, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 184, pbuffer, 815, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 244, pbuffer, 1115, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 254, pbuffer, 1155, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 284, pbuffer, 1275, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 344, pbuffer, 1515, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 444, pbuffer, 1915, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 459, pbuffer, 1960, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 504, pbuffer, 2095, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 594, pbuffer, 2365, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {515, 521});

                pbuffer.scale(2.0 * a_exp, {545, 563});

                pbuffer.scale(2.0 * a_exp, {635, 671});

                pbuffer.scale(2.0 * a_exp, {815, 875});

                pbuffer.scale(2.0 * a_exp, {1115, 1125});

                pbuffer.scale(2.0 * a_exp, {1155, 1185});

                pbuffer.scale(2.0 * a_exp, {1275, 1335});

                pbuffer.scale(2.0 * a_exp, {1515, 1615});

                pbuffer.scale(2.0 * a_exp, {1915, 1930});

                pbuffer.scale(2.0 * a_exp, {1960, 2005});

                pbuffer.scale(2.0 * a_exp, {2095, 2185});

                pbuffer.scale(2.0 * a_exp, {2365, 2515});

                pbuffer.scale(pfactors, 0, 2.0, {2815, 2836});

                pbuffer.scale(pfactors, 0, 2.0, {2857, 2920});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2983, 3109});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3235, 3445});

                pbuffer.scale(pfactors, 0, 2.0, {3655, 3683});

                pbuffer.scale(pfactors, 0, 2.0, {3683, 3767});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3767, 3935});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3935, 4215});

                t2cfunc::reduce(cbuffer, 1064, pbuffer, 515, 6, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1070, pbuffer, 545, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1088, pbuffer, 635, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1124, pbuffer, 815, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1184, pbuffer, 1115, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1194, pbuffer, 1155, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1224, pbuffer, 1275, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1284, pbuffer, 1515, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1384, pbuffer, 1915, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1399, pbuffer, 1960, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1444, pbuffer, 2095, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1534, pbuffer, 2365, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1684, pbuffer, 2815, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1705, pbuffer, 2857, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1768, pbuffer, 2983, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1894, pbuffer, 3235, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2104, pbuffer, 3655, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2132, pbuffer, 3683, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2216, pbuffer, 3767, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2384, pbuffer, 3935, 280, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 180, cbuffer, 0, 6, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 822, cbuffer, 24, 34, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1842, cbuffer, 64, 79, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2877, cbuffer, 744, 750, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 3519, cbuffer, 768, 778, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 4539, cbuffer, 808, 823, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 6024, cbuffer, 868, 889, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 8061, cbuffer, 952, 980, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 124, 130, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 18, cbuffer, 130, 148, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 72, cbuffer, 148, 184, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 198, cbuffer, 0, 0, 18, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 252, cbuffer, 6, 18, 72, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 414, 180, 198, 252, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 522, cbuffer, 244, 254, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 552, cbuffer, 254, 284, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 642, cbuffer, 284, 344, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 852, cbuffer, 24, 522, 552, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 942, cbuffer, 34, 552, 642, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1212, 822, 852, 942, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1392, cbuffer, 444, 459, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1437, cbuffer, 459, 504, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1572, cbuffer, 504, 594, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1887, cbuffer, 64, 1392, 1437, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2022, cbuffer, 79, 1437, 1572, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 2427, 1842, 1887, 2022, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2697, cbuffer, 1064, 1070, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2715, cbuffer, 1070, 1088, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2769, cbuffer, 1088, 1124, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2895, cbuffer, 744, 2697, 2715, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2949, cbuffer, 750, 2715, 2769, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 3111, 2877, 2895, 2949, cfactors, 6, 0, 2);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 3219, cbuffer, 1184, 1194, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3249, cbuffer, 1194, 1224, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3339, cbuffer, 1224, 1284, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 3549, cbuffer, 768, 3219, 3249, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3639, cbuffer, 778, 3249, 3339, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 3909, 3519, 3549, 3639, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 4089, cbuffer, 1384, 1399, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4134, cbuffer, 1399, 1444, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4269, cbuffer, 1444, 1534, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 4584, cbuffer, 808, 4089, 4134, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 4719, cbuffer, 823, 4134, 4269, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 5124, 4539, 4584, 4719, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 5394, cbuffer, 1684, 1705, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 5457, cbuffer, 1705, 1768, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 5646, cbuffer, 1768, 1894, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 6087, cbuffer, 868, 5394, 5457, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 6276, cbuffer, 889, 5457, 5646, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 6843, 6024, 6087, 6276, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 7221, cbuffer, 2104, 2132, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 7305, cbuffer, 2132, 2216, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 7557, cbuffer, 2216, 2384, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 8145, cbuffer, 952, 7221, 7305, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 8397, cbuffer, 980, 7305, 7557, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 9153, 8061, 8145, 8397, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<2, 0>(skbuffer, 0, ckbuffer, 414, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 30, ckbuffer, 450, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 60, ckbuffer, 486, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 360, ckbuffer, 1212, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 410, ckbuffer, 1272, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 460, ckbuffer, 1332, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 960, ckbuffer, 2427, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1035, ckbuffer, 2517, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 1110, ckbuffer, 2607, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15270, ckbuffer, 3111, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15300, ckbuffer, 3147, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15330, ckbuffer, 3183, 0, 2);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15360, ckbuffer, 3909, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15410, ckbuffer, 3969, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15460, ckbuffer, 4029, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15510, ckbuffer, 5124, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15585, ckbuffer, 5214, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15660, ckbuffer, 5304, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15735, ckbuffer, 6843, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15840, ckbuffer, 6969, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 15945, ckbuffer, 7095, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 16050, ckbuffer, 9153, 0, 6);

            t4cfunc::ket_transform<2, 0>(skbuffer, 16190, ckbuffer, 9321, 0, 6);

            t4cfunc::ket_transform<2, 0>(skbuffer, 16330, ckbuffer, 9489, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2805, 0, 360, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2895, 30, 410, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 2985, 60, 460, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 3885, 360, 960, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 4035, 410, 1035, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 4185, 460, 1110, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 7710, 2805, 3885, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 7890, 2895, 4035, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_ddxx(skbuffer, 8070, 2985, 4185, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sdxx(skbuffer, 90, 15270, 15360, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 510, 15360, 15510, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 1185, 15510, 15735, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 1860, 15735, 16050, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pdxx(skbuffer, 3075, 0, 90, 510, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 4335, 360, 510, 1185, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 5685, 960, 1185, 1860, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ddxx(skbuffer, 8250, 2805, 3075, 4335, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 9870, 3885, 4335, 5685, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_fdxx(skbuffer, 12570, 7710, 8250, 9870, r_ab, 2, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 0, skbuffer, 12570, 2, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 175, skbuffer, 12870, 2, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 350, skbuffer, 13170, 2, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 525, skbuffer, 13470, 2, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 700, skbuffer, 13770, 2, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 875, skbuffer, 14070, 2, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1050, skbuffer, 14370, 2, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1225, skbuffer, 14670, 2, 0);

            t4cfunc::bra_transform<3, 2>(sbuffer, 1400, skbuffer, 14970, 2, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 2, 2, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFDDS_hpp */
