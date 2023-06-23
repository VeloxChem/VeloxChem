#include "OverlapRecGF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "MathConst.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGF(CSubMatrix*      matrix,
              const CGtoBlock& bra_gto_block,
              const CGtoBlock& ket_gto_block,
              const bool       ang_order,
              const int64_t    bra_first,
              const int64_t    bra_last) -> void

{
    // spherical transformation factors

    const double f4_35 = 4.0 * std::sqrt(35);

    const double f4_17 = 4.0 * std::sqrt(17.5);

    const double f4_5 = 4.0 * std::sqrt(5.0);

    const double f4_2 = 4.0 * std::sqrt(2.5);

    const double f3_5 = std::sqrt(2.5);

    const double f3_15 = 2.0 * std::sqrt(15.0);

    const double f3_3 = std::sqrt(1.5);

    // intialize GTOs data on bra side

    const auto bra_gto_coords = bra_gto_block.getCoordinates();

    const auto bra_gto_exps = bra_gto_block.getExponents();

    const auto bra_gto_norms = bra_gto_block.getNormalizationFactors();

    const auto bra_gto_indexes = bra_gto_block.getOrbitalIndexes();

    const auto bra_ncgtos = bra_gto_block.getNumberOfBasisFunctions();

    const auto bra_npgtos = bra_gto_block.getNumberOfPrimitives();

    // intialize GTOs data on ket side

    const auto ket_gto_coords = ket_gto_block.getCoordinates();

    const auto ket_gto_exps = ket_gto_block.getExponents();

    const auto ket_gto_norms = ket_gto_block.getNormalizationFactors();

    const auto ket_gto_indexes = ket_gto_block.getOrbitalIndexes();

    const auto ket_ncgtos = ket_gto_block.getNumberOfBasisFunctions();

    const auto ket_npgtos = ket_gto_block.getNumberOfPrimitives();

    // initialize aligned arrays for ket side

    alignas(64) TDoubleArray ket_coords_x;

    alignas(64) TDoubleArray ket_coords_y;

    alignas(64) TDoubleArray ket_coords_z;

    alignas(64) TDoubleArray ket_exps;

    alignas(64) TDoubleArray ket_norms;

    // initialize contracted integrals buffer

    alignas(64) TDoubleArray buffer_xxx;

    alignas(64) TDoubleArray buffer_xxy;

    alignas(64) TDoubleArray buffer_xxz;

    alignas(64) TDoubleArray buffer_xyy;

    alignas(64) TDoubleArray buffer_xyz;

    alignas(64) TDoubleArray buffer_xzz;

    alignas(64) TDoubleArray buffer_yyy;

    alignas(64) TDoubleArray buffer_yyz;

    alignas(64) TDoubleArray buffer_yzz;

    alignas(64) TDoubleArray buffer_zzz;

    // loop over integral batches

    const auto nbatches = batch::getNumberOfBatches(ket_ncgtos, simd_width);

    for (int64_t i = 0; i < nbatches; i++)
    {
        const auto [ket_first, ket_last] = batch::getBatchRange(i, ket_ncgtos, simd_width);

        const auto ket_dim = ket_last - ket_first;

        simd::loadCoordinates(ket_coords_x, ket_coords_y, ket_coords_z, ket_gto_coords, ket_first, ket_last);

        for (int64_t j = bra_first; j < bra_last; j++)
        {
            const auto bra_coord = bra_gto_coords[j];

            // compute primitive integrals block (XXXX)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_XXXX_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXY)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_XXXY_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXXZ)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_XXXZ_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYY)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_XXYY_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 6.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 6.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -1.50 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -6.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXYZ)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_XXYZ_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZZ)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_XXZZ_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYY)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_XYYY_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYYZ)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_XYYZ_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 7, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZZ)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_XYZZ_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 6.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZZ)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_XZZZ_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYY)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_YYYY_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYYZ)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_YYYZ_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZZ)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_YYZZ_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZZ)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_YZZZ_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZZ)

            simd::zero(buffer_xxx);

            simd::zero(buffer_xxy);

            simd::zero(buffer_xxz);

            simd::zero(buffer_xyy);

            simd::zero(buffer_xyz);

            simd::zero(buffer_xzz);

            simd::zero(buffer_yyy);

            simd::zero(buffer_yyz);

            simd::zero(buffer_yzz);

            simd::zero(buffer_zzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGF_ZZZZ_T(buffer_xxx,
                                                          buffer_xxy,
                                                          buffer_xxz,
                                                          buffer_xyy,
                                                          buffer_xyz,
                                                          buffer_xzz,
                                                          buffer_yyy,
                                                          buffer_yyz,
                                                          buffer_yzz,
                                                          buffer_zzz,
                                                          bra_exp,
                                                          bra_norm,
                                                          bra_coord,
                                                          ket_exps,
                                                          ket_norms,
                                                          ket_coords_x,
                                                          ket_coords_y,
                                                          ket_coords_z,
                                                          ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_xxx, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 8.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 8.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -8.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);
        }
    }
}

auto
compPrimitiveOverlapGF_XXXX_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x + 6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x +
                               (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x + 9.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (27.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x);

        fints_xxx[i] += fss * (3.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x +
                               (15.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_x * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x + 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y + 6.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y);

        fints_xxy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x + 4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z + 6.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z);

        fints_xxz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x + 2.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x + 3.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x);

        fints_xyy[i] += fss * (fe_0 * fe_0 * rpa_x * rpa_x * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x + 2.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y +
                               3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x +
                               rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x + 2.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_x + 3.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x);

        fints_xzz[i] += fss * (fe_0 * fe_0 * rpa_x * rpa_x * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_yyy[i] += fss * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y;

        fints_yyz[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z);

        fints_yyz[i] += fss * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y;

        fints_yzz[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_yzz[i] += fss * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y;

        fints_zzz[i] += fss * (3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z);

        fints_zzz[i] += fss * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z;
    }
}

auto
compPrimitiveOverlapGF_XXXY_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x + (9.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x + (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x);

        fints_xxx[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_y * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x + 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y);

        fints_xxy[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x + 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_xxz[i] += fss * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x;

        fints_xyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x + fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x);

        fints_xyy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xyy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x);

        fints_xyz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x);

        fints_xzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y);

        fints_yyy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_y * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z + fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y);

        fints_yyz[i] += fss * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y;

        fints_yzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z);

        fints_yzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z + rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveOverlapGF_XXXZ_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x + (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x + (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x);

        fints_xxx[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x + 3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x);

        fints_xxy[i] += fss * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x;

        fints_xxz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x + 3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z);

        fints_xxz[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x);

        fints_xyy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x);

        fints_xyz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_x + fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x);

        fints_xzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_xzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y + rpa_z * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y);

        fints_yyz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_y + fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y);

        fints_yzz[i] += fss * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y;

        fints_zzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z +
                   (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z);

        fints_zzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveOverlapGF_XXYY_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x);

        fints_xxx[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxx[i] += fss * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_y * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] +=
            fss * (fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x + 2.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxy[i] += fss * (2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                               fe_0 * fe_0 * rpa_x * rpb_y * rpb_x);

        fints_xxy[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (2.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z);

        fints_xxz[i] += fss * (fe_0 * fe_0 * rpa_x * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] +=
            fss * (2.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x + fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyy[i] += fss * (2.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y + fe_0 * fe_0 * rpa_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y);

        fints_xyy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x + fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x + fe_0 * fe_0 * rpa_y * rpa_x * rpb_z);

        fints_xyz[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x + rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x);

        fints_xzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xzz[i] += fss * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * (3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x);

        fints_yyy[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyy[i] += fss * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (2.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y + fe_0 * fe_0 * rpa_y * rpb_z * rpb_y);

        fints_yyz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x);

        fints_yzz[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yzz[i] += fss * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z);

        fints_zzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveOverlapGF_XXYZ_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] += fss * (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x);

        fints_xxx[i] += fss * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x;

        fints_xxy[i] += fss * (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_xxy[i] += fss * (fe_0 * fe_0 * rpa_z * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_xxz[i] += fss * (fe_0 * fe_0 * rpa_y * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x + fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x);

        fints_xyy[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + fe_0 * fe_0 * rpa_z * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z);

        fints_xyz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xyz[i] += fss * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x + fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x);

        fints_xzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x + rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x);

        fints_yyy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y + fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_yyz[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] +=
            fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z + fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y +
                   (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_yzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_yzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x);

        fints_zzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveOverlapGF_XXZZ_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x);

        fints_xxx[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxx[i] += fss * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y);

        fints_xxy[i] += fss * (fe_0 * fe_0 * rpa_x * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] +=
            fss * (fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x + 2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               fe_0 * fe_0 * rpa_x * rpb_z * rpb_x);

        fints_xxz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x);

        fints_xyy[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xyy[i] += fss * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x + fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x + fe_0 * fe_0 * rpa_z * rpa_x * rpb_y);

        fints_xyz[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x + rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] +=
            fss * (2.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x + fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_xzz[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z + fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z);

        fints_xzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y);

        fints_yyy[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x);

        fints_yyz[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_yyz[i] += fss * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (2.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y + fe_0 * fe_0 * rpa_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * (3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x);

        fints_zzz[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_zzz[i] += fss * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveOverlapGF_XYYY_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_x +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_xxx[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_y * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y);

        fints_xxy[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_xxz[i] += fss * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x;

        fints_xyy[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x + 3.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x);

        fints_xyy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xyy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z);

        fints_xyz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_xzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y + (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x);

        fints_yyy[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_y * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y);

        fints_yyz[i] += fss * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y;

        fints_yzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x);

        fints_yzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z + rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveOverlapGF_XYYZ_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x);

        fints_xxx[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y +
                               fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x);

        fints_xxy[i] += fss * (fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z + fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z);

        fints_xxz[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x + fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_xyy[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_xyz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xyz[i] += fss * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] +=
            fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x + fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x +
                   (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y);

        fints_xzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_xzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y);

        fints_yyy[i] += fss * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y;

        fints_yyz[i] += fss * (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z);

        fints_yyz[i] += fss * (fe_0 * fe_0 * rpa_y * rpa_x * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y + fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x);

        fints_yzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y + fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y + rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x);

        fints_zzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveOverlapGF_XYZZ_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x);

        fints_xxx[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y + fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x);

        fints_xxy[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxy[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z +
                               fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x);

        fints_xxz[i] += fss * (fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x + rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y + fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y);

        fints_xyy[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xyy[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_xyz[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_xyz[i] += fss * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x + fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_xzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y);

        fints_yyy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z +
                               fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x);

        fints_yyz[i] += fss * (fe_0 * fe_0 * rpa_z * rpa_x * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y + rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y + fe_0 * fe_0 * rpa_z * rpa_x * rpb_z);

        fints_yzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z);

        fints_zzz[i] += fss * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z;
    }
}

auto
compPrimitiveOverlapGF_XZZZ_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_x +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x);

        fints_xxx[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y + fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x);

        fints_xxy[i] += fss * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x;

        fints_xxz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z + fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z);

        fints_xxz[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y);

        fints_xyy[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y);

        fints_xyz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x + 3.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x);

        fints_xzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_xzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y + rpa_z * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x);

        fints_yyz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y + 3.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y);

        fints_yzz[i] += fss * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y;

        fints_zzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z + (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x);

        fints_zzz[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveOverlapGF_YYYY_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpb_x +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_xxx[i] += fss * rpa_y * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x;

        fints_xxy[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x + 2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpb_y + 3.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y);

        fints_xxy[i] += fss * (fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_y * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpb_z +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z);

        fints_xxz[i] += fss * rpa_y * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x;

        fints_xyy[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x + 4.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpb_x + 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x);

        fints_xyy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_y * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x + 2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x +
                               3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x +
                               rpa_y * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_xzz[i] += fss * rpa_y * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x;

        fints_yyy[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y + 6.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y +
                               (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpb_y + 9.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y +
                               (27.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y);

        fints_yyy[i] += fss * (3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y +
                               (15.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_y * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y + 4.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpb_z + 6.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z);

        fints_yyz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y + 2.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpb_y + 3.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y);

        fints_yzz[i] += fss * (fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_y * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * (3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpb_z +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z);

        fints_zzz[i] += fss * rpa_y * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z;
    }
}

auto
compPrimitiveOverlapGF_YYYZ_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + rpa_z * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y);

        fints_xxy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x + 3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x);

        fints_xyy[i] += fss * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x;

        fints_xyz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_x + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_xzz[i] += fss * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x;

        fints_yyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y + (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y);

        fints_yyy[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_yyz[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_y + fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_yzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_yzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_zzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveOverlapGF_YYZZ_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x);

        fints_xxx[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y);

        fints_xxy[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xxy[i] += fss * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y);

        fints_xxz[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_xxz[i] += fss * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x);

        fints_xyy[i] += fss * (fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x + fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x + fe_0 * fe_0 * rpa_z * rpa_y * rpb_x);

        fints_xyz[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x + rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (2.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + fe_0 * fe_0 * rpa_z * rpb_z * rpb_x);

        fints_xzz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y);

        fints_yyy[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyy[i] += fss * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] +=
            fss * (fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y + 2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               fe_0 * fe_0 * rpa_y * rpb_z * rpb_y);

        fints_yyz[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] +=
            fss * (2.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y + fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z + fe_0 * fe_0 * rpa_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_yzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * (3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y);

        fints_zzz[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_zzz[i] += fss * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveOverlapGF_YZZZ_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + rpa_z * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x);

        fints_xxy[i] += fss * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y);

        fints_xxz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_x + fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x);

        fints_xyy[i] += fss * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x;

        fints_xyz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x);

        fints_xyz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x + 3.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_xzz[i] += fss * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x;

        fints_yyy[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y);

        fints_yyy[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z + fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_yyz[i] += fss * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y + 3.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_yzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_yzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] +=
            fss * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z + (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y);

        fints_zzz[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveOverlapGF_ZZZZ_T(TDoubleArray&       buffer_xxx,
                              TDoubleArray&       buffer_xxy,
                              TDoubleArray&       buffer_xxz,
                              TDoubleArray&       buffer_xyy,
                              TDoubleArray&       buffer_xyz,
                              TDoubleArray&       buffer_xzz,
                              TDoubleArray&       buffer_yyy,
                              TDoubleArray&       buffer_yyz,
                              TDoubleArray&       buffer_yzz,
                              TDoubleArray&       buffer_zzz,
                              const double        bra_exp,
                              const double        bra_norm,
                              const TPoint3D&     bra_coord,
                              const TDoubleArray& ket_exps,
                              const TDoubleArray& ket_norms,
                              const TDoubleArray& ket_coords_x,
                              const TDoubleArray& ket_coords_y,
                              const TDoubleArray& ket_coords_z,
                              const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_x +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_xxx[i] += fss * rpa_z * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x;

        fints_xxy[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_xxy[i] += fss * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x;

        fints_xxz[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x + 2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_z + 3.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z);

        fints_xxz[i] += fss * (fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x);

        fints_xyy[i] += fss * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x;

        fints_xyz[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x + 2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x +
                               3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x +
                               rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x + 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_x + 6.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x);

        fints_xzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y);

        fints_yyy[i] += fss * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y;

        fints_yyz[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y + 2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_z + 3.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z);

        fints_yyz[i] += fss * (fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y + 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_y + 6.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y);

        fints_yzz[i] += fss * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * (3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z + 6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z +
                               (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpb_z + 9.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z +
                               (27.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z);

        fints_zzz[i] += fss * (3.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z +
                               (15.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpa_z + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_z * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z);
    }
}

}  // namespace ovlrec
