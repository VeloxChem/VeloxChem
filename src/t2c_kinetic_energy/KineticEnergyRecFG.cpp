#include "KineticEnergyRecFG.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "MathConst.hpp"
#include "T2CDistributor.hpp"

namespace kinrec {  // kinrec namespace

auto
compKineticEnergyFG(CSubMatrix*      matrix,
                    const CGtoBlock& bra_gto_block,
                    const CGtoBlock& ket_gto_block,
                    const bool       ang_order,
                    const int64_t    bra_first,
                    const int64_t    bra_last) -> void
{
    // spherical transformation factors

    const double f3_5 = std::sqrt(2.5);

    const double f3_15 = 2.0 * std::sqrt(15.0);

    const double f3_3 = std::sqrt(1.5);

    const double f4_35 = 4.0 * std::sqrt(35);

    const double f4_17 = 4.0 * std::sqrt(17.5);

    const double f4_5 = 4.0 * std::sqrt(5.0);

    const double f4_2 = 4.0 * std::sqrt(2.5);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_XXXX(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_XXXY(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_XXXZ(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_XXYY(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 6.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -6.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 1.50 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 6.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -1.50 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -6.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 1.50 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -6.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 1.50 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -6.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -6.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 1.50 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, 6.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -1.50 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 6.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -1.50 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_XXYZ(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_XXZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_XYYY(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_XYYZ(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_XYZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -6.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 6.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -6.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -6.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -6.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 6.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 6.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_XZZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 4.0 * f4_2 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_YYYY(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -0.5 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -0.25 * f4_35 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 0.5 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 0.25 * f4_35 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.25 * f4_35 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -0.25 * f4_35 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -0.5 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -0.25 * f4_35 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, 3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 0.5 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 0.25 * f4_35 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 0.5 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 0.25 * f4_35 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_YYYZ(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, f4_17 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -f4_17 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_17 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, f4_17 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, f4_17 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -f4_17 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -f4_17 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_YYZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -24.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, 24.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f4_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -24.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f4_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 24.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 24.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f4_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, 24.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f4_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, -24.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f4_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -24.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f4_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_YZZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxy, 4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xxz, 4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xyy, -4.0 * f4_2 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f4_2 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_xzz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f4_2 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f4_2 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yyz, -4.0 * f4_2 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix, buffer_yzz, 4.0 * f4_2 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

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

                    kinrec::compPrimitiveKineticEnergyFG_T_ZZZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, 8.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -8.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, 8.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -8.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -8.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -8.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -8.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, 8.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 8.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);
        }
    }
}

auto
compPrimitiveKineticEnergyFG_T_XXXX(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-9.0 * fbe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 - 6.0 * fbe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 -
                               (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 - 9.0 * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               3.0 * fbe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += fss * (-9.0 * fe_0 * fke_0 * rpa_x * rpb_x * rpb_x * fz_0 - 18.0 * fe_0 * fke_0 * rpa_x * rpa_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_x * rpa_x * rpa_x * fz_0 + 18.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 +
                               72.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += fss * (36.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 - (27.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 -
                               9.0 * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + 135.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 +
                               90.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * fz_0);

        fints_xxx[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * fz_0 + 30.0 * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 +
                               45.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 + 60.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               6.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += fss * 14.0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x + 6.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x +
                               3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x + (27.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               9.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x);

        fints_xxx[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x + 3.0 * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x +
                               (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x + (15.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 -
                               fbe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - 12.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpa_x * fz_0);

        fints_xxy[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_y * rpb_x * rpb_x * fz_0 + 48.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 +
                               36.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0);

        fints_xxy[i] += fss * (60.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 + 15.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 -
                               6.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += fss * 14.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxy[i] += ftt * (4.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x + 3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x + 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x);

        fints_xxy[i] += ftt * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x);

        fints_xxz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 -
                               fbe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - 12.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpa_x * fz_0);

        fints_xxz[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_z * rpb_x * rpb_x * fz_0 + 48.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 +
                               36.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0);

        fints_xxz[i] += fss * (60.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 + 15.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 -
                               6.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += fss * 14.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxz[i] += ftt * (4.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x + 3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x + 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x);

        fints_xxz[i] += ftt * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x);

        fints_xyy[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 - 2.0 * fbe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 - 3.0 * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xyy[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_y * rpa_y * rpa_x * fz_0 - 6.0 * fe_0 * fke_0 * rpa_y * rpa_y * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_x * rpb_x * rpb_x * fz_0 + 36.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * fz_0 +
                               24.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xyy[i] += fss * (6.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 -
                               3.0 * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * fz_0 +
                               30.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * fz_0);

        fints_xyy[i] += fss * (15.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 + 10.0 * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 + 12.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               6.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xyy[i] += ftt * (3.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x + 2.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x +
                               3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x);

        fints_xyy[i] += ftt * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x + fe_0 * fe_0 * rpb_x * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x);

        fints_xyz[i] +=
            fss * (-3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpa_x * fz_0 - 6.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_x * fz_0 +
                   36.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * fz_0 + 24.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * fz_0 +
                   (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * fz_0);

        fints_xyz[i] += fss * (30.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * fz_0 - 6.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xyz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x + 2.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x + 3.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x +
                               rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x);

        fints_xzz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 - 2.0 * fbe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 - 3.0 * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xzz[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpa_x * fz_0 - 6.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_x * rpb_x * rpb_x * fz_0 + 36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * fz_0 +
                               24.0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xzz[i] += fss * (6.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 -
                               3.0 * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * fz_0 +
                               30.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * fz_0);

        fints_xzz[i] += fss * (15.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 + 10.0 * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 + 12.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               6.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xzz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x + 2.0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x +
                               3.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x);

        fints_xzz[i] += ftt * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x + fe_0 * fe_0 * rpb_x * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpb_x);

        fints_yyy[i] += fss * (-9.0 * fbe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 - (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 -
                               3.0 * fbe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - 9.0 * fe_0 * fke_0 * rpa_y * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_y * rpa_y * fz_0);

        fints_yyy[i] += fss * (18.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 +
                               36.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * fz_0);

        fints_yyy[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 - 6.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * fz_0 +
                               14.0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * fz_0);

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x + 3.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyy[i] += ftt * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x;

        fints_yyz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 -
                               fbe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpa_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpb_x * rpb_x * fz_0);

        fints_yyz[i] += fss * (36.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0);

        fints_yyz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 - 6.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * fz_0);

        fints_yyz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_yyz[i] += ftt * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x;

        fints_yzz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 -
                               fbe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpa_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpb_x * rpb_x * fz_0);

        fints_yzz[i] += fss * (36.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * fz_0 + 15.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0);

        fints_yzz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 - 6.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x * fz_0);

        fints_yzz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yzz[i] += ftt * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpb_x;

        fints_zzz[i] += fss * (-9.0 * fbe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 - (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 -
                               3.0 * fbe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 - 9.0 * fe_0 * fke_0 * rpa_z * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpa_z * fz_0);

        fints_zzz[i] += fss * (18.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x * fz_0 +
                               36.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * fz_0);

        fints_zzz[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 - 6.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x * fz_0);

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x + 3.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_zzz[i] += ftt * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x * rpb_x;
    }
}

auto
compPrimitiveKineticEnergyFG_T_XXXY(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 - (9.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xxx[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_y * fz_0 + 18.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 +
                   54.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 + 18.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0);

        fints_xxx[i] += fss * ((135.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 +
                               15.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 - 3.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xxx[i] += fss * 14.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxx[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x +
                   (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x + (27.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y);

        fints_xxx[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_y * fz_0);

        fints_xxy[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                   36.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += fss * (6.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0);

        fints_xxy[i] += fss * (15.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 +
                               9.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xxy[i] +=
            fss * (-3.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 + 14.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += ftt * (3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y);

        fints_xxy[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxy[i] += ftt * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_xxz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - fbe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_x * fz_0 +
                               36.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxz[i] +=
            fss * (18.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 -
                   3.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xxz[i] += fss * 14.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x);

        fints_xxz[i] += ftt * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x;

        fints_xyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_x * fz_0);

        fints_xyy[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_x * fz_0 +
                   12.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xyy[i] += fss * (6.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + 15.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * fz_0 +
                               15.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0);

        fints_xyy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_xyy[i] +=
            fss * (-3.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 + 14.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xyy[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x +
                               (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x);

        fints_xyy[i] += ftt * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xyy[i] += ftt * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_xyz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_x * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xyz[i] += fss * (-(3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0);

        fints_xyz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 + 14.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xyz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x);

        fints_xyz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_xzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_y * fz_0);

        fints_xzz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0);

        fints_xzz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x);

        fints_xzz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_yyy[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_x * fz_0);

        fints_yyy[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0);

        fints_yyy[i] += fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               3.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * fz_0);

        fints_yyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * fz_0;

        fints_yyy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x);

        fints_yyy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_yyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - fbe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_yyz[i] +=
            fss * (18.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 -
                   3.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * fz_0);

        fints_yyz[i] += fss * 14.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * fz_0;

        fints_yyz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x);

        fints_yyz[i] += ftt * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x;

        fints_yzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_x * fz_0);

        fints_yzz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0);

        fints_yzz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * fz_0);

        fints_yzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x * fz_0;

        fints_yzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x);

        fints_yzz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpb_x);

        fints_zzz[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * fz_0);

        fints_zzz[i] +=
            fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - 3.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * fz_0 +
                   14.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_zzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpb_x);
    }
}

auto
compPrimitiveKineticEnergyFG_T_XXXZ(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 - (9.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xxx[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_z * fz_0 + 18.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 +
                   54.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * fz_0 + 18.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0);

        fints_xxx[i] += fss * ((135.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 +
                               15.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 - 3.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xxx[i] += fss * 14.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxx[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x +
                   (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x + (27.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z);

        fints_xxx[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_x * fz_0 +
                               36.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxy[i] +=
            fss * (18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 -
                   3.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xxy[i] += fss * 14.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxy[i] += ftt * (3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_xxy[i] += ftt * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x;

        fints_xxz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_z * fz_0);

        fints_xxz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                   36.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += fss * (6.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0);

        fints_xxz[i] += fss * (15.0 * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 +
                               9.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xxz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 + 14.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z);

        fints_xxz[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxz[i] += ftt * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x);

        fints_xyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_z * fz_0);

        fints_xyy[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0);

        fints_xyy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               3.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xyy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x);

        fints_xyy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x);

        fints_xyz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_y * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_x * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xyz[i] += fss * (-(3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0);

        fints_xyz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 + 14.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xyz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x);

        fints_xyz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x);

        fints_xzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_x * fz_0);

        fints_xzz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_x * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xzz[i] += fss * (6.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * fz_0 +
                               15.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0);

        fints_xzz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_xzz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * fz_0 + 14.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xzz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x +
                               (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x);

        fints_xzz[i] += ftt * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_xzz[i] += ftt * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpb_x);

        fints_yyy[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * fz_0);

        fints_yyy[i] +=
            fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - 3.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * fz_0 +
                   14.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * fz_0);

        fints_yyy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x + rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x);

        fints_yyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_x * fz_0);

        fints_yyz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0);

        fints_yyz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * fz_0);

        fints_yyz[i] += fss * 14.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * fz_0;

        fints_yyz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x);

        fints_yyz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x);

        fints_yzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_yzz[i] +=
            fss * (18.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 -
                   3.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * fz_0);

        fints_yzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x * fz_0;

        fints_yzz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_yzz[i] += ftt * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpb_x;

        fints_zzz[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_x * fz_0);

        fints_zzz[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0);

        fints_zzz[i] += fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * fz_0);

        fints_zzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x * fz_0;

        fints_zzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x);

        fints_zzz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x * rpb_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x * rpb_x);
    }
}

auto
compPrimitiveKineticEnergyFG_T_XXYY(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xxx[i] += fss * (-3.0 * fbe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_x * rpa_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxx[i] +=
            fss * (18.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 + 36.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 -
                   3.0 * fe_0 * fe_0 * fke_0 * rpa_x * fz_0);

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxx[i] += fss * (15.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 +
                               6.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 - fke_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 -
                               fke_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += fss * 14.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0;

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x + 3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y);

        fints_xxx[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxx[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_xxy[i] += fss * (-fbe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 - 2.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_x * fz_0 -
                               fe_0 * fke_0 * rpa_y * rpa_x * rpa_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += fss * (-fe_0 * fke_0 * rpa_x * rpa_x * rpb_y * fz_0 + 24.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += fss * (12.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 - fe_0 * fe_0 * fke_0 * rpa_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + 10.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * fz_0);

        fints_xxy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 +
                               20.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fke_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 - fke_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 +
                               14.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += ftt * (2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x + fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxy[i] += ftt * (fe_0 * fe_0 * rpa_y * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x +
                               2.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x);

        fints_xxy[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 -
                               fbe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 - 2.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_x * fz_0);

        fints_xxz[i] +=
            fss * (-fe_0 * fke_0 * rpa_z * rpa_x * rpa_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_y * fz_0 -
                   (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * rpb_x * fz_0 + 24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xxz[i] += fss * (6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 -
                               fe_0 * fe_0 * fke_0 * rpa_z * fz_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * fz_0);

        fints_xxz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 - fke_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 -
                               fke_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += fss * 14.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0;

        fints_xxz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x + fe_0 * fe_0 * rpa_z * rpa_x * rpb_x);

        fints_xxz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xyy[i] += fss * (-fbe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 - 2.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_y * fz_0 -
                               fe_0 * fke_0 * rpa_y * rpa_y * rpa_x * fz_0 - fe_0 * fke_0 * rpa_y * rpa_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xyy[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * rpb_x * fz_0 + 24.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * fz_0 +
                   12.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += fss * (6.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 - fe_0 * fe_0 * fke_0 * rpa_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + 10.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * fz_0 +
                               20.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * fz_0 + 5.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fke_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * fz_0 - fke_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * fz_0 +
                               14.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xyy[i] += ftt * (2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x + fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x);

        fints_xyy[i] += ftt * (fe_0 * fe_0 * rpa_y * rpa_x * rpb_y + 2.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y);

        fints_xyy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x);

        fints_xyz[i] += fss * (-fe_0 * fke_0 * rpa_z * rpa_y * rpa_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_y * rpb_x * fz_0 -
                               fe_0 * fke_0 * rpa_z * rpa_x * rpb_y * fz_0 + 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xyz[i] += fss * (12.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * fz_0 + 5.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * fz_0);

        fints_xyz[i] +=
            fss * (10.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - fke_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * fz_0 -
                   fke_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * fz_0 + 14.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xyz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x + fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x +
                               fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x);

        fints_xyz[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y +
                               fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x);

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xzz[i] += fss * (-fbe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_z * rpa_x * fz_0 -
                               fe_0 * fke_0 * rpa_z * rpa_z * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xzz[i] += fss * (6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0);

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xzz[i] +=
            fss * (5.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 + fe_0 * fe_0 * fe_0 * rpa_x * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                   fke_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * fz_0 - fke_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x * fz_0;

        fints_xzz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x + fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x);

        fints_xzz[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xzz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpb_x);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_yyy[i] += fss * (-3.0 * fbe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_y * rpb_y * fz_0 - fe_0 * fke_0 * rpa_y * rpa_y * rpa_y * fz_0);

        fints_yyy[i] +=
            fss * (18.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 + 36.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * fz_0 -
                   3.0 * fe_0 * fe_0 * fke_0 * rpa_y * fz_0);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * fz_0);

        fints_yyy[i] += fss * (15.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 +
                               6.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 - fke_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * fz_0 -
                               fke_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * fz_0);

        fints_yyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * fz_0;

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y);

        fints_yyy[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 -
                               fbe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 - 2.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_y * fz_0);

        fints_yyz[i] +=
            fss * (-fe_0 * fke_0 * rpa_z * rpa_y * rpa_y * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_y * fz_0 -
                   (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * rpb_x * fz_0 + 24.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += fss * (6.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 -
                               fe_0 * fe_0 * fke_0 * rpa_z * fz_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * fz_0);

        fints_yyz[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 - fke_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * fz_0 -
                               fke_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * fz_0);

        fints_yyz[i] += fss * 14.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * fz_0;

        fints_yyz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x + fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_yyz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_yzz[i] += fss * (-fbe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_z * rpa_y * fz_0 -
                               fe_0 * fke_0 * rpa_z * rpa_z * rpb_y * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * rpb_x * fz_0);

        fints_yzz[i] += fss * (6.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0);

        fints_yzz[i] +=
            fss * (5.0 * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 + fe_0 * fe_0 * fe_0 * rpa_y * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                   fke_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * fz_0 - fke_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * fz_0);

        fints_yzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x * fz_0;

        fints_yzz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x + fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y);

        fints_yzz[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yzz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpb_x);

        fints_zzz[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 -
                   (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 -
                   (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_y * fz_0);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_z * rpa_z * fz_0 +
                               18.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * fz_0);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0);

        fints_zzz[i] += fss * (-fke_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * fz_0 - fke_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x * fz_0);

        fints_zzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x);

        fints_zzz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * rpb_x);
    }
}

auto
compPrimitiveKineticEnergyFG_T_XXYZ(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 -
                   3.0 * fbe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_y * fz_0 +
                   18.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += fss * (36.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 +
                               15.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 - fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xxx[i] += fss * 14.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0;

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x + 3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x);

        fints_xxx[i] += ftt * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x;

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_y * fz_0);

        fints_xxy[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_z * fz_0 + 24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 +
                               10.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 - fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 +
                               14.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += ftt * (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y);

        fints_xxy[i] += ftt * (fe_0 * fe_0 * rpa_x * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_y * fz_0);

        fints_xxz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_y * fz_0 + 24.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 +
                               10.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 - fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 +
                               14.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y);

        fints_xxz[i] += ftt * (fe_0 * fe_0 * rpa_x * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_y * rpa_x * rpb_z * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xyy[i] +=
            fss * (12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 +
                   12.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 +
                   5.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * fz_0);

        fints_xyy[i] += fss * (10.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 - fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 +
                               14.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xyy[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y +
                               fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z);

        fints_xyy[i] += ftt * (fe_0 * fe_0 * rpa_y * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x + rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x);

        fints_xyz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_x * rpb_z * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_x * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 + 12.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xyz[i] += fss * (6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 - (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * fz_0 + 5.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * fz_0);

        fints_xyz[i] += fss * (5.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 +
                               fe_0 * fe_0 * fe_0 * rpa_x * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xyz[i] += fss * 14.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0;

        fints_xyz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y + fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z);

        fints_xyz[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xyz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x);

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_x * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xzz[i] +=
            fss * (12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 +
                   5.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * fz_0);

        fints_xzz[i] += fss * (10.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 - fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * fz_0 +
                               14.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xzz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y +
                               fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y);

        fints_xzz[i] += ftt * (fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x + rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpb_x);

        fints_yyy[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 -
                   (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 -
                   (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_y * fz_0);

        fints_yyy[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_z * fz_0 + 18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0);

        fints_yyy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * fz_0);

        fints_yyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * fz_0;

        fints_yyy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z);

        fints_yyy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_y * rpb_z * fz_0);

        fints_yyz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_y * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_y * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * fz_0);

        fints_yyz[i] += fss * (6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + 5.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0);

        fints_yyz[i] += fss * (5.0 * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 +
                               fe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_yyz[i] +=
            fss * (-fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 + 14.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * fz_0);

        fints_yyz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_yyz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyz[i] += ftt * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_y * rpb_y * fz_0);

        fints_yzz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_z * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_y * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x * fz_0);

        fints_yzz[i] += fss * (6.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + 5.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0);

        fints_yzz[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 +
                               fe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_yzz[i] +=
            fss * (-fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * fz_0 + 14.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x * fz_0);

        fints_yzz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_yzz[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_yzz[i] += ftt * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpb_x);

        fints_zzz[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 -
                   (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 -
                   (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_y * fz_0);

        fints_zzz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_y * fz_0 + 18.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0);

        fints_zzz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * fz_0);

        fints_zzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x * fz_0;

        fints_zzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y);

        fints_zzz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * rpb_x);
    }
}

auto
compPrimitiveKineticEnergyFG_T_XXZZ(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xxx[i] += fss * (-3.0 * fbe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_x * rpa_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxx[i] +=
            fss * (18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 + 36.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 -
                   3.0 * fe_0 * fe_0 * fke_0 * rpa_x * fz_0);

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxx[i] += fss * (15.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 +
                               6.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 - fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 -
                               fke_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += fss * 14.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0;

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x + 3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z);

        fints_xxx[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxx[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 - 2.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_x * fz_0);

        fints_xxy[i] +=
            fss * (-fe_0 * fke_0 * rpa_y * rpa_x * rpa_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_z * fz_0 -
                   (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * rpb_x * fz_0 + 24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * fz_0);

        fints_xxy[i] += fss * (6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 -
                               fe_0 * fe_0 * fke_0 * rpa_y * fz_0 + 10.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * fz_0);

        fints_xxy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 - fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 -
                               fke_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += fss * 14.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0;

        fints_xxy[i] += ftt * (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + fe_0 * fe_0 * rpa_y * rpa_x * rpb_x);

        fints_xxy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_xxz[i] += fss * (-fbe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 - 2.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_x * fz_0 -
                               fe_0 * fke_0 * rpa_z * rpa_x * rpa_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += fss * (-fe_0 * fke_0 * rpa_x * rpa_x * rpb_z * fz_0 + 24.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += fss * (12.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * fz_0 - fe_0 * fe_0 * fke_0 * rpa_z * fz_0 -
                               (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * fz_0);

        fints_xxz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 +
                               20.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 + 5.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 - fke_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x + fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] += ftt * (fe_0 * fe_0 * rpa_z * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x +
                               2.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x);

        fints_xxz[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xyy[i] += fss * (-fbe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_y * rpa_y * rpa_x * fz_0 -
                               fe_0 * fke_0 * rpa_y * rpa_y * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xyy[i] += fss * (6.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xyy[i] +=
            fss * (5.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 + fe_0 * fe_0 * fe_0 * rpa_x * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                   fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 - fke_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0;

        fints_xyy[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x + fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x);

        fints_xyy[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xyy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xyz[i] += fss * (-fe_0 * fke_0 * rpa_z * rpa_y * rpa_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_y * rpb_x * fz_0 -
                               fe_0 * fke_0 * rpa_y * rpa_x * rpb_z * fz_0 + 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 +
                               6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xyz[i] += fss * (12.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * fz_0 + 5.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * fz_0);

        fints_xyz[i] +=
            fss * (10.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 -
                   fke_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * fz_0 + 14.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xyz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x + fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x +
                               fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x);

        fints_xyz[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               fe_0 * fe_0 * rpa_y * rpb_z * rpb_x + rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xzz[i] += fss * (-fbe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 - 2.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_z * fz_0 -
                               fe_0 * fke_0 * rpa_z * rpa_z * rpa_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_z * fz_0);

        fints_xzz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * rpb_x * fz_0 + 24.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xzz[i] += fss * (6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 - fe_0 * fe_0 * fke_0 * rpa_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * fz_0 +
                               20.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0);

        fints_xzz[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * fz_0 + 5.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xzz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * fz_0 - fke_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xzz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x + fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_xzz[i] += ftt * (fe_0 * fe_0 * rpa_z * rpa_x * rpb_z + 2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z);

        fints_xzz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_yyy[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 -
                   (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 -
                   (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_z * fz_0);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_y * rpa_y * rpa_y * fz_0 +
                               18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * fz_0);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0);

        fints_yyy[i] += fss * (-fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 - fke_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * fz_0 +
                               14.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * fz_0);

        fints_yyy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x);

        fints_yyy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_yyz[i] += fss * (-fbe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_y * rpa_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_z * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * rpb_x * fz_0 -
                               fe_0 * fke_0 * rpa_y * rpa_y * rpb_z * fz_0);

        fints_yyz[i] += fss * (6.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * fz_0);

        fints_yyz[i] +=
            fss * (5.0 * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 + fe_0 * fe_0 * fe_0 * rpa_z * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                   fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 - fke_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * fz_0);

        fints_yyz[i] += fss * 14.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * fz_0;

        fints_yyz[i] +=
            ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x + fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x +
                   (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y);

        fints_yyz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_yyz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 - 2.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_z * fz_0);

        fints_yzz[i] +=
            fss * (-fe_0 * fke_0 * rpa_z * rpa_z * rpa_y * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_z * fz_0 -
                   (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * rpb_x * fz_0 + 24.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * fz_0);

        fints_yzz[i] += fss * (6.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 -
                               fe_0 * fe_0 * fke_0 * rpa_y * fz_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * fz_0);

        fints_yzz[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 - fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * fz_0);

        fints_yzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x * fz_0;

        fints_yzz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x + fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_yzz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpb_x);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_zzz[i] += fss * (-3.0 * fbe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpb_z * fz_0 - fe_0 * fke_0 * rpa_z * rpa_z * rpa_z * fz_0);

        fints_zzz[i] +=
            fss * (18.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x * fz_0 + 36.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * fz_0 -
                   3.0 * fe_0 * fe_0 * fke_0 * rpa_z * fz_0);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * fz_0);

        fints_zzz[i] += fss * (15.0 * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 +
                               6.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 - fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * fz_0);

        fints_zzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x * fz_0;

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x + 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z);

        fints_zzz[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x * rpb_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_zzz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x * rpb_x);
    }
}

auto
compPrimitiveKineticEnergyFG_T_XYYY(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xxx[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_y * fz_0 + 18.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * fz_0 + 18.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0);

        fints_xxx[i] += fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               3.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xxx[i] += fss * 14.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0;

        fints_xxx[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y);

        fints_xxx[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x);

        fints_xxy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_y * fz_0);

        fints_xxy[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                   12.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * fz_0 + 18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xxy[i] += fss * (18.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0);

        fints_xxy[i] += fss * (15.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xxy[i] +=
            fss * (-3.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 + 14.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xxy[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x +
                               (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y);

        fints_xxy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxy[i] += ftt * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x);

        fints_xxz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - fbe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xxz[i] +=
            fss * (18.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 -
                   3.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xxz[i] += fss * 14.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0;

        fints_xxz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x);

        fints_xxz[i] += ftt * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x;

        fints_xyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_x * fz_0);

        fints_xyy[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_x * fz_0 +
                   36.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xyy[i] += fss * (6.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + 15.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * fz_0 +
                               15.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0);

        fints_xyy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 +
                               9.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_xyy[i] +=
            fss * (-3.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 + 14.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += ftt * (3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x);

        fints_xyy[i] += ftt * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xyy[i] += ftt * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_x * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += fss * (-(3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0);

        fints_xyz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 + 14.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x);

        fints_xyz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x);

        fints_xzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_y * fz_0);

        fints_xzz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0);

        fints_xzz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x * fz_0;

        fints_xzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x);

        fints_xzz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_x);

        fints_yyy[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 - (9.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_x * fz_0);

        fints_yyy[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 +
                   54.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0);

        fints_yyy[i] += fss * ((135.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 +
                               15.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 - 3.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * fz_0);

        fints_yyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * fz_0;

        fints_yyy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x + (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x);

        fints_yyy[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x);

        fints_yyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - fbe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_x * fz_0 +
                               36.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * fz_0);

        fints_yyz[i] +=
            fss * (18.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 -
                   3.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * fz_0);

        fints_yyz[i] += fss * 14.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * fz_0;

        fints_yyz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x);

        fints_yyz[i] += ftt * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x;

        fints_yzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_x * fz_0);

        fints_yzz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0);

        fints_yzz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * fz_0);

        fints_yzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x * fz_0;

        fints_yzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x);

        fints_yzz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_x);

        fints_zzz[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * fz_0);

        fints_zzz[i] +=
            fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - 3.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * fz_0 +
                   14.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x * fz_0);

        fints_zzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_x);
    }
}

auto
compPrimitiveKineticEnergyFG_T_XYYZ(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 -
                   (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 -
                   (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xxx[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_z * fz_0 + 18.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0);

        fints_xxx[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xxx[i] += fss * 14.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0;

        fints_xxx[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z);

        fints_xxx[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 - fe_0 * fke_0 * rpa_y * rpa_x * rpb_z * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_x * fz_0);

        fints_xxy[i] += fss * (12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * fz_0);

        fints_xxy[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 + 10.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 - fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 +
                               14.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xxy[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x + fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z);

        fints_xxy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x + fe_0 * fe_0 * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x + rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_x * rpb_z * fz_0);

        fints_xxz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xxz[i] += fss * (6.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + 5.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0);

        fints_xxz[i] += fss * (5.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 +
                               fe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xxz[i] +=
            fss * (-fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 + 14.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xxz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z);

        fints_xxz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxz[i] += ftt * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_z * fz_0);

        fints_xyy[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_x * fz_0 + 24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + 10.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0);

        fints_xyy[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 - fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 +
                               14.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += ftt * (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x + fe_0 * fe_0 * rpa_y * rpb_z * rpb_y);

        fints_xyy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_y * rpb_z * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_x * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += fss * (6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 - (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * fz_0 + 5.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * fz_0);

        fints_xyz[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 + 5.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 +
                               fe_0 * fe_0 * fe_0 * rpa_y * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xyz[i] += fss * 14.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0;

        fints_xyz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y + fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_xyz[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xyz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x);

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_x * rpb_x * fz_0);

        fints_xzz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_z * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_x * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * fz_0);

        fints_xzz[i] += fss * (6.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + 5.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0);

        fints_xzz[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 +
                               fe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_xzz[i] +=
            fss * (-fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * fz_0 + 14.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xzz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x);

        fints_xzz[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_xzz[i] += ftt * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_x);

        fints_yyy[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 -
                   3.0 * fbe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * fz_0);

        fints_yyy[i] += fss * (36.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 +
                               15.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 - fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * fz_0);

        fints_yyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * fz_0;

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x);

        fints_yyy[i] += ftt * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x;

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_x * fz_0);

        fints_yyz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_x * fz_0 + 24.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * fz_0);

        fints_yyz[i] += fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 +
                               10.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0);

        fints_yyz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 - fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * fz_0);

        fints_yyz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x);

        fints_yyz[i] += ftt * (fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_x * fz_0);

        fints_yzz[i] +=
            fss * (12.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 +
                   5.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * fz_0);

        fints_yzz[i] += fss * (10.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 - fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x * fz_0);

        fints_yzz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x +
                               fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x);

        fints_yzz[i] += ftt * (fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x + rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_x);

        fints_zzz[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 -
                   (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 -
                   (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_x * fz_0);

        fints_zzz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0);

        fints_zzz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * fz_0);

        fints_zzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x * fz_0;

        fints_zzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x);

        fints_zzz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_x);
    }
}

auto
compPrimitiveKineticEnergyFG_T_XYZZ(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 -
                   (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 -
                   (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xxx[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_y * fz_0 + 18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0);

        fints_xxx[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fke_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xxx[i] += fss * 14.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0;

        fints_xxx[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y);

        fints_xxx[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 - fe_0 * fke_0 * rpa_y * rpa_x * rpb_y * fz_0);

        fints_xxy[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                   12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xxy[i] += fss * (6.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + 5.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0);

        fints_xxy[i] += fss * (5.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 +
                               fe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xxy[i] +=
            fss * (-fke_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 + 14.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xxy[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y);

        fints_xxy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxy[i] += ftt * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_x * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_x * fz_0);

        fints_xxz[i] += fss * (12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * fz_0);

        fints_xxz[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 + 10.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 - fke_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xxz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x + fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y);

        fints_xxz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + fe_0 * fe_0 * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x + rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 - fe_0 * fke_0 * rpa_y * rpa_x * rpb_x * fz_0);

        fints_xyy[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_y * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_x * fz_0 +
                   12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * fz_0);

        fints_xyy[i] += fss * (6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + 5.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0);

        fints_xyy[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 +
                               fe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_xyy[i] +=
            fss * (-fke_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 + 14.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x);

        fints_xyy[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xyy[i] += ftt * ((1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x);

        fints_xyz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_x * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xyz[i] += fss * (12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 - (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0);

        fints_xyz[i] += fss * (5.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 + 5.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 +
                               fe_0 * fe_0 * fe_0 * rpa_z * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fke_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += fss * 14.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0;

        fints_xyz[i] +=
            ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x + fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x +
                   (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_xyz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_xyz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_y * fz_0);

        fints_xzz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_x * fz_0 + 24.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xzz[i] += fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0);

        fints_xzz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 - fke_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xzz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x + fe_0 * fe_0 * rpa_z * rpb_z * rpb_y);

        fints_xzz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_x);

        fints_yyy[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 -
                   (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 -
                   (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_x * fz_0);

        fints_yyy[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0);

        fints_yyy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fke_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * fz_0);

        fints_yyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * fz_0;

        fints_yyy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x);

        fints_yyy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_x * fz_0);

        fints_yyz[i] += fss * (12.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * fz_0);

        fints_yyz[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 + 10.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 - fke_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * fz_0);

        fints_yyz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x + fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x);

        fints_yyz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x + fe_0 * fe_0 * rpa_y * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x + rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_x * fz_0);

        fints_yzz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_x * fz_0 + 24.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * fz_0);

        fints_yzz[i] += fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0);

        fints_yzz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 - fke_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * fz_0 +
                               14.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x * fz_0);

        fints_yzz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x + fe_0 * fe_0 * rpa_z * rpb_z * rpb_x);

        fints_yzz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_x);

        fints_zzz[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 -
                   3.0 * fbe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * fz_0);

        fints_zzz[i] += fss * (36.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * fz_0 +
                               6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x * fz_0 +
                               15.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x * fz_0 - fke_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * fz_0);

        fints_zzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x * fz_0;

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x + 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_x);

        fints_zzz[i] += ftt * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * rpb_x;
    }
}

auto
compPrimitiveKineticEnergyFG_T_XZZZ(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xxx[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_z * fz_0 + 18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 + 18.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0);

        fints_xxx[i] += fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               3.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xxx[i] += fss * 14.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0;

        fints_xxx[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z +
                   (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z);

        fints_xxx[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x);

        fints_xxy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_x * fz_0 +
                               12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * fz_0);

        fints_xxy[i] +=
            fss * (18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 -
                   3.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xxy[i] += fss * 14.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0;

        fints_xxy[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_xxy[i] += ftt * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x;

        fints_xxz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_z * fz_0);

        fints_xxz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 + 18.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xxz[i] += fss * (18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0);

        fints_xxz[i] += fss * (15.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0);

        fints_xxz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 + 14.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xxz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x +
                               (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z);

        fints_xxz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxz[i] += ftt * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x + rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x);

        fints_xyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_z * fz_0);

        fints_xyy[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0);

        fints_xyy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               3.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0;

        fints_xyy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x);

        fints_xyy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x);

        fints_xyz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_y * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_x * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * fz_0 +
                   18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xyz[i] += fss * (-(3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0);

        fints_xyz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 + 14.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xyz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_x);

        fints_xyz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x);

        fints_xzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_x * fz_0);

        fints_xzz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_x * fz_0 +
                   36.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_xzz[i] += fss * (6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x * fz_0 +
                               15.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0);

        fints_xzz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 +
                               9.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_xzz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * fz_0 + 14.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xzz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_x);

        fints_xzz[i] += ftt * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_xzz[i] += ftt * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 +
                   18.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * fz_0);

        fints_yyy[i] +=
            fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - 3.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * fz_0 +
                   14.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * fz_0);

        fints_yyy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x);

        fints_yyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_x * fz_0);

        fints_yyz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * fz_0 +
                   6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 + 18.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0);

        fints_yyz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * fz_0);

        fints_yyz[i] += fss * 14.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * fz_0;

        fints_yyz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_x);

        fints_yyz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x);

        fints_yzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_x * fz_0 +
                               36.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * fz_0);

        fints_yzz[i] +=
            fss * (18.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x * fz_0 -
                   3.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * fz_0);

        fints_yzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x * fz_0;

        fints_yzz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_x);

        fints_yzz[i] += ftt * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_x;

        fints_zzz[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 - (9.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_x * fz_0);

        fints_zzz[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x * fz_0 +
                   54.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_x * fz_0);

        fints_zzz[i] += fss * ((135.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x * fz_0 +
                               15.0 * fe_0 * fe_0 * fe_0 * rpb_x * fz_0 - 3.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * fz_0);

        fints_zzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x * fz_0;

        fints_zzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x + (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_x +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_x);

        fints_zzz[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_x +
                               rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * rpb_x);
    }
}

auto
compPrimitiveKineticEnergyFG_T_YYYY(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-9.0 * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 - (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 -
                               3.0 * fbe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 - 9.0 * fe_0 * fke_0 * rpa_x * rpb_y * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxx[i] += fss * (18.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 +
                               36.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxx[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 - 6.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 +
                               14.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y + 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxx[i] += ftt * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y;

        fints_xxy[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 - 2.0 * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 - 3.0 * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xxy[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpa_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_y * rpb_y * rpb_y * fz_0 -
                               6.0 * fe_0 * fke_0 * rpa_x * rpa_x * rpb_y * fz_0 + 36.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xxy[i] += fss * (24.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 -
                               3.0 * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * fz_0 +
                               15.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0);

        fints_xxy[i] += fss * (30.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * fz_0 + 10.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 + 12.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               6.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xxy[i] += fss * 14.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0;

        fints_xxy[i] += ftt * (3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y +
                               2.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y);

        fints_xxy[i] +=
            ftt * (3.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y + fe_0 * fe_0 * rpb_y * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                   (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y);

        fints_xxz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 -
                               fbe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpa_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpb_y * rpb_y * fz_0);

        fints_xxz[i] += fss * (36.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0);

        fints_xxz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 - 6.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 +
                               14.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xxz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_xxz[i] += ftt * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y;

        fints_xyy[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 -
                               fbe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 - 12.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_y * rpa_x * fz_0);

        fints_xyy[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_x * rpb_y * rpb_y * fz_0 + 48.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * fz_0 +
                               36.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0);

        fints_xyy[i] += fss * (60.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 + 15.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 -
                               6.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0;

        fints_xyy[i] += ftt * (4.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y + 3.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y + 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x);

        fints_xyy[i] += ftt * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y);

        fints_xyz[i] +=
            fss * (-3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpa_x * fz_0 - 6.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_y * fz_0 +
                   36.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * fz_0 + 24.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * fz_0 +
                   (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * fz_0);

        fints_xyz[i] += fss * (30.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * fz_0 - 6.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * fz_0 +
                               14.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xyz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y + 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x + 3.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y +
                               rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y);

        fints_xzz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 -
                               fbe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpa_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xzz[i] += fss * (36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xzz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 - 6.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * fz_0 +
                               14.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xzz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xzz[i] += ftt * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpb_y;

        fints_yyy[i] += fss * (-9.0 * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 - 6.0 * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 -
                               (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 - 9.0 * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               3.0 * fbe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += fss * (-9.0 * fe_0 * fke_0 * rpa_y * rpb_y * rpb_y * fz_0 - 18.0 * fe_0 * fke_0 * rpa_y * rpa_y * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_y * rpa_y * fz_0 + 18.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 +
                               72.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += fss * (36.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * fz_0 - (27.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 -
                               9.0 * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + 135.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 +
                               90.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * fz_0);

        fints_yyy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * fz_0 + 30.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 +
                               45.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 + 60.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               6.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * fz_0;

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y + 6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y +
                               3.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y + (27.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y +
                               9.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y);

        fints_yyy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + 3.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y +
                               (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y + (15.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 -
                               fbe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 - 12.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpa_y * fz_0);

        fints_yyz[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_z * rpb_y * rpb_y * fz_0 + 48.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * fz_0 +
                               36.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0);

        fints_yyz[i] += fss * (60.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 + 15.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 -
                               6.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += fss * 14.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * fz_0;

        fints_yyz[i] += ftt * (4.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y + 3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y + 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y);

        fints_yyz[i] += ftt * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                               rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y);

        fints_yzz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 - 2.0 * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 - 3.0 * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yzz[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpa_y * fz_0 - 6.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpb_y * rpb_y * fz_0 + 36.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * fz_0 +
                               24.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yzz[i] += fss * (6.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 -
                               3.0 * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * fz_0 +
                               30.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * fz_0);

        fints_yzz[i] += fss * (15.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 + 10.0 * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 + 12.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               6.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * fz_0);

        fints_yzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y * fz_0;

        fints_yzz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y + 2.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y +
                               3.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y);

        fints_yzz[i] += ftt * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + fe_0 * fe_0 * rpb_y * rpb_y * rpb_y +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpb_y);

        fints_zzz[i] += fss * (-9.0 * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 - (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 -
                               3.0 * fbe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 - 9.0 * fe_0 * fke_0 * rpa_z * rpb_y * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpa_z * fz_0);

        fints_zzz[i] += fss * (18.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * fz_0 +
                               36.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * fz_0);

        fints_zzz[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 - 6.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * fz_0 +
                               14.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y * fz_0);

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y + 3.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_zzz[i] += ftt * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpb_y;
    }
}

auto
compPrimitiveKineticEnergyFG_T_YYYZ(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_y * fz_0 + 18.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 +
                   18.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xxx[i] +=
            fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - 3.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 +
                   14.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xxx[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y + rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y);

        fints_xxy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_y * fz_0);

        fints_xxy[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_z * fz_0 + 18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 + 18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0);

        fints_xxy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               3.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xxy[i] += fss * 14.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0;

        fints_xxy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z);

        fints_xxy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y);

        fints_xxz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_y * fz_0);

        fints_xxz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_y * fz_0 + 18.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0);

        fints_xxz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xxz[i] += fss * 14.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0;

        fints_xxz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y);

        fints_xxz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y);

        fints_xyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_y * fz_0 +
                               36.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * fz_0);

        fints_xyy[i] +=
            fss * (18.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 -
                   3.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0;

        fints_xyy[i] += ftt * (3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y);

        fints_xyy[i] += ftt * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y;

        fints_xyz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_x * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_x * rpb_y * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 + 18.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xyz[i] += fss * (-(3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0);

        fints_xyz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 + 14.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xyz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y);

        fints_xyz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y);

        fints_xzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_y * fz_0 +
                               12.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * fz_0);

        fints_xzz[i] +=
            fss * (18.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 -
                   3.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y * fz_0;

        fints_xzz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y);

        fints_xzz[i] += ftt * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpb_y;

        fints_yyy[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 - (9.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_y * fz_0);

        fints_yyy[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_z * fz_0 + 18.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 +
                   54.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * fz_0 + 18.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0);

        fints_yyy[i] += fss * ((135.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 +
                               15.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 - 3.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * fz_0);

        fints_yyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * fz_0;

        fints_yyy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y + (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z);

        fints_yyy[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_z * fz_0);

        fints_yyz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_y * fz_0 +
                   36.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * fz_0 + 18.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += fss * (6.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0);

        fints_yyz[i] += fss * (15.0 * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 +
                               9.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_yyz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 + 14.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_yyz[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyz[i] += ftt * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y);

        fints_yzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_y * fz_0);

        fints_yzz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_y * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * fz_0);

        fints_yzz[i] += fss * (6.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * fz_0 +
                               15.0 * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0);

        fints_yzz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_yzz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * fz_0 + 14.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yzz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y +
                               (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_yzz[i] += ftt * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_yzz[i] += ftt * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpb_y);

        fints_zzz[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_y * fz_0);

        fints_zzz[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_y * fz_0 + 18.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0);

        fints_zzz[i] += fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * fz_0);

        fints_zzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y * fz_0;

        fints_zzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y);

        fints_zzz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y * rpb_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpb_y);
    }
}

auto
compPrimitiveKineticEnergyFG_T_YYZZ(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] +=
            fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 -
                   (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 -
                   (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_z * fz_0);

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_y * fz_0 - fe_0 * fke_0 * rpa_x * rpa_x * rpa_x * fz_0 +
                               18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0);

        fints_xxx[i] += fss * (-fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 - fke_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 +
                               14.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0);

        fints_xxx[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y);

        fints_xxx[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_xxy[i] += fss * (-fbe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 - fe_0 * fke_0 * rpa_y * rpa_x * rpa_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_z * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_y * fz_0 -
                               fe_0 * fke_0 * rpa_x * rpa_x * rpb_y * fz_0);

        fints_xxy[i] += fss * (6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 +
                               12.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * fz_0);

        fints_xxy[i] +=
            fss * (5.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 + fe_0 * fe_0 * fe_0 * rpa_y * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                   fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 - fke_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xxy[i] += fss * 14.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0;

        fints_xxy[i] +=
            ftt * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y +
                   (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x);

        fints_xxy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xxy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_xxz[i] += fss * (-fbe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 - fe_0 * fke_0 * rpa_z * rpa_x * rpa_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_z * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_y * fz_0 -
                               fe_0 * fke_0 * rpa_x * rpa_x * rpb_z * fz_0);

        fints_xxz[i] += fss * (6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 +
                               12.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 - (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * fz_0);

        fints_xxz[i] +=
            fss * (5.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 + fe_0 * fe_0 * fe_0 * rpa_z * fz_0 + 2.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                   fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 - fke_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xxz[i] += fss * 14.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0;

        fints_xxz[i] +=
            ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y + fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y +
                   (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x);

        fints_xxz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y +
                               (1.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_xxz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 - 2.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_y * fz_0);

        fints_xyy[i] +=
            fss * (-fe_0 * fke_0 * rpa_y * rpa_y * rpa_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_z * fz_0 -
                   (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_y * fz_0 + 24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * fz_0);

        fints_xyy[i] += fss * (6.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 -
                               fe_0 * fe_0 * fke_0 * rpa_x * fz_0 + 10.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * fz_0);

        fints_xyy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 - fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 -
                               fke_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0;

        fints_xyy[i] += ftt * (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y + fe_0 * fe_0 * rpa_y * rpa_x * rpb_y);

        fints_xyy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_xyz[i] += fss * (-fe_0 * fke_0 * rpa_z * rpa_y * rpa_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_x * rpb_y * fz_0 -
                               fe_0 * fke_0 * rpa_y * rpa_x * rpb_z * fz_0 + 6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 +
                               6.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xyz[i] += fss * (12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 +
                               12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * fz_0 +
                               5.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * fz_0 + 5.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * fz_0);

        fints_xyz[i] +=
            fss * (10.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 -
                   fke_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * fz_0 + 14.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0);

        fints_xyz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y + fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y +
                               fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x);

        fints_xyz[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               fe_0 * fe_0 * rpa_x * rpb_z * rpb_y + rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 - (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 - 2.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_z * fz_0);

        fints_xzz[i] +=
            fss * (-fe_0 * fke_0 * rpa_z * rpa_z * rpa_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_z * fz_0 -
                   (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * rpb_y * fz_0 + 24.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * fz_0);

        fints_xzz[i] += fss * (6.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 -
                               fe_0 * fe_0 * fke_0 * rpa_x * fz_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * fz_0);

        fints_xzz[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 - fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y * fz_0;

        fints_xzz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y + fe_0 * fe_0 * rpa_z * rpa_x * rpb_z);

        fints_xzz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_yyy[i] += fss * (-3.0 * fbe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_y * rpb_y * fz_0 - fe_0 * fke_0 * rpa_y * rpa_y * rpa_y * fz_0);

        fints_yyy[i] +=
            fss * (18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 + 36.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * fz_0 -
                   3.0 * fe_0 * fe_0 * fke_0 * rpa_y * fz_0);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 + 15.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * fz_0);

        fints_yyy[i] += fss * (15.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 +
                               6.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 - fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 -
                               fke_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * fz_0;

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_yyy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_yyz[i] += fss * (-fbe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 - 2.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_y * fz_0 -
                               fe_0 * fke_0 * rpa_z * rpa_y * rpa_y * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += fss * (-fe_0 * fke_0 * rpa_y * rpa_y * rpb_z * fz_0 + 24.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * fz_0 +
                               6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += fss * (12.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * fz_0 - fe_0 * fe_0 * fke_0 * rpa_z * fz_0 -
                               (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * fz_0);

        fints_yyz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 + (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 +
                               20.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 + 5.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 - fke_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * fz_0 +
                               14.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y + fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += ftt * (fe_0 * fe_0 * rpa_z * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y +
                               2.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y);

        fints_yyz[i] += ftt * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_yzz[i] += fss * (-fbe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 - 2.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_z * fz_0 -
                               fe_0 * fke_0 * rpa_z * rpa_z * rpa_y * fz_0 - fe_0 * fke_0 * rpa_z * rpa_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_z * fz_0);

        fints_yzz[i] +=
            fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * rpb_y * fz_0 + 24.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yzz[i] += fss * (6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 - fe_0 * fe_0 * fke_0 * rpa_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + 10.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * fz_0 +
                               20.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0);

        fints_yzz[i] += fss * ((5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * fz_0 + 5.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y * fz_0 +
                               5.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yzz[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * fz_0 - fke_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * fz_0 +
                               14.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y * fz_0);

        fints_yzz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y + fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += ftt * (fe_0 * fe_0 * rpa_z * rpa_y * rpb_z + 2.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z);

        fints_yzz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpb_y);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_zzz[i] += fss * (-3.0 * fbe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpb_z * fz_0 - fe_0 * fke_0 * rpa_z * rpa_z * rpa_z * fz_0);

        fints_zzz[i] +=
            fss * (18.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * fz_0 + 36.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * fz_0 -
                   3.0 * fe_0 * fe_0 * fke_0 * rpa_z * fz_0);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * fz_0);

        fints_zzz[i] += fss * (15.0 * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 +
                               6.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 - fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * fz_0);

        fints_zzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y * fz_0;

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z);

        fints_zzz[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y * rpb_y +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_zzz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * rpb_y);
    }
}

auto
compPrimitiveKineticEnergyFG_T_YZZZ(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_y * fz_0 + 18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 +
                   18.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xxx[i] +=
            fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - 3.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 +
                   14.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0);

        fints_xxx[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y + rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y);

        fints_xxy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_y * fz_0);

        fints_xxy[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_z * fz_0 + 18.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0);

        fints_xxy[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               3.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xxy[i] += fss * 14.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0;

        fints_xxy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z);

        fints_xxy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y);

        fints_xxz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_y * fz_0);

        fints_xxz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * rpb_y * fz_0 + 18.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 + 18.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 -
                   (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0);

        fints_xxz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 + 3.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xxz[i] += fss * 14.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0;

        fints_xxz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpb_y);

        fints_xxz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y);

        fints_xyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_y * fz_0 +
                               12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * fz_0);

        fints_xyy[i] +=
            fss * (18.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 -
                   3.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0;

        fints_xyy[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y);

        fints_xyy[i] += ftt * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y;

        fints_xyz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_x * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_x * rpb_y * fz_0 +
                   18.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 +
                   18.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * fz_0);

        fints_xyz[i] += fss * (-(3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0);

        fints_xyz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 + 14.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0);

        fints_xyz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z +
                   (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpb_y);

        fints_xyz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y);

        fints_xzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * rpb_y * fz_0 +
                               36.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * fz_0);

        fints_xzz[i] +=
            fss * (18.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 +
                   15.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y * fz_0 -
                   3.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * fz_0);

        fints_xzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y * fz_0;

        fints_xzz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpb_y +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_y);

        fints_xzz[i] += ftt * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpb_y;

        fints_yyy[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_y * fz_0);

        fints_yyy[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_z * fz_0 + 18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 +
                   18.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * fz_0 + 18.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0);

        fints_yyy[i] += fss * ((45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 + 9.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               3.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * fz_0);

        fints_yyy[i] += fss * 14.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * fz_0;

        fints_yyy[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z +
                   (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z);

        fints_yyy[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y);

        fints_yyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 -
                               (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_z * fz_0);

        fints_yyz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * rpb_y * fz_0 +
                   12.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * fz_0 + 18.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yyz[i] += fss * (18.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0);

        fints_yyz[i] += fss * (15.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0);

        fints_yyz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 + 14.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yyz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y +
                               (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z);

        fints_yyz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyz[i] += ftt * ((3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y + rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_y * fz_0);

        fints_yzz[i] +=
            fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * rpb_y * fz_0 +
                   36.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * fz_0 +
                   6.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_yzz[i] += fss * (6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + 15.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y * fz_0 +
                               15.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0);

        fints_yzz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y * fz_0 +
                               (5.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 +
                               9.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0);

        fints_yzz[i] +=
            fss * (-3.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * fz_0 + 14.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yzz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpb_y);

        fints_yzz[i] += ftt * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * rpa_z);

        fints_yzz[i] += ftt * ((9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] +=
            fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 - (9.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 -
                   (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 -
                   (9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * rpb_y * fz_0);

        fints_zzz[i] +=
            fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * rpb_y * fz_0 + 18.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y * fz_0 +
                   54.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * fz_0 -
                   (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * rpb_y * fz_0);

        fints_zzz[i] += fss * ((135.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y * fz_0 +
                               (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y * fz_0 + (45.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y * fz_0 +
                               15.0 * fe_0 * fe_0 * fe_0 * rpb_y * fz_0 - 3.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * fz_0);

        fints_zzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y * fz_0;

        fints_zzz[i] +=
            ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y +
                   (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y + (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_y +
                   (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpb_y);

        fints_zzz[i] += ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpb_y +
                               rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * rpb_y);
    }
}

auto
compPrimitiveKineticEnergyFG_T_ZZZZ(TDoubleArray&       buffer_xxx,
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-9.0 * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 - (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 -
                               3.0 * fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 - 9.0 * fe_0 * fke_0 * rpa_x * rpb_z * rpb_z * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxx[i] += fss * (18.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 +
                               36.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxx[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 - 6.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 +
                               14.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z + 3.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xxx[i] += ftt * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z;

        fints_xxy[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 - 3.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpa_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpb_z * rpb_z * fz_0);

        fints_xxy[i] += fss * (36.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 +
                               6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0);

        fints_xxy[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 - 6.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 +
                               14.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_xxy[i] += ftt * (3.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_xxy[i] += ftt * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z;

        fints_xxz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 - 2.0 * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 - 3.0 * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_xxz[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpa_x * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpb_z * rpb_z * fz_0 -
                               6.0 * fe_0 * fke_0 * rpa_x * rpa_x * rpb_z * fz_0 + 36.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 +
                               6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_xxz[i] += fss * (24.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 -
                               3.0 * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * fz_0 +
                               15.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0);

        fints_xxz[i] += fss * (30.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * fz_0 + 10.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 + 12.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               6.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * fz_0);

        fints_xxz[i] += fss * 14.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0;

        fints_xxz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z +
                               2.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z);

        fints_xxz[i] +=
            ftt * (3.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z + fe_0 * fe_0 * rpb_z * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                   (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z);

        fints_xyy[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 - 3.0 * fe_0 * fke_0 * rpa_y * rpa_y * rpa_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_x * rpb_z * rpb_z * fz_0);

        fints_xyy[i] += fss * (36.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 +
                               6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0 +
                               (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * fz_0 + 15.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0);

        fints_xyy[i] += fss * (3.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 - 6.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 +
                               14.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_xyy[i] += ftt * (3.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x);

        fints_xyy[i] += ftt * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z;

        fints_xyz[i] +=
            fss * (-3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpa_x * fz_0 - 6.0 * fe_0 * fke_0 * rpa_y * rpa_x * rpb_z * fz_0 +
                   36.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 + 24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 +
                   (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * fz_0);

        fints_xyz[i] += fss * (30.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * fz_0 - 6.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 +
                               14.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_xyz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z + 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x + 3.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z +
                               rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z);

        fints_xzz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_x * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 - 12.0 * fe_0 * fke_0 * rpa_z * rpa_x * rpb_z * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpa_x * fz_0);

        fints_xzz[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_x * rpb_z * rpb_z * fz_0 + 48.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 +
                               36.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * fz_0 +
                               6.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_x * fz_0);

        fints_xzz[i] += fss * (60.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z * fz_0 + 15.0 * fe_0 * fe_0 * fe_0 * rpa_x * fz_0 -
                               6.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * fz_0);

        fints_xzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z * fz_0;

        fints_xzz[i] += ftt * (4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z + 3.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z + 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x);

        fints_xzz[i] += ftt * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z * rpb_z + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_x +
                               rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpb_z);

        fints_yyy[i] += fss * (-9.0 * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 - (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 -
                               3.0 * fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 - 9.0 * fe_0 * fke_0 * rpa_y * rpb_z * rpb_z * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_y * rpa_y * fz_0);

        fints_yyy[i] += fss * (18.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 +
                               36.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * fz_0);

        fints_yyy[i] += fss * (9.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 - 6.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 +
                               14.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z + 3.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z +
                               (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y +
                               (9.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y);

        fints_yyy[i] += ftt * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z;

        fints_yyz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 - 2.0 * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 - 3.0 * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_yyz[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpa_y * fz_0 - 3.0 * fe_0 * fke_0 * rpa_z * rpb_z * rpb_z * fz_0 -
                               6.0 * fe_0 * fke_0 * rpa_y * rpa_y * rpb_z * fz_0 + 36.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 +
                               6.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_yyz[i] += fss * (24.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * fz_0 - (3.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 -
                               3.0 * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * fz_0 +
                               15.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0);

        fints_yyz[i] += fss * (30.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * fz_0 + 10.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 +
                               3.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 + 12.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               6.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * fz_0);

        fints_yyz[i] += fss * 14.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * fz_0;

        fints_yyz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z +
                               2.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z);

        fints_yyz[i] +=
            ftt * (3.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z + fe_0 * fe_0 * rpb_z * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z +
                   (3.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z + rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z);

        fints_yzz[i] += fss * (-3.0 * fbe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 - (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_y * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 - 12.0 * fe_0 * fke_0 * rpa_z * rpa_y * rpb_z * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpa_y * fz_0);

        fints_yzz[i] += fss * (-3.0 * fe_0 * fke_0 * rpa_y * rpb_z * rpb_z * fz_0 + 48.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * fz_0 +
                               36.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * fz_0 +
                               6.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 - (9.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_y * fz_0);

        fints_yzz[i] += fss * (60.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * fz_0 + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * fz_0 +
                               45.0 * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z * fz_0 + 15.0 * fe_0 * fe_0 * fe_0 * rpa_y * fz_0 -
                               6.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * fz_0);

        fints_yzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z * fz_0;

        fints_yzz[i] += ftt * (4.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z + 3.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z + 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y);

        fints_yzz[i] += ftt * ((9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z * rpb_z + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_y +
                               rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * (-9.0 * fbe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 - 6.0 * fbe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 -
                               (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * rpa_z * fz_0 - 9.0 * fbe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               3.0 * fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_zzz[i] += fss * (-9.0 * fe_0 * fke_0 * rpa_z * rpb_z * rpb_z * fz_0 - 18.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpb_z * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_z * rpa_z * fz_0 + 18.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z * fz_0 +
                               72.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * fz_0);

        fints_zzz[i] += fss * (36.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * fz_0 - (27.0 / 2.0) * fe_0 * fe_0 * fke_0 * rpa_z * fz_0 -
                               9.0 * fe_0 * fe_0 * fke_0 * rpb_z * fz_0 + 135.0 * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z * fz_0 +
                               90.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * fz_0);

        fints_zzz[i] += fss * ((15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * fz_0 + 30.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z * fz_0 +
                               45.0 * fe_0 * fe_0 * fe_0 * rpa_z * fz_0 + 60.0 * fe_0 * fe_0 * fe_0 * rpb_z * fz_0 -
                               6.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * fz_0);

        fints_zzz[i] += fss * 14.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z * fz_0;

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z + 6.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z +
                               3.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z + (27.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z * rpb_z +
                               9.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z);

        fints_zzz[i] += ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z + 3.0 * fe_0 * fe_0 * rpb_z * rpb_z * rpb_z +
                               (45.0 / 8.0) * fe_0 * fe_0 * fe_0 * rpa_z + (15.0 / 2.0) * fe_0 * fe_0 * fe_0 * rpb_z +
                               rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * rpb_z);
    }
}

}  // namespace kinrec
