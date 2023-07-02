#include "NuclearPotentialRecFG.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveNuclearPotentialFG_T_XXXX.hpp"
#include "PrimitiveNuclearPotentialFG_T_XXXY.hpp"
#include "PrimitiveNuclearPotentialFG_T_XXXZ.hpp"
#include "PrimitiveNuclearPotentialFG_T_XXYY.hpp"
#include "PrimitiveNuclearPotentialFG_T_XXYZ.hpp"
#include "PrimitiveNuclearPotentialFG_T_XXZZ.hpp"
#include "PrimitiveNuclearPotentialFG_T_XYYY.hpp"
#include "PrimitiveNuclearPotentialFG_T_XYYZ.hpp"
#include "PrimitiveNuclearPotentialFG_T_XYZZ.hpp"
#include "PrimitiveNuclearPotentialFG_T_XZZZ.hpp"
#include "PrimitiveNuclearPotentialFG_T_YYYY.hpp"
#include "PrimitiveNuclearPotentialFG_T_YYYZ.hpp"
#include "PrimitiveNuclearPotentialFG_T_YYZZ.hpp"
#include "PrimitiveNuclearPotentialFG_T_YZZZ.hpp"
#include "PrimitiveNuclearPotentialFG_T_ZZZZ.hpp"
#include "T2CDistributor.hpp"

namespace npotrec {  // npotrec namespace

auto
compNuclearPotentialFG(CSubMatrix*      matrix,
                       const double     charge,
                       const TPoint3D&  point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_XXXX(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_XXXY(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_XXXZ(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_XXYY(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_XXYZ(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_XXZZ(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_XYYY(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_XYYZ(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_XYZZ(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_XZZZ(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_YYYY(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_YYYZ(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_YYZZ(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_YZZZ(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

                    npotrec::compPrimitiveNuclearPotentialFG_T_ZZZZ(buffer_xxx,
                                                                    buffer_xxy,
                                                                    buffer_xxz,
                                                                    buffer_xyy,
                                                                    buffer_xyz,
                                                                    buffer_xzz,
                                                                    buffer_yyy,
                                                                    buffer_yyz,
                                                                    buffer_yzz,
                                                                    buffer_zzz,
                                                                    charge,
                                                                    point,
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

}  // namespace npotrec
