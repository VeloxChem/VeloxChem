#include "OverlapGeom400RecDF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveOverlapGeom400DF_XX_XXX.hpp"
#include "PrimitiveOverlapGeom400DF_XX_XXY.hpp"
#include "PrimitiveOverlapGeom400DF_XX_XXZ.hpp"
#include "PrimitiveOverlapGeom400DF_XX_XYY.hpp"
#include "PrimitiveOverlapGeom400DF_XX_XYZ.hpp"
#include "PrimitiveOverlapGeom400DF_XX_XZZ.hpp"
#include "PrimitiveOverlapGeom400DF_XX_YYY.hpp"
#include "PrimitiveOverlapGeom400DF_XX_YYZ.hpp"
#include "PrimitiveOverlapGeom400DF_XX_YZZ.hpp"
#include "PrimitiveOverlapGeom400DF_XX_ZZZ.hpp"
#include "PrimitiveOverlapGeom400DF_XY_XXX.hpp"
#include "PrimitiveOverlapGeom400DF_XY_XXY.hpp"
#include "PrimitiveOverlapGeom400DF_XY_XXZ.hpp"
#include "PrimitiveOverlapGeom400DF_XY_XYY.hpp"
#include "PrimitiveOverlapGeom400DF_XY_XYZ.hpp"
#include "PrimitiveOverlapGeom400DF_XY_XZZ.hpp"
#include "PrimitiveOverlapGeom400DF_XY_YYY.hpp"
#include "PrimitiveOverlapGeom400DF_XY_YYZ.hpp"
#include "PrimitiveOverlapGeom400DF_XY_YZZ.hpp"
#include "PrimitiveOverlapGeom400DF_XY_ZZZ.hpp"
#include "PrimitiveOverlapGeom400DF_XZ_XXX.hpp"
#include "PrimitiveOverlapGeom400DF_XZ_XXY.hpp"
#include "PrimitiveOverlapGeom400DF_XZ_XXZ.hpp"
#include "PrimitiveOverlapGeom400DF_XZ_XYY.hpp"
#include "PrimitiveOverlapGeom400DF_XZ_XYZ.hpp"
#include "PrimitiveOverlapGeom400DF_XZ_XZZ.hpp"
#include "PrimitiveOverlapGeom400DF_XZ_YYY.hpp"
#include "PrimitiveOverlapGeom400DF_XZ_YYZ.hpp"
#include "PrimitiveOverlapGeom400DF_XZ_YZZ.hpp"
#include "PrimitiveOverlapGeom400DF_XZ_ZZZ.hpp"
#include "PrimitiveOverlapGeom400DF_YY_XXX.hpp"
#include "PrimitiveOverlapGeom400DF_YY_XXY.hpp"
#include "PrimitiveOverlapGeom400DF_YY_XXZ.hpp"
#include "PrimitiveOverlapGeom400DF_YY_XYY.hpp"
#include "PrimitiveOverlapGeom400DF_YY_XYZ.hpp"
#include "PrimitiveOverlapGeom400DF_YY_XZZ.hpp"
#include "PrimitiveOverlapGeom400DF_YY_YYY.hpp"
#include "PrimitiveOverlapGeom400DF_YY_YYZ.hpp"
#include "PrimitiveOverlapGeom400DF_YY_YZZ.hpp"
#include "PrimitiveOverlapGeom400DF_YY_ZZZ.hpp"
#include "PrimitiveOverlapGeom400DF_YZ_XXX.hpp"
#include "PrimitiveOverlapGeom400DF_YZ_XXY.hpp"
#include "PrimitiveOverlapGeom400DF_YZ_XXZ.hpp"
#include "PrimitiveOverlapGeom400DF_YZ_XYY.hpp"
#include "PrimitiveOverlapGeom400DF_YZ_XYZ.hpp"
#include "PrimitiveOverlapGeom400DF_YZ_XZZ.hpp"
#include "PrimitiveOverlapGeom400DF_YZ_YYY.hpp"
#include "PrimitiveOverlapGeom400DF_YZ_YYZ.hpp"
#include "PrimitiveOverlapGeom400DF_YZ_YZZ.hpp"
#include "PrimitiveOverlapGeom400DF_YZ_ZZZ.hpp"
#include "PrimitiveOverlapGeom400DF_ZZ_XXX.hpp"
#include "PrimitiveOverlapGeom400DF_ZZ_XXY.hpp"
#include "PrimitiveOverlapGeom400DF_ZZ_XXZ.hpp"
#include "PrimitiveOverlapGeom400DF_ZZ_XYY.hpp"
#include "PrimitiveOverlapGeom400DF_ZZ_XYZ.hpp"
#include "PrimitiveOverlapGeom400DF_ZZ_XZZ.hpp"
#include "PrimitiveOverlapGeom400DF_ZZ_YYY.hpp"
#include "PrimitiveOverlapGeom400DF_ZZ_YYZ.hpp"
#include "PrimitiveOverlapGeom400DF_ZZ_YZZ.hpp"
#include "PrimitiveOverlapGeom400DF_ZZ_ZZZ.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGeom400DF(CSubMatrix*      matrix_xxxx,
                     CSubMatrix*      matrix_xxxy,
                     CSubMatrix*      matrix_xxxz,
                     CSubMatrix*      matrix_xxyy,
                     CSubMatrix*      matrix_xxyz,
                     CSubMatrix*      matrix_xxzz,
                     CSubMatrix*      matrix_xyyy,
                     CSubMatrix*      matrix_xyyz,
                     CSubMatrix*      matrix_xyzz,
                     CSubMatrix*      matrix_xzzz,
                     CSubMatrix*      matrix_yyyy,
                     CSubMatrix*      matrix_yyyz,
                     CSubMatrix*      matrix_yyzz,
                     CSubMatrix*      matrix_yzzz,
                     CSubMatrix*      matrix_zzzz,
                     const CGtoBlock& bra_gto_block,
                     const CGtoBlock& ket_gto_block,
                     const bool       ang_order,
                     const int64_t    bra_first,
                     const int64_t    bra_last) -> void
{
    // spherical transformation factors

    const double f2_3 = 2.0 * std::sqrt(3.0);

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

    alignas(64) TDoubleArray buffer_xxxx;

    alignas(64) TDoubleArray buffer_xxxy;

    alignas(64) TDoubleArray buffer_xxxz;

    alignas(64) TDoubleArray buffer_xxyy;

    alignas(64) TDoubleArray buffer_xxyz;

    alignas(64) TDoubleArray buffer_xxzz;

    alignas(64) TDoubleArray buffer_xyyy;

    alignas(64) TDoubleArray buffer_xyyz;

    alignas(64) TDoubleArray buffer_xyzz;

    alignas(64) TDoubleArray buffer_xzzz;

    alignas(64) TDoubleArray buffer_yyyy;

    alignas(64) TDoubleArray buffer_yyyz;

    alignas(64) TDoubleArray buffer_yyzz;

    alignas(64) TDoubleArray buffer_yzzz;

    alignas(64) TDoubleArray buffer_zzzz;

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

            // compute primitive integrals block (XX_XXX)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XX_XXX(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XXY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XX_XXY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XXZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XX_XXZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XX_XYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XX_XYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_XZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XX_XZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_YYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XX_YYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_YYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XX_YYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_YZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XX_YZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XX_ZZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XX_ZZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XXX)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XY_XXX(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XXY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XY_XXY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XXZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XY_XXZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XY_XYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XY_XYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_XZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XY_XZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_YYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XY_YYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_YYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XY_YYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_YZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XY_YZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XY_ZZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XY_ZZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XXX)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XZ_XXX(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XXY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XZ_XXY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XXZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XZ_XXZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XZ_XYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XZ_XYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_XZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XZ_XZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_YYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XZ_YYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_YYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XZ_YYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_YZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XZ_YZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZ_ZZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_XZ_ZZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XXX)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YY_XXX(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XXY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YY_XXY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XXZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YY_XXZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YY_XYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YY_XYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_XZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YY_XZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_YYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YY_YYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_YYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YY_YYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_YZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YY_YZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YY_ZZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YY_ZZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XXX)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YZ_XXX(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XXY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YZ_XXY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XXZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YZ_XXZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YZ_XYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f2_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YZ_XYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_XZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YZ_XZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_YYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YZ_YYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_YYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YZ_YYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f2_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_YZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YZ_YZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZ_ZZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_YZ_ZZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XXX)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_ZZ_XXX(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XXY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_ZZ_XXY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XXZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_ZZ_XXZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_ZZ_XYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_ZZ_XYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_XZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_ZZ_XZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_YYY)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_ZZ_YYY(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_YYZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_ZZ_YYZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_YZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_ZZ_YZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZ_ZZZ)

            simd::zero(buffer_xxxx);

            simd::zero(buffer_xxxy);

            simd::zero(buffer_xxxz);

            simd::zero(buffer_xxyy);

            simd::zero(buffer_xxyz);

            simd::zero(buffer_xxzz);

            simd::zero(buffer_xyyy);

            simd::zero(buffer_xyyz);

            simd::zero(buffer_xyzz);

            simd::zero(buffer_xzzz);

            simd::zero(buffer_yyyy);

            simd::zero(buffer_yyyz);

            simd::zero(buffer_yyzz);

            simd::zero(buffer_yzzz);

            simd::zero(buffer_zzzz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom400DF_ZZ_ZZZ(buffer_xxxx,
                                                                 buffer_xxxy,
                                                                 buffer_xxxz,
                                                                 buffer_xxyy,
                                                                 buffer_xxyz,
                                                                 buffer_xxzz,
                                                                 buffer_xyyy,
                                                                 buffer_xyyz,
                                                                 buffer_xyzz,
                                                                 buffer_xzzz,
                                                                 buffer_yyyy,
                                                                 buffer_yyyz,
                                                                 buffer_yyzz,
                                                                 buffer_yzzz,
                                                                 buffer_zzzz,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace ovlrec
