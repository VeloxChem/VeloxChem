#include "OverlapGeom400RecFF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveOverlapGeom400FF_XXX_XXX.hpp"
#include "PrimitiveOverlapGeom400FF_XXX_XXY.hpp"
#include "PrimitiveOverlapGeom400FF_XXX_XXZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXX_XYY.hpp"
#include "PrimitiveOverlapGeom400FF_XXX_XYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXX_XZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXX_YYY.hpp"
#include "PrimitiveOverlapGeom400FF_XXX_YYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXX_YZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXX_ZZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXY_XXX.hpp"
#include "PrimitiveOverlapGeom400FF_XXY_XXY.hpp"
#include "PrimitiveOverlapGeom400FF_XXY_XXZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXY_XYY.hpp"
#include "PrimitiveOverlapGeom400FF_XXY_XYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXY_XZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXY_YYY.hpp"
#include "PrimitiveOverlapGeom400FF_XXY_YYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXY_YZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXY_ZZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXZ_XXX.hpp"
#include "PrimitiveOverlapGeom400FF_XXZ_XXY.hpp"
#include "PrimitiveOverlapGeom400FF_XXZ_XXZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXZ_XYY.hpp"
#include "PrimitiveOverlapGeom400FF_XXZ_XYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXZ_XZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXZ_YYY.hpp"
#include "PrimitiveOverlapGeom400FF_XXZ_YYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXZ_YZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XXZ_ZZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYY_XXX.hpp"
#include "PrimitiveOverlapGeom400FF_XYY_XXY.hpp"
#include "PrimitiveOverlapGeom400FF_XYY_XXZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYY_XYY.hpp"
#include "PrimitiveOverlapGeom400FF_XYY_XYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYY_XZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYY_YYY.hpp"
#include "PrimitiveOverlapGeom400FF_XYY_YYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYY_YZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYY_ZZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYZ_XXX.hpp"
#include "PrimitiveOverlapGeom400FF_XYZ_XXY.hpp"
#include "PrimitiveOverlapGeom400FF_XYZ_XXZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYZ_XYY.hpp"
#include "PrimitiveOverlapGeom400FF_XYZ_XYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYZ_XZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYZ_YYY.hpp"
#include "PrimitiveOverlapGeom400FF_XYZ_YYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYZ_YZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XYZ_ZZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XZZ_XXX.hpp"
#include "PrimitiveOverlapGeom400FF_XZZ_XXY.hpp"
#include "PrimitiveOverlapGeom400FF_XZZ_XXZ.hpp"
#include "PrimitiveOverlapGeom400FF_XZZ_XYY.hpp"
#include "PrimitiveOverlapGeom400FF_XZZ_XYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XZZ_XZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XZZ_YYY.hpp"
#include "PrimitiveOverlapGeom400FF_XZZ_YYZ.hpp"
#include "PrimitiveOverlapGeom400FF_XZZ_YZZ.hpp"
#include "PrimitiveOverlapGeom400FF_XZZ_ZZZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYY_XXX.hpp"
#include "PrimitiveOverlapGeom400FF_YYY_XXY.hpp"
#include "PrimitiveOverlapGeom400FF_YYY_XXZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYY_XYY.hpp"
#include "PrimitiveOverlapGeom400FF_YYY_XYZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYY_XZZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYY_YYY.hpp"
#include "PrimitiveOverlapGeom400FF_YYY_YYZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYY_YZZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYY_ZZZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYZ_XXX.hpp"
#include "PrimitiveOverlapGeom400FF_YYZ_XXY.hpp"
#include "PrimitiveOverlapGeom400FF_YYZ_XXZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYZ_XYY.hpp"
#include "PrimitiveOverlapGeom400FF_YYZ_XYZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYZ_XZZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYZ_YYY.hpp"
#include "PrimitiveOverlapGeom400FF_YYZ_YYZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYZ_YZZ.hpp"
#include "PrimitiveOverlapGeom400FF_YYZ_ZZZ.hpp"
#include "PrimitiveOverlapGeom400FF_YZZ_XXX.hpp"
#include "PrimitiveOverlapGeom400FF_YZZ_XXY.hpp"
#include "PrimitiveOverlapGeom400FF_YZZ_XXZ.hpp"
#include "PrimitiveOverlapGeom400FF_YZZ_XYY.hpp"
#include "PrimitiveOverlapGeom400FF_YZZ_XYZ.hpp"
#include "PrimitiveOverlapGeom400FF_YZZ_XZZ.hpp"
#include "PrimitiveOverlapGeom400FF_YZZ_YYY.hpp"
#include "PrimitiveOverlapGeom400FF_YZZ_YYZ.hpp"
#include "PrimitiveOverlapGeom400FF_YZZ_YZZ.hpp"
#include "PrimitiveOverlapGeom400FF_YZZ_ZZZ.hpp"
#include "PrimitiveOverlapGeom400FF_ZZZ_XXX.hpp"
#include "PrimitiveOverlapGeom400FF_ZZZ_XXY.hpp"
#include "PrimitiveOverlapGeom400FF_ZZZ_XXZ.hpp"
#include "PrimitiveOverlapGeom400FF_ZZZ_XYY.hpp"
#include "PrimitiveOverlapGeom400FF_ZZZ_XYZ.hpp"
#include "PrimitiveOverlapGeom400FF_ZZZ_XZZ.hpp"
#include "PrimitiveOverlapGeom400FF_ZZZ_YYY.hpp"
#include "PrimitiveOverlapGeom400FF_ZZZ_YYZ.hpp"
#include "PrimitiveOverlapGeom400FF_ZZZ_YZZ.hpp"
#include "PrimitiveOverlapGeom400FF_ZZZ_ZZZ.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGeom400FF(CSubMatrix*      matrix_xxxx,
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
                     const int64_t    bra_first,
                     const int64_t    bra_last,
                     const mat_t      mat_type) -> void
{
    // spherical transformation factors

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

            // compute primitive integrals block (XXX_XXX)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXX_XXX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_XXY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXX_XXY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_XXZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXX_XXZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_XYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXX_XYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_XYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXX_XYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_XZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXX_XZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_YYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXX_YYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_YYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXX_YYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_YZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXX_YZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_ZZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXX_ZZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XXX)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXY_XXX(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XXY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXY_XXY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XXZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXY_XXZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXY_XYY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXY_XYZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXY_XZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_YYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXY_YYY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_YYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXY_YYZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_YZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXY_YZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_ZZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXY_ZZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XXX)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXZ_XXX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XXY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXZ_XXY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XXZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXZ_XXZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXZ_XYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXZ_XYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXZ_XZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_YYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXZ_YYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_YYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXZ_YYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_YZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXZ_YZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_ZZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XXZ_ZZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XXX)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYY_XXX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XXY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYY_XXY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XXZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYY_XXZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYY_XYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYY_XYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYY_XZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_YYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYY_YYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_YYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYY_YYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_YZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYY_YZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_ZZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYY_ZZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XXX)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYZ_XXX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XXY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYZ_XXY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XXZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYZ_XXZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYZ_XYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYZ_XYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYZ_XZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_YYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYZ_YYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_YYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYZ_YYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_YZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYZ_YZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_ZZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XYZ_ZZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XXX)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XZZ_XXX(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XXY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XZZ_XXY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XXZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XZZ_XXZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XZZ_XYY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XZZ_XYZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XZZ_XZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_YYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XZZ_YYY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_YYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XZZ_YYZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_YZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XZZ_YZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_ZZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_XZZ_ZZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XXX)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYY_XXX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XXY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYY_XXY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XXZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYY_XXZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYY_XYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYY_XYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYY_XZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_YYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYY_YYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_YYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYY_YYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_YZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYY_YZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_ZZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYY_ZZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XXX)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYZ_XXX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XXY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYZ_XXY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XXZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYZ_XXZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYZ_XYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYZ_XYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYZ_XZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_YYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYZ_YYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_YYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYZ_YYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_YZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYZ_YZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_ZZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YYZ_ZZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XXX)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YZZ_XXX(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XXY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YZZ_XXY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XXZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YZZ_XXZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YZZ_XYY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YZZ_XYZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YZZ_XZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_YYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YZZ_YYY(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_YYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YZZ_YYZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_YZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YZZ_YZZ(buffer_xxxx,
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
                matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_ZZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_YZZ_ZZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XXX)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_ZZZ_XXX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XXY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_ZZZ_XXY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XXZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_ZZZ_XXZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_ZZZ_XYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_ZZZ_XYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_ZZZ_XZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_YYY)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_ZZZ_YYY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_YYZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_ZZZ_YYZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_YZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_ZZZ_YZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_ZZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400FF_ZZZ_ZZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);
        }
    }
}

}  // namespace ovlrec
