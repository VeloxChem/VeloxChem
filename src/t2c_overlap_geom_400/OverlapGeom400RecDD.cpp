#include "OverlapGeom400RecDD.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveOverlapGeom400DD_XX_XX.hpp"
#include "PrimitiveOverlapGeom400DD_XX_XY.hpp"
#include "PrimitiveOverlapGeom400DD_XX_XZ.hpp"
#include "PrimitiveOverlapGeom400DD_XX_YY.hpp"
#include "PrimitiveOverlapGeom400DD_XX_YZ.hpp"
#include "PrimitiveOverlapGeom400DD_XX_ZZ.hpp"
#include "PrimitiveOverlapGeom400DD_XY_XX.hpp"
#include "PrimitiveOverlapGeom400DD_XY_XY.hpp"
#include "PrimitiveOverlapGeom400DD_XY_XZ.hpp"
#include "PrimitiveOverlapGeom400DD_XY_YY.hpp"
#include "PrimitiveOverlapGeom400DD_XY_YZ.hpp"
#include "PrimitiveOverlapGeom400DD_XY_ZZ.hpp"
#include "PrimitiveOverlapGeom400DD_XZ_XX.hpp"
#include "PrimitiveOverlapGeom400DD_XZ_XY.hpp"
#include "PrimitiveOverlapGeom400DD_XZ_XZ.hpp"
#include "PrimitiveOverlapGeom400DD_XZ_YY.hpp"
#include "PrimitiveOverlapGeom400DD_XZ_YZ.hpp"
#include "PrimitiveOverlapGeom400DD_XZ_ZZ.hpp"
#include "PrimitiveOverlapGeom400DD_YY_XX.hpp"
#include "PrimitiveOverlapGeom400DD_YY_XY.hpp"
#include "PrimitiveOverlapGeom400DD_YY_XZ.hpp"
#include "PrimitiveOverlapGeom400DD_YY_YY.hpp"
#include "PrimitiveOverlapGeom400DD_YY_YZ.hpp"
#include "PrimitiveOverlapGeom400DD_YY_ZZ.hpp"
#include "PrimitiveOverlapGeom400DD_YZ_XX.hpp"
#include "PrimitiveOverlapGeom400DD_YZ_XY.hpp"
#include "PrimitiveOverlapGeom400DD_YZ_XZ.hpp"
#include "PrimitiveOverlapGeom400DD_YZ_YY.hpp"
#include "PrimitiveOverlapGeom400DD_YZ_YZ.hpp"
#include "PrimitiveOverlapGeom400DD_YZ_ZZ.hpp"
#include "PrimitiveOverlapGeom400DD_ZZ_XX.hpp"
#include "PrimitiveOverlapGeom400DD_ZZ_XY.hpp"
#include "PrimitiveOverlapGeom400DD_ZZ_XZ.hpp"
#include "PrimitiveOverlapGeom400DD_ZZ_YY.hpp"
#include "PrimitiveOverlapGeom400DD_ZZ_YZ.hpp"
#include "PrimitiveOverlapGeom400DD_ZZ_ZZ.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGeom400DD(CSubMatrix*      matrix_xxxx,
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

    const double f2_3 = 2.0 * std::sqrt(3.0);

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

            // compute primitive integrals block (XX_XX)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XX_XX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XX_XY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XX_XY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XX_XZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XX_XZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XX_YY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XX_YY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XX_YZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XX_YZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XX_ZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XX_ZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_XX)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XY_XX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_XY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XY_XY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_XZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XY_XZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_YY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XY_YY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_YZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XY_YZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_ZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XY_ZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_XX)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XZ_XX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_XY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XZ_XY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_XZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XZ_XZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_YY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XZ_YY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_YZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XZ_YZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_ZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_XZ_ZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_XX)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YY_XX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_XY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YY_XY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_XZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YY_XZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_YY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YY_YY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_YZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YY_YZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_ZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YY_ZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_XX)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YZ_XX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_XY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YZ_XY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_XZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YZ_XZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_YY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YZ_YY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_YZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YZ_YZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_ZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_YZ_ZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_XX)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_ZZ_XX(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_XY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_ZZ_XY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_XZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_ZZ_XZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_YY)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_ZZ_YY(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxx, buffer_xxxx, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxy, buffer_xxxy, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxxz, buffer_xxxz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyy, buffer_xxyy, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxyz, buffer_xxyz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxzz, buffer_xxzz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyy, buffer_xyyy, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyyz, buffer_xyyz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyzz, buffer_xyzz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzzz, buffer_xzzz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyy, buffer_yyyy, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyyz, buffer_yyyz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyzz, buffer_yyzz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzzz, buffer_yzzz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzzz, buffer_zzzz, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_YZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_ZZ_YZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_ZZ)

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

                    ovlrec::compPrimitiveOverlapGeom400DD_ZZ_ZZ(buffer_xxxx,
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

            t2cfunc::distribute(matrix_xxxx, buffer_xxxx, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxy, buffer_xxxy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxxz, buffer_xxxz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyy, buffer_xxyy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxyz, buffer_xxyz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxzz, buffer_xxzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyy, buffer_xyyy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyyz, buffer_xyyz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyzz, buffer_xyzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzzz, buffer_xzzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyy, buffer_yyyy, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyyz, buffer_yyyz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyzz, buffer_yyzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzzz, buffer_yzzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzzz, buffer_zzzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);
        }
    }
}

}  // namespace ovlrec
