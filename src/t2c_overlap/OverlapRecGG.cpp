#include "OverlapRecGG.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveOverlapGG_XXXX_T.hpp"
#include "PrimitiveOverlapGG_XXXY_T.hpp"
#include "PrimitiveOverlapGG_XXXZ_T.hpp"
#include "PrimitiveOverlapGG_XXYY_T.hpp"
#include "PrimitiveOverlapGG_XXYZ_T.hpp"
#include "PrimitiveOverlapGG_XXZZ_T.hpp"
#include "PrimitiveOverlapGG_XYYY_T.hpp"
#include "PrimitiveOverlapGG_XYYZ_T.hpp"
#include "PrimitiveOverlapGG_XYZZ_T.hpp"
#include "PrimitiveOverlapGG_XZZZ_T.hpp"
#include "PrimitiveOverlapGG_YYYY_T.hpp"
#include "PrimitiveOverlapGG_YYYZ_T.hpp"
#include "PrimitiveOverlapGG_YYZZ_T.hpp"
#include "PrimitiveOverlapGG_YZZZ_T.hpp"
#include "PrimitiveOverlapGG_ZZZZ_T.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGG(CSubMatrix* matrix, const CGtoBlock& gto_block, const int64_t bra_first, const int64_t bra_last) -> void
{
    // spherical transformation factors

    const double f4_35 = 4.0 * std::sqrt(35);

    const double f4_17 = 4.0 * std::sqrt(17.5);

    const double f4_5 = 4.0 * std::sqrt(5.0);

    const double f4_2 = 4.0 * std::sqrt(2.5);

    // intialize GTOs data

    const auto gto_coords = gto_block.getCoordinates();

    const auto gto_exps = gto_block.getExponents();

    const auto gto_norms = gto_block.getNormalizationFactors();

    const auto gto_indexes = gto_block.getOrbitalIndexes();

    const auto ncgtos = gto_block.getNumberOfBasisFunctions();

    const auto npgtos = gto_block.getNumberOfPrimitives();

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

    const auto nbatches = batch::getNumberOfBatches(ncgtos, simd_width);

    for (int64_t i = 0; i < nbatches; i++)
    {
        const auto [ket_first, ket_last] = batch::getBatchRange(i, ncgtos, simd_width);

        const auto ket_dim = ket_last - ket_first;

        simd::loadCoordinates(ket_coords_x, ket_coords_y, ket_coords_z, gto_coords, ket_first, ket_last);

        for (int64_t j = bra_first; j < bra_last; j++)
        {
            const auto bra_coord = gto_coords[j];

            // compute primitive integrals block (XXXX)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_XXXX_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -0.5 * f4_5 * 3.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.25 * f4_35 * 3.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.5 * f4_5 * 0.5 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -0.25 * f4_35 * 0.5 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -0.5 * f4_5 * 0.25 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.25 * f4_35 * 0.25 * f4_35, gto_indexes, 8, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -0.5 * f4_5 * f4_35, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 0.25 * f4_35 * f4_35, gto_indexes, 8, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 0.5 * f4_5 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -0.25 * f4_35 * f4_5, gto_indexes, 8, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 0.5 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -0.25 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -0.5 * f4_5 * f4_17, gto_indexes, 6, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 0.25 * f4_35 * f4_17, gto_indexes, 8, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * 6.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -0.5 * f4_5 * 6.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 0.25 * f4_35 * 6.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * 1.50 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 0.5 * f4_5 * 1.50 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -0.25 * f4_35 * 1.50 * f4_35, gto_indexes, 8, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * 3.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -0.5 * f4_5 * 3.0 * f4_17, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 0.25 * f4_35 * 3.0 * f4_17, gto_indexes, 8, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 0.5 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -0.25 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 0.5 * f4_5 * 24.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -0.25 * f4_35 * 24.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -0.5 * f4_5 * 3.0 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 0.25 * f4_35 * 3.0 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 0.5 * f4_5 * f4_35, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.25 * f4_35 * f4_35, gto_indexes, 8, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 0.5 * f4_5 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.25 * f4_35 * f4_5, gto_indexes, 8, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 0.5 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -0.25 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * 3.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 0.5 * f4_5 * 3.0 * f4_17, gto_indexes, 6, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -0.25 * f4_35 * 3.0 * f4_17, gto_indexes, 8, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 3.0 * 6.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -0.5 * f4_5 * 6.0 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 0.25 * f4_35 * 6.0 * f4_5, gto_indexes, 8, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 3.0 * 4.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -0.5 * f4_5 * 4.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 0.25 * f4_35 * 4.0 * f4_2, gto_indexes, 8, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -0.5 * f4_5 * 3.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.25 * f4_35 * 3.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -0.5 * f4_5 * 0.5 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.25 * f4_35 * 0.5 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -0.5 * f4_5 * 0.25 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.25 * f4_35 * 0.25 * f4_35, gto_indexes, 8, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 0.5 * f4_5 * f4_17, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -0.25 * f4_35 * f4_17, gto_indexes, 8, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 0.5 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -0.25 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 0.5 * f4_5 * 24.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -0.25 * f4_35 * 24.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 0.5 * f4_5 * 3.0 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -0.25 * f4_35 * 3.0 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 3.0 * 4.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -0.5 * f4_5 * 4.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 0.25 * f4_35 * 4.0 * f4_2, gto_indexes, 8, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 3.0 * 8.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -0.5 * f4_5 * 8.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 0.25 * f4_35 * 8.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XXXY)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_XXXY_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, f4_35 * 3.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_5 * 3.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_35 * 0.5 * f4_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_5 * 0.5 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_35 * 0.25 * f4_35, gto_indexes, 0, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_5 * 0.25 * f4_35, gto_indexes, 2, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_35 * f4_35, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_5 * f4_35, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_35 * f4_5, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_5 * f4_5, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_35 * 3.0 * f4_2, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_5 * 3.0 * f4_2, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_35 * f4_17, gto_indexes, 0, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_5 * f4_17, gto_indexes, 2, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_35 * 6.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_5 * 6.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_35 * 1.50 * f4_35, gto_indexes, 0, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_5 * 1.50 * f4_35, gto_indexes, 2, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_35 * 3.0 * f4_17, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_5 * 3.0 * f4_17, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_35 * 3.0 * f4_2, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_5 * 3.0 * f4_2, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_35 * 24.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_5 * 24.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_35 * 3.0 * f4_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_5 * 3.0 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_35 * f4_35, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_5 * f4_35, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_35 * f4_5, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_5 * f4_5, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -f4_35 * 3.0 * f4_2, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_5 * 3.0 * f4_2, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -f4_35 * 3.0 * f4_17, gto_indexes, 0, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_5 * 3.0 * f4_17, gto_indexes, 2, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, f4_35 * 6.0 * f4_5, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -f4_5 * 6.0 * f4_5, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, f4_35 * 4.0 * f4_2, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -f4_5 * 4.0 * f4_2, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_35 * 3.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 3.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_35 * 0.5 * f4_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 0.5 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_35 * 0.25 * f4_35, gto_indexes, 0, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 0.25 * f4_35, gto_indexes, 2, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -f4_35 * f4_17, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_5 * f4_17, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -f4_35 * 3.0 * f4_2, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_5 * 3.0 * f4_2, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -f4_35 * 24.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_5 * 24.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -f4_35 * 3.0 * f4_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_5 * 3.0 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, f4_35 * 4.0 * f4_2, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -f4_5 * 4.0 * f4_2, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, f4_35 * 8.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -f4_5 * 8.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XXXZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_XXXZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 3.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_17 * 3.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_2 * 0.5 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_17 * 0.5 * f4_5, gto_indexes, 7, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 0.25 * f4_35, gto_indexes, 5, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_17 * 0.25 * f4_35, gto_indexes, 7, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_2 * f4_35, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_17 * f4_35, gto_indexes, 7, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_2 * f4_5, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_17 * f4_5, gto_indexes, 7, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_17 * 3.0 * f4_2, gto_indexes, 7, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_2 * f4_17, gto_indexes, 5, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_17 * f4_17, gto_indexes, 7, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_2 * 6.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_17 * 6.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * f4_2 * 1.50 * f4_35, gto_indexes, 5, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_17 * 1.50 * f4_35, gto_indexes, 7, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * f4_2 * 3.0 * f4_17, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_17 * 3.0 * f4_17, gto_indexes, 7, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_17 * 3.0 * f4_2, gto_indexes, 7, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_2 * 24.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_17 * 24.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * f4_2 * 3.0 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_17 * 3.0 * f4_5, gto_indexes, 7, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_35, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_17 * f4_35, gto_indexes, 7, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_5, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_17 * f4_5, gto_indexes, 7, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -f4_17 * 3.0 * f4_2, gto_indexes, 7, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_17, gto_indexes, 5, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -f4_17 * 3.0 * f4_17, gto_indexes, 7, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -3.0 * f4_2 * 6.0 * f4_5, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, f4_17 * 6.0 * f4_5, gto_indexes, 7, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -3.0 * f4_2 * 4.0 * f4_2, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, f4_17 * 4.0 * f4_2, gto_indexes, 7, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 3.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_17 * 3.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 0.5 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_17 * 0.5 * f4_5, gto_indexes, 7, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 0.25 * f4_35, gto_indexes, 5, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_17 * 0.25 * f4_35, gto_indexes, 7, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * f4_17, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -f4_17 * f4_17, gto_indexes, 7, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -f4_17 * 3.0 * f4_2, gto_indexes, 7, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 24.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -f4_17 * 24.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 3.0 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -f4_17 * 3.0 * f4_5, gto_indexes, 7, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -3.0 * f4_2 * 4.0 * f4_2, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, f4_17 * 4.0 * f4_2, gto_indexes, 7, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_2 * 8.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, f4_17 * 8.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XXYY)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_XXYY_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 6.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -1.50 * f4_35 * 3.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -6.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 1.50 * f4_35 * 0.5 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 6.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -1.50 * f4_35 * 0.25 * f4_35, gto_indexes, 8, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 6.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -1.50 * f4_35 * f4_35, gto_indexes, 8, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -6.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 1.50 * f4_35 * f4_5, gto_indexes, 8, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -6.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 1.50 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 6.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -1.50 * f4_35 * f4_17, gto_indexes, 8, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 6.0 * 6.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -1.50 * f4_35 * 6.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -6.0 * 1.50 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 1.50 * f4_35 * 1.50 * f4_35, gto_indexes, 8, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 6.0 * 3.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -1.50 * f4_35 * 3.0 * f4_17, gto_indexes, 8, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -6.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 1.50 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -6.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 1.50 * f4_35 * 24.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 6.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -1.50 * f4_35 * 3.0 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -6.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 1.50 * f4_35 * f4_35, gto_indexes, 8, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -6.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 1.50 * f4_35 * f4_5, gto_indexes, 8, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -6.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 1.50 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -6.0 * 3.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 1.50 * f4_35 * 3.0 * f4_17, gto_indexes, 8, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 6.0 * 6.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -1.50 * f4_35 * 6.0 * f4_5, gto_indexes, 8, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 6.0 * 4.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -1.50 * f4_35 * 4.0 * f4_2, gto_indexes, 8, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 6.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -1.50 * f4_35 * 3.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 6.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -1.50 * f4_35 * 0.5 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 6.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -1.50 * f4_35 * 0.25 * f4_35, gto_indexes, 8, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -6.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 1.50 * f4_35 * f4_17, gto_indexes, 8, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -6.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 1.50 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -6.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 1.50 * f4_35 * 24.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -6.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 1.50 * f4_35 * 3.0 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 6.0 * 4.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -1.50 * f4_35 * 4.0 * f4_2, gto_indexes, 8, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 6.0 * 8.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -1.50 * f4_35 * 8.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XXYZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_XXYZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_17 * 3.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 3.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_17 * 0.5 * f4_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_2 * 0.5 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_17 * 0.25 * f4_35, gto_indexes, 1, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 0.25 * f4_35, gto_indexes, 3, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_17 * f4_35, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_2 * f4_35, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_17 * f4_5, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_2 * f4_5, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_17 * 3.0 * f4_2, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_17 * f4_17, gto_indexes, 1, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_2 * f4_17, gto_indexes, 3, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * f4_17 * 6.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_2 * 6.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_17 * 1.50 * f4_35, gto_indexes, 1, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * f4_2 * 1.50 * f4_35, gto_indexes, 3, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * f4_17 * 3.0 * f4_17, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * f4_2 * 3.0 * f4_17, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * f4_17 * 3.0 * f4_2, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * f4_17 * 24.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_2 * 24.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_17 * 3.0 * f4_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * f4_2 * 3.0 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_17 * f4_35, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_35, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_17 * f4_5, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_5, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * f4_17 * 3.0 * f4_2, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * f4_17 * 3.0 * f4_17, gto_indexes, 1, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_17, gto_indexes, 3, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 3.0 * f4_17 * 6.0 * f4_5, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -3.0 * f4_2 * 6.0 * f4_5, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 3.0 * f4_17 * 4.0 * f4_2, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -3.0 * f4_2 * 4.0 * f4_2, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * f4_17 * 3.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 3.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * f4_17 * 0.5 * f4_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 0.5 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * f4_17 * 0.25 * f4_35, gto_indexes, 1, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 0.25 * f4_35, gto_indexes, 3, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_17 * f4_17, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * f4_17, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_17 * 3.0 * f4_2, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * f4_17 * 24.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 24.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * f4_17 * 3.0 * f4_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 3.0 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 3.0 * f4_17 * 4.0 * f4_2, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -3.0 * f4_2 * 4.0 * f4_2, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 3.0 * f4_17 * 8.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_2 * 8.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XXZZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_XXZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -24.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_5 * 3.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 24.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_5 * 0.5 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -24.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_5 * 0.25 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -24.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_5 * f4_35, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 24.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_5 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 24.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -24.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_5 * f4_17, gto_indexes, 6, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -24.0 * 6.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * f4_5 * 6.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 24.0 * 1.50 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_5 * 1.50 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -24.0 * 3.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * f4_5 * 3.0 * f4_17, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 24.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 24.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * f4_5 * 24.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -24.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_5 * 3.0 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 24.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_5 * f4_35, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 24.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_5 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 24.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 24.0 * 3.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * f4_5 * 3.0 * f4_17, gto_indexes, 6, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -24.0 * 6.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 3.0 * f4_5 * 6.0 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -24.0 * 4.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 3.0 * f4_5 * 4.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * f4_5 * 3.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * f4_5 * 0.5 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * f4_5 * 0.25 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 24.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_5 * f4_17, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 24.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 24.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * f4_5 * 24.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 24.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * f4_5 * 3.0 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -24.0 * 4.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 3.0 * f4_5 * 4.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -24.0 * 8.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 3.0 * f4_5 * 8.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XYYY)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_XYYY_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_35 * 3.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_5 * 3.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_35 * 0.5 * f4_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_5 * 0.5 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_35 * 0.25 * f4_35, gto_indexes, 0, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_5 * 0.25 * f4_35, gto_indexes, 2, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_35 * f4_35, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_5 * f4_35, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_35 * f4_5, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_5 * f4_5, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_35 * 3.0 * f4_2, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_5 * 3.0 * f4_2, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_35 * f4_17, gto_indexes, 0, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_5 * f4_17, gto_indexes, 2, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_35 * 6.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_5 * 6.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_35 * 1.50 * f4_35, gto_indexes, 0, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_5 * 1.50 * f4_35, gto_indexes, 2, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_35 * 3.0 * f4_17, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_5 * 3.0 * f4_17, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_35 * 3.0 * f4_2, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_5 * 3.0 * f4_2, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_35 * 24.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_5 * 24.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_35 * 3.0 * f4_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_5 * 3.0 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_35 * f4_35, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_5 * f4_35, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_35 * f4_5, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_5 * f4_5, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_35 * 3.0 * f4_2, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_5 * 3.0 * f4_2, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_35 * 3.0 * f4_17, gto_indexes, 0, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_5 * 3.0 * f4_17, gto_indexes, 2, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -f4_35 * 6.0 * f4_5, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -f4_5 * 6.0 * f4_5, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -f4_35 * 4.0 * f4_2, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -f4_5 * 4.0 * f4_2, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_35 * 3.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 3.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_35 * 0.5 * f4_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 0.5 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_35 * 0.25 * f4_35, gto_indexes, 0, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 0.25 * f4_35, gto_indexes, 2, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_35 * f4_17, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_5 * f4_17, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_35 * 3.0 * f4_2, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_5 * 3.0 * f4_2, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_35 * 24.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_5 * 24.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_35 * 3.0 * f4_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_5 * 3.0 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -f4_35 * 4.0 * f4_2, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -f4_5 * 4.0 * f4_2, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -f4_35 * 8.0, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -f4_5 * 8.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XYYZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_XYYZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 3.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_17 * 3.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_2 * 0.5 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_17 * 0.5 * f4_5, gto_indexes, 7, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 0.25 * f4_35, gto_indexes, 5, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_17 * 0.25 * f4_35, gto_indexes, 7, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_2 * f4_35, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_17 * f4_35, gto_indexes, 7, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_2 * f4_5, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_17 * f4_5, gto_indexes, 7, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_17 * 3.0 * f4_2, gto_indexes, 7, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_2 * f4_17, gto_indexes, 5, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_17 * f4_17, gto_indexes, 7, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_2 * 6.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_17 * 6.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * f4_2 * 1.50 * f4_35, gto_indexes, 5, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * f4_17 * 1.50 * f4_35, gto_indexes, 7, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * f4_2 * 3.0 * f4_17, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * f4_17 * 3.0 * f4_17, gto_indexes, 7, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * f4_17 * 3.0 * f4_2, gto_indexes, 7, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_2 * 24.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_17 * 24.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * f4_2 * 3.0 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * f4_17 * 3.0 * f4_5, gto_indexes, 7, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_35, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_17 * f4_35, gto_indexes, 7, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_5, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_17 * f4_5, gto_indexes, 7, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_17 * 3.0 * f4_2, gto_indexes, 7, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_17, gto_indexes, 5, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_17 * 3.0 * f4_17, gto_indexes, 7, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -3.0 * f4_2 * 6.0 * f4_5, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -3.0 * f4_17 * 6.0 * f4_5, gto_indexes, 7, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -3.0 * f4_2 * 4.0 * f4_2, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -3.0 * f4_17 * 4.0 * f4_2, gto_indexes, 7, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 3.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_17 * 3.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 0.5 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_17 * 0.5 * f4_5, gto_indexes, 7, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 0.25 * f4_35, gto_indexes, 5, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_17 * 0.25 * f4_35, gto_indexes, 7, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * f4_17, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_17 * f4_17, gto_indexes, 7, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_17 * 3.0 * f4_2, gto_indexes, 7, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 24.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_17 * 24.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 3.0 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_17 * 3.0 * f4_5, gto_indexes, 7, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -3.0 * f4_2 * 4.0 * f4_2, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -3.0 * f4_17 * 4.0 * f4_2, gto_indexes, 7, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_2 * 8.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_17 * 8.0, gto_indexes, 7, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XYZZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_XYZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 6.0 * f4_5 * 3.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -6.0 * f4_5 * 0.5 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 6.0 * f4_5 * 0.25 * f4_35, gto_indexes, 2, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 6.0 * f4_5 * f4_35, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -6.0 * f4_5 * f4_5, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -6.0 * f4_5 * 3.0 * f4_2, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 6.0 * f4_5 * f4_17, gto_indexes, 2, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 6.0 * f4_5 * 6.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -6.0 * f4_5 * 1.50 * f4_35, gto_indexes, 2, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 6.0 * f4_5 * 3.0 * f4_17, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -6.0 * f4_5 * 3.0 * f4_2, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -6.0 * f4_5 * 24.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 6.0 * f4_5 * 3.0 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -6.0 * f4_5 * f4_35, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -6.0 * f4_5 * f4_5, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -6.0 * f4_5 * 3.0 * f4_2, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -6.0 * f4_5 * 3.0 * f4_17, gto_indexes, 2, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 6.0 * f4_5 * 6.0 * f4_5, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 6.0 * f4_5 * 4.0 * f4_2, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 6.0 * f4_5 * 3.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 6.0 * f4_5 * 0.5 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 6.0 * f4_5 * 0.25 * f4_35, gto_indexes, 2, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -6.0 * f4_5 * f4_17, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -6.0 * f4_5 * 3.0 * f4_2, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -6.0 * f4_5 * 24.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -6.0 * f4_5 * 3.0 * f4_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 6.0 * f4_5 * 4.0 * f4_2, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 6.0 * f4_5 * 8.0, gto_indexes, 2, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XZZZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_XZZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 4.0 * f4_2 * 3.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -4.0 * f4_2 * 0.5 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 4.0 * f4_2 * 0.25 * f4_35, gto_indexes, 5, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 4.0 * f4_2 * f4_35, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -4.0 * f4_2 * f4_5, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -4.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 4.0 * f4_2 * f4_17, gto_indexes, 5, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 4.0 * f4_2 * 6.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -4.0 * f4_2 * 1.50 * f4_35, gto_indexes, 5, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 4.0 * f4_2 * 3.0 * f4_17, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -4.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -4.0 * f4_2 * 24.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 4.0 * f4_2 * 3.0 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -4.0 * f4_2 * f4_35, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -4.0 * f4_2 * f4_5, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -4.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -4.0 * f4_2 * 3.0 * f4_17, gto_indexes, 5, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 4.0 * f4_2 * 6.0 * f4_5, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 4.0 * f4_2 * 4.0 * f4_2, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 4.0 * f4_2 * 3.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 4.0 * f4_2 * 0.5 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 4.0 * f4_2 * 0.25 * f4_35, gto_indexes, 5, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -4.0 * f4_2 * f4_17, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -4.0 * f4_2 * 3.0 * f4_2, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -4.0 * f4_2 * 24.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -4.0 * f4_2 * 3.0 * f4_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 4.0 * f4_2 * 4.0 * f4_2, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 4.0 * f4_2 * 8.0, gto_indexes, 5, 4, j, ket_first, ket_last);

            // compute primitive integrals block (YYYY)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_YYYY_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.5 * f4_5 * 3.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.25 * f4_35 * 3.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -0.5 * f4_5 * 0.5 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -0.25 * f4_35 * 0.5 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.5 * f4_5 * 0.25 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.25 * f4_35 * 0.25 * f4_35, gto_indexes, 8, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 0.5 * f4_5 * f4_35, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 0.25 * f4_35 * f4_35, gto_indexes, 8, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -0.5 * f4_5 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -0.25 * f4_35 * f4_5, gto_indexes, 8, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -0.5 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -0.25 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 0.5 * f4_5 * f4_17, gto_indexes, 6, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 0.25 * f4_35 * f4_17, gto_indexes, 8, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * 6.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 0.5 * f4_5 * 6.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 0.25 * f4_35 * 6.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * 1.50 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -0.5 * f4_5 * 1.50 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -0.25 * f4_35 * 1.50 * f4_35, gto_indexes, 8, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * 3.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 0.5 * f4_5 * 3.0 * f4_17, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 0.25 * f4_35 * 3.0 * f4_17, gto_indexes, 8, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -0.5 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -0.25 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -0.5 * f4_5 * 24.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -0.25 * f4_35 * 24.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 0.5 * f4_5 * 3.0 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 0.25 * f4_35 * 3.0 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.5 * f4_5 * f4_35, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.25 * f4_35 * f4_35, gto_indexes, 8, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.5 * f4_5 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.25 * f4_35 * f4_5, gto_indexes, 8, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -0.5 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -0.25 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * 3.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -0.5 * f4_5 * 3.0 * f4_17, gto_indexes, 6, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -0.25 * f4_35 * 3.0 * f4_17, gto_indexes, 8, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 3.0 * 6.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 0.5 * f4_5 * 6.0 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 0.25 * f4_35 * 6.0 * f4_5, gto_indexes, 8, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 3.0 * 4.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 0.5 * f4_5 * 4.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 0.25 * f4_35 * 4.0 * f4_2, gto_indexes, 8, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.5 * f4_5 * 3.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.25 * f4_35 * 3.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.5 * f4_5 * 0.5 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.25 * f4_35 * 0.5 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.5 * f4_5 * 0.25 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.25 * f4_35 * 0.25 * f4_35, gto_indexes, 8, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -0.5 * f4_5 * f4_17, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -0.25 * f4_35 * f4_17, gto_indexes, 8, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -0.5 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -0.25 * f4_35 * 3.0 * f4_2, gto_indexes, 8, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -0.5 * f4_5 * 24.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -0.25 * f4_35 * 24.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -0.5 * f4_5 * 3.0 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -0.25 * f4_35 * 3.0 * f4_5, gto_indexes, 8, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 3.0 * 4.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 0.5 * f4_5 * 4.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 0.25 * f4_35 * 4.0 * f4_2, gto_indexes, 8, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 3.0 * 8.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 0.5 * f4_5 * 8.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 0.25 * f4_35 * 8.0, gto_indexes, 8, 4, j, ket_first, ket_last);

            // compute primitive integrals block (YYYZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_YYYZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_17 * 3.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 3.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_17 * 0.5 * f4_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_2 * 0.5 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_17 * 0.25 * f4_35, gto_indexes, 1, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 0.25 * f4_35, gto_indexes, 3, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_17 * f4_35, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_2 * f4_35, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_17 * f4_5, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_2 * f4_5, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_17 * 3.0 * f4_2, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_17 * f4_17, gto_indexes, 1, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_2 * f4_17, gto_indexes, 3, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_17 * 6.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_2 * 6.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_17 * 1.50 * f4_35, gto_indexes, 1, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * f4_2 * 1.50 * f4_35, gto_indexes, 3, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_17 * 3.0 * f4_17, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * f4_2 * 3.0 * f4_17, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_17 * 3.0 * f4_2, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_17 * 24.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_2 * 24.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_17 * 3.0 * f4_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * f4_2 * 3.0 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_17 * f4_35, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_35, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_17 * f4_5, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_5, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_17 * 3.0 * f4_2, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_17 * 3.0 * f4_17, gto_indexes, 1, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_17, gto_indexes, 3, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -f4_17 * 6.0 * f4_5, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -3.0 * f4_2 * 6.0 * f4_5, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -f4_17 * 4.0 * f4_2, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -3.0 * f4_2 * 4.0 * f4_2, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_17 * 3.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 3.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_17 * 0.5 * f4_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 0.5 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_17 * 0.25 * f4_35, gto_indexes, 1, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 0.25 * f4_35, gto_indexes, 3, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_17 * f4_17, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * f4_17, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_17 * 3.0 * f4_2, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_17 * 24.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 24.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_17 * 3.0 * f4_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 3.0 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -f4_17 * 4.0 * f4_2, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -3.0 * f4_2 * 4.0 * f4_2, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -f4_17 * 8.0, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_2 * 8.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            // compute primitive integrals block (YYZZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_YYZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -24.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_5 * 3.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 24.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_5 * 0.5 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -24.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_5 * 0.25 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -24.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_5 * f4_35, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 24.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_5 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 24.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -24.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_5 * f4_17, gto_indexes, 6, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -24.0 * 6.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_5 * 6.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 24.0 * 1.50 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * f4_5 * 1.50 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -24.0 * 3.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * f4_5 * 3.0 * f4_17, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 24.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 24.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_5 * 24.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -24.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * f4_5 * 3.0 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 24.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_5 * f4_35, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 24.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_5 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 24.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 24.0 * 3.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, 3.0 * f4_5 * 3.0 * f4_17, gto_indexes, 6, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -24.0 * 6.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, -3.0 * f4_5 * 6.0 * f4_5, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -24.0 * 4.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, -3.0 * f4_5 * 4.0 * f4_2, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_5 * 3.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_5 * 0.5 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_5 * 0.25 * f4_35, gto_indexes, 6, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 24.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_5 * f4_17, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 24.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_5 * 3.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 24.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_5 * 24.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 24.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_5 * 3.0 * f4_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -24.0 * 4.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, -3.0 * f4_5 * 4.0 * f4_2, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -24.0 * 8.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_5 * 8.0, gto_indexes, 6, 4, j, ket_first, ket_last);

            // compute primitive integrals block (YZZZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_YZZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 4.0 * f4_2 * 3.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -4.0 * f4_2 * 0.5 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 4.0 * f4_2 * 0.25 * f4_35, gto_indexes, 3, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 4.0 * f4_2 * f4_35, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -4.0 * f4_2 * f4_5, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -4.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 4.0 * f4_2 * f4_17, gto_indexes, 3, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 4.0 * f4_2 * 6.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -4.0 * f4_2 * 1.50 * f4_35, gto_indexes, 3, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 4.0 * f4_2 * 3.0 * f4_17, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -4.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -4.0 * f4_2 * 24.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 4.0 * f4_2 * 3.0 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -4.0 * f4_2 * f4_35, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -4.0 * f4_2 * f4_5, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -4.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -4.0 * f4_2 * 3.0 * f4_17, gto_indexes, 3, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 4.0 * f4_2 * 6.0 * f4_5, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 4.0 * f4_2 * 4.0 * f4_2, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 4.0 * f4_2 * 3.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 4.0 * f4_2 * 0.5 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 4.0 * f4_2 * 0.25 * f4_35, gto_indexes, 3, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -4.0 * f4_2 * f4_17, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -4.0 * f4_2 * 3.0 * f4_2, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -4.0 * f4_2 * 24.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -4.0 * f4_2 * 3.0 * f4_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 4.0 * f4_2 * 4.0 * f4_2, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 4.0 * f4_2 * 8.0, gto_indexes, 3, 4, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGG_ZZZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 8.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, -8.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxx, 8.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, 8.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxy, -8.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, -8.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxxz, 8.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, 8.0 * 6.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyy, -8.0 * 1.50 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, 8.0 * 3.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxyz, -8.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, -8.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxzz, 8.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -8.0 * f4_35, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyy, -8.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -8.0 * 3.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyyz, -8.0 * 3.0 * f4_17, gto_indexes, 4, 7, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyzz, 8.0 * 6.0 * f4_5, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzzz, 8.0 * 4.0 * f4_2, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 8.0 * 3.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 8.0 * 0.5 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyy, 8.0 * 0.25 * f4_35, gto_indexes, 4, 8, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -8.0 * f4_17, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyyz, -8.0 * 3.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -8.0 * 24.0, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyzz, -8.0 * 3.0 * f4_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzzz, 8.0 * 4.0 * f4_2, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzzz, 8.0 * 8.0, gto_indexes, 4, 4, j, ket_first, ket_last);
        }
    }
}

auto
compOverlapGG(CSubMatrix*      matrix,
              const CGtoBlock& bra_gto_block,
              const CGtoBlock& ket_gto_block,
              const int64_t    bra_first,
              const int64_t    bra_last,
              const mat_t      mat_type) -> void
{
    // spherical transformation factors

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

            // compute primitive integrals block (XXXX)

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

                    ovlrec::compPrimitiveOverlapGG_XXXX_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 0.5 * f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -0.25 * f4_35 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -0.5 * f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 0.25 * f4_35 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -0.5 * f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 0.25 * f4_35 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 0.5 * f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -0.25 * f4_35 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, 0.5 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, -0.25 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -0.5 * f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 0.25 * f4_35 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -0.5 * f4_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 0.25 * f4_35 * 6.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, 0.5 * f4_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, -0.25 * f4_35 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -0.5 * f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 0.25 * f4_35 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 0.5 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -0.25 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 0.5 * f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -0.25 * f4_35 * 24.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, -0.5 * f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, 0.25 * f4_35 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 0.5 * f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.25 * f4_35 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 0.5 * f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.25 * f4_35 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 0.5 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -0.25 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 0.5 * f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -0.25 * f4_35 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, 3.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, -0.5 * f4_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, 0.25 * f4_35 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, 3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, -0.5 * f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, 0.25 * f4_35 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -0.5 * f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 0.25 * f4_35 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -0.5 * f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 0.25 * f4_35 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 0.5 * f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -0.25 * f4_35 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, 0.5 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, -0.25 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 0.5 * f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -0.25 * f4_35 * 24.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, 0.5 * f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, -0.25 * f4_35 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, 3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, -0.5 * f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, 0.25 * f4_35 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 3.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -0.5 * f4_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 0.25 * f4_35 * 8.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXXY)

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

                    ovlrec::compPrimitiveOverlapGG_XXXY_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_35 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_35 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_35 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_35 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_35 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_35 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_35 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_35 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_35 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_35 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_35 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_35 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -f4_35 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, f4_35 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, -f4_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, f4_35 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, -f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_35 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_35 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -f4_35 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -f4_35 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -f4_35 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, f4_35 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, -f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, f4_35 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -f4_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXXZ)

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

                    ovlrec::compPrimitiveOverlapGG_XXXZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 3.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_17 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -3.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_17 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 7, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_17 * f4_35, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_17 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_17 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_2 * 6.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_17 * 6.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, 3.0 * f4_2 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_17 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 7, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -3.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_17 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_17 * 24.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, -3.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_17 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_17 * f4_35, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -f4_17 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -f4_17 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, -3.0 * f4_2 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, f4_17 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, -3.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, f4_17 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_17 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, f4_17 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 7, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -f4_17 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -f4_17 * 24.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, 3.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -f4_17 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, -3.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, f4_17 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_2 * 8.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, f4_17 * 8.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXYY)

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

                    ovlrec::compPrimitiveOverlapGG_XXYY_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 6.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -6.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 1.50 * f4_35 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, 6.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -1.50 * f4_35 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 6.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -1.50 * f4_35 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 1.50 * f4_35 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -6.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, 1.50 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 6.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -1.50 * f4_35 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 6.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -1.50 * f4_35 * 6.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -6.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, 1.50 * f4_35 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, 6.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -1.50 * f4_35 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -6.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 1.50 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -6.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 1.50 * f4_35 * 24.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 6.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, -1.50 * f4_35 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -6.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 1.50 * f4_35 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 1.50 * f4_35 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -6.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 1.50 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -6.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 1.50 * f4_35 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, 6.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, -1.50 * f4_35 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, 6.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, -1.50 * f4_35 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 6.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -1.50 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 6.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -1.50 * f4_35 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 6.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -1.50 * f4_35 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -6.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 1.50 * f4_35 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -6.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, 1.50 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -6.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 1.50 * f4_35 * 24.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -6.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, 1.50 * f4_35 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, 6.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, -1.50 * f4_35 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 6.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -1.50 * f4_35 * 8.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXYZ)

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

                    ovlrec::compPrimitiveOverlapGG_XXYZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -3.0 * f4_17 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 3.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 3.0 * f4_17 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -3.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_17 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_17 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, -3.0 * f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_17 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * f4_17 * 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_2 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, -3.0 * f4_17 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, 3.0 * f4_2 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 3.0 * f4_17 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -3.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -3.0 * f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * f4_17 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, 3.0 * f4_17 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, -3.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_17 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_17 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -3.0 * f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -3.0 * f4_17 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, 3.0 * f4_17 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, -3.0 * f4_2 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, 3.0 * f4_17 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, -3.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 3.0 * f4_17 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 3.0 * f4_17 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_17 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, -3.0 * f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * f4_17 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, -3.0 * f4_17 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, 3.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, 3.0 * f4_17 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, -3.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 3.0 * f4_17 * 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_2 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZZ)

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

                    ovlrec::compPrimitiveOverlapGG_XXZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, 24.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -3.0 * f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -24.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 3.0 * f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -24.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 24.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 24.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, -3.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -24.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -24.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * f4_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 24.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, -3.0 * f4_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -24.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 3.0 * f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, 24.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -3.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 24.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -24.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, 3.0 * f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 24.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 24.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, 24.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -3.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, 24.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -3.0 * f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, -24.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, 3.0 * f4_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, -24.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, 3.0 * f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 3.0 * f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 3.0 * f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 24.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 24.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, -3.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 24.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 24.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, -3.0 * f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, -24.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, 3.0 * f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -24.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 3.0 * f4_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYYY)

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

                    ovlrec::compPrimitiveOverlapGG_XYYY_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_35 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_35 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_35 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_35 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_35 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_35 * 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_35 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_35 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_35 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_35 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_35 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_35 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_35 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, -f4_35 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, -f4_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, -f4_35 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, -f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_35 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_35 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_35 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_35 * 24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_35 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, -f4_35 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, -f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -f4_35 * 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -f4_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYYZ)

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

                    ovlrec::compPrimitiveOverlapGG_XYYZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 3.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 3.0 * f4_17 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -3.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -3.0 * f4_17 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 7, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_17 * f4_35, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_17 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, 3.0 * f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_17 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_2 * 6.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_17 * 6.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, 3.0 * f4_2 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, 3.0 * f4_17 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 7, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -3.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -3.0 * f4_17 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 3.0 * f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_17 * 24.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, -3.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, -3.0 * f4_17 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_17 * f4_35, bra_gto_indexes, ket_gto_indexes, 7, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_17 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_17 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, -3.0 * f4_2 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, -3.0 * f4_17 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, -3.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, -3.0 * f4_17 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_17 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_17 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 7, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_17 * f4_17, bra_gto_indexes, ket_gto_indexes, 7, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, 3.0 * f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_17 * 24.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, 3.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, 3.0 * f4_17 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 7, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, -3.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, -3.0 * f4_17 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 7, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_2 * 8.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_17 * 8.0, bra_gto_indexes, ket_gto_indexes, 7, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZZ)

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

                    ovlrec::compPrimitiveOverlapGG_XYZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -6.0 * f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 6.0 * f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 6.0 * f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -6.0 * f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, -6.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 6.0 * f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 6.0 * f4_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, -6.0 * f4_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 6.0 * f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -6.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -6.0 * f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, 6.0 * f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -6.0 * f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -6.0 * f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -6.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -6.0 * f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, 6.0 * f4_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, 6.0 * f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 6.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 6.0 * f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 6.0 * f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -6.0 * f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, -6.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -6.0 * f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, -6.0 * f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, 6.0 * f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 6.0 * f4_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZZ)

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

                    ovlrec::compPrimitiveOverlapGG_XZZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -4.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 4.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 4.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -4.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, -4.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 4.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 4.0 * f4_2 * 6.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, -4.0 * f4_2 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 4.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -4.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -4.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, 4.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -4.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -4.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -4.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -4.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, 4.0 * f4_2 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, 4.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 4.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 4.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 5, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -4.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, -4.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -4.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, -4.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, 4.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 4.0 * f4_2 * 8.0, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYYY)

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

                    ovlrec::compPrimitiveOverlapGG_YYYY_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, 0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -0.5 * f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -0.25 * f4_35 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, 3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 0.5 * f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 0.25 * f4_35 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 0.5 * f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 0.25 * f4_35 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -0.5 * f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -0.25 * f4_35 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, -0.5 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, -0.25 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 0.5 * f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 0.25 * f4_35 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 3.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 0.5 * f4_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 0.25 * f4_35 * 6.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, -0.5 * f4_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, -0.25 * f4_35 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, 3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 0.5 * f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 0.25 * f4_35 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -0.5 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -0.25 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -0.5 * f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -0.25 * f4_35 * 24.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, 0.5 * f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, 0.25 * f4_35 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.5 * f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.25 * f4_35 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.5 * f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -0.25 * f4_35 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -0.5 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -0.25 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -3.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -0.5 * f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -0.25 * f4_35 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, 3.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, 0.5 * f4_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, 0.25 * f4_35 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, 3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, 0.5 * f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, 0.25 * f4_35 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.5 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 0.25 * f4_35 * 3.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 0.5 * f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 0.25 * f4_35 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 3.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 0.5 * f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 0.25 * f4_35 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 8, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -0.5 * f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -0.25 * f4_35 * f4_17, bra_gto_indexes, ket_gto_indexes, 8, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -3.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, -0.5 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, -0.25 * f4_35 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -0.5 * f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -0.25 * f4_35 * 24.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -3.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, -0.5 * f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, -0.25 * f4_35 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 8, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, 3.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, 0.5 * f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, 0.25 * f4_35 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 8, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 3.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 0.5 * f4_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 0.25 * f4_35 * 8.0, bra_gto_indexes, ket_gto_indexes, 8, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYYZ)

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

                    ovlrec::compPrimitiveOverlapGG_YYYZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, f4_17 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 3.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -f4_17 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -3.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -f4_17 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, f4_17 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -f4_17 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -f4_17 * 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_2 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, f4_17 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, 3.0 * f4_2 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -f4_17 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -3.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, f4_17 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -f4_17 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, -3.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_17 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, f4_17 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, f4_17 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, -f4_17 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, -3.0 * f4_2 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, -f4_17 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, -3.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_17 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_17 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -f4_17 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_17 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, f4_17 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, 3.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_17 * 24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, f4_17 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, 3.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, -f4_17 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, -3.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -f4_17 * 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_2 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZZ)

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

                    ovlrec::compPrimitiveOverlapGG_YYZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, -24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, 24.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 3.0 * f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -24.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -3.0 * f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -24.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -3.0 * f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 24.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 3.0 * f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 24.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, 3.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -24.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -3.0 * f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -24.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -3.0 * f4_5 * 6.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 24.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, 3.0 * f4_5 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -24.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -3.0 * f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, 24.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 3.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 24.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 3.0 * f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -24.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, -3.0 * f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 24.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_5 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 24.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, 3.0 * f4_5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, 24.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, 24.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, 3.0 * f4_5 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, -24.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, -3.0 * f4_5 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, -24.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, -3.0 * f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -3.0 * f4_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_5 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, -24.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, -3.0 * f4_5 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 6, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 24.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 3.0 * f4_5 * f4_17, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, 24.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, 3.0 * f4_5 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 24.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 3.0 * f4_5 * 24.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, 24.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, 3.0 * f4_5 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, -24.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, -3.0 * f4_5 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -24.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, -3.0 * f4_5 * 8.0, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZZ)

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

                    ovlrec::compPrimitiveOverlapGG_YZZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, -4.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxx, 4.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 4.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -4.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxxz, -4.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 4.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 4.0 * f4_2 * 6.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyy, -4.0 * f4_2 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, 4.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxyz, -4.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -4.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxzz, 4.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -4.0 * f4_2 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -4.0 * f4_2 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -4.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyyz, -4.0 * f4_2 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyzz, 4.0 * f4_2 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzzz, 4.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 4.0 * f4_2 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 4.0 * f4_2 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyy, 4.0 * f4_2 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 3, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -4.0 * f4_2 * f4_17, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyyz, -4.0 * f4_2 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -4.0 * f4_2 * 24.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyzz, -4.0 * f4_2 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzzz, 4.0 * f4_2 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 4.0 * f4_2 * 8.0, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZZ)

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

                    ovlrec::compPrimitiveOverlapGG_ZZZZ_T(buffer_xxxx,
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

            t2cfunc::distribute(matrix, buffer_xxxx, 8.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, -8.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxx, 8.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, 8.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxy, -8.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, -8.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxxz, 8.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, 8.0 * 6.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyy, -8.0 * 1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, 8.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxyz, -8.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, -8.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxzz, 8.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -8.0 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyy, -8.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -8.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyyz, -8.0 * 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 7, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyzz, 8.0 * 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzzz, 8.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 8.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 8.0 * 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyy, 8.0 * 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 4, 8, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -8.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyyz, -8.0 * 3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -8.0 * 24.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyzz, -8.0 * 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzzz, 8.0 * 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzzz, 8.0 * 8.0, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);
        }
    }
}

}  // namespace ovlrec
