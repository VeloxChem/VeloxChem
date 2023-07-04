#include "OctupoleRecPG.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveOctupolePG_X_XXXX.hpp"
#include "PrimitiveOctupolePG_X_XXXY.hpp"
#include "PrimitiveOctupolePG_X_XXXZ.hpp"
#include "PrimitiveOctupolePG_X_XXYY.hpp"
#include "PrimitiveOctupolePG_X_XXYZ.hpp"
#include "PrimitiveOctupolePG_X_XXZZ.hpp"
#include "PrimitiveOctupolePG_X_XYYY.hpp"
#include "PrimitiveOctupolePG_X_XYYZ.hpp"
#include "PrimitiveOctupolePG_X_XYZZ.hpp"
#include "PrimitiveOctupolePG_X_XZZZ.hpp"
#include "PrimitiveOctupolePG_X_YYYY.hpp"
#include "PrimitiveOctupolePG_X_YYYZ.hpp"
#include "PrimitiveOctupolePG_X_YYZZ.hpp"
#include "PrimitiveOctupolePG_X_YZZZ.hpp"
#include "PrimitiveOctupolePG_X_ZZZZ.hpp"
#include "PrimitiveOctupolePG_Y_XXXX.hpp"
#include "PrimitiveOctupolePG_Y_XXXY.hpp"
#include "PrimitiveOctupolePG_Y_XXXZ.hpp"
#include "PrimitiveOctupolePG_Y_XXYY.hpp"
#include "PrimitiveOctupolePG_Y_XXYZ.hpp"
#include "PrimitiveOctupolePG_Y_XXZZ.hpp"
#include "PrimitiveOctupolePG_Y_XYYY.hpp"
#include "PrimitiveOctupolePG_Y_XYYZ.hpp"
#include "PrimitiveOctupolePG_Y_XYZZ.hpp"
#include "PrimitiveOctupolePG_Y_XZZZ.hpp"
#include "PrimitiveOctupolePG_Y_YYYY.hpp"
#include "PrimitiveOctupolePG_Y_YYYZ.hpp"
#include "PrimitiveOctupolePG_Y_YYZZ.hpp"
#include "PrimitiveOctupolePG_Y_YZZZ.hpp"
#include "PrimitiveOctupolePG_Y_ZZZZ.hpp"
#include "PrimitiveOctupolePG_Z_XXXX.hpp"
#include "PrimitiveOctupolePG_Z_XXXY.hpp"
#include "PrimitiveOctupolePG_Z_XXXZ.hpp"
#include "PrimitiveOctupolePG_Z_XXYY.hpp"
#include "PrimitiveOctupolePG_Z_XXYZ.hpp"
#include "PrimitiveOctupolePG_Z_XXZZ.hpp"
#include "PrimitiveOctupolePG_Z_XYYY.hpp"
#include "PrimitiveOctupolePG_Z_XYYZ.hpp"
#include "PrimitiveOctupolePG_Z_XYZZ.hpp"
#include "PrimitiveOctupolePG_Z_XZZZ.hpp"
#include "PrimitiveOctupolePG_Z_YYYY.hpp"
#include "PrimitiveOctupolePG_Z_YYYZ.hpp"
#include "PrimitiveOctupolePG_Z_YYZZ.hpp"
#include "PrimitiveOctupolePG_Z_YZZZ.hpp"
#include "PrimitiveOctupolePG_Z_ZZZZ.hpp"
#include "T2CDistributor.hpp"

namespace mpol {  // mpol namespace

auto
compOctupolePG(CSubMatrix*      matrix_xxx,
               CSubMatrix*      matrix_xxy,
               CSubMatrix*      matrix_xxz,
               CSubMatrix*      matrix_xyy,
               CSubMatrix*      matrix_xyz,
               CSubMatrix*      matrix_xzz,
               CSubMatrix*      matrix_yyy,
               CSubMatrix*      matrix_yyz,
               CSubMatrix*      matrix_yzz,
               CSubMatrix*      matrix_zzz,
               const TPoint3D&  point,
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

            // compute primitive integrals block (X_XXXX)

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

                    mpol::compPrimitiveOctupolePG_X_XXXX(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXXY)

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

                    mpol::compPrimitiveOctupolePG_X_XXXY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXXZ)

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

                    mpol::compPrimitiveOctupolePG_X_XXXZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXYY)

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

                    mpol::compPrimitiveOctupolePG_X_XXYY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 6.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXYZ)

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

                    mpol::compPrimitiveOctupolePG_X_XXYZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXZZ)

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

                    mpol::compPrimitiveOctupolePG_X_XXZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XYYY)

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

                    mpol::compPrimitiveOctupolePG_X_XYYY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_35, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XYYZ)

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

                    mpol::compPrimitiveOctupolePG_X_XYYZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 2, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XYZZ)

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

                    mpol::compPrimitiveOctupolePG_X_XYZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XZZZ)

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

                    mpol::compPrimitiveOctupolePG_X_XZZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YYYY)

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

                    mpol::compPrimitiveOctupolePG_X_YYYY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 2, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YYYZ)

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

                    mpol::compPrimitiveOctupolePG_X_YYYZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_17, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YYZZ)

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

                    mpol::compPrimitiveOctupolePG_X_YYZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -24.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YZZZ)

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

                    mpol::compPrimitiveOctupolePG_X_YZZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_ZZZZ)

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

                    mpol::compPrimitiveOctupolePG_X_ZZZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 8.0, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXXX)

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

                    mpol::compPrimitiveOctupolePG_Y_XXXX(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXXY)

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

                    mpol::compPrimitiveOctupolePG_Y_XXXY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXXZ)

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

                    mpol::compPrimitiveOctupolePG_Y_XXXZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXYY)

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

                    mpol::compPrimitiveOctupolePG_Y_XXYY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 6.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXYZ)

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

                    mpol::compPrimitiveOctupolePG_Y_XXYZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXZZ)

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

                    mpol::compPrimitiveOctupolePG_Y_XXZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XYYY)

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

                    mpol::compPrimitiveOctupolePG_Y_XYYY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_35, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XYYZ)

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

                    mpol::compPrimitiveOctupolePG_Y_XYYZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 0, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XYZZ)

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

                    mpol::compPrimitiveOctupolePG_Y_XYZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XZZZ)

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

                    mpol::compPrimitiveOctupolePG_Y_XZZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YYYY)

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

                    mpol::compPrimitiveOctupolePG_Y_YYYY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 0, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YYYZ)

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

                    mpol::compPrimitiveOctupolePG_Y_YYYZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_17, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YYZZ)

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

                    mpol::compPrimitiveOctupolePG_Y_YYZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -24.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YZZZ)

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

                    mpol::compPrimitiveOctupolePG_Y_YZZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_ZZZZ)

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

                    mpol::compPrimitiveOctupolePG_Y_ZZZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 8.0, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXXX)

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

                    mpol::compPrimitiveOctupolePG_Z_XXXX(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXXY)

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

                    mpol::compPrimitiveOctupolePG_Z_XXXY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXXZ)

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

                    mpol::compPrimitiveOctupolePG_Z_XXXZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXYY)

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

                    mpol::compPrimitiveOctupolePG_Z_XXYY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 6.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -1.50 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXYZ)

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

                    mpol::compPrimitiveOctupolePG_Z_XXYZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXZZ)

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

                    mpol::compPrimitiveOctupolePG_Z_XXZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XYYY)

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

                    mpol::compPrimitiveOctupolePG_Z_XYYY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_35, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XYYZ)

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

                    mpol::compPrimitiveOctupolePG_Z_XYYZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_17, bra_gto_indexes, ket_gto_indexes, 1, 7, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XYZZ)

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

                    mpol::compPrimitiveOctupolePG_Z_XYZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 6.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XZZZ)

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

                    mpol::compPrimitiveOctupolePG_Z_XZZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YYYY)

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

                    mpol::compPrimitiveOctupolePG_Z_YYYY(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.5 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.25 * f4_35, bra_gto_indexes, ket_gto_indexes, 1, 8, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YYYZ)

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

                    mpol::compPrimitiveOctupolePG_Z_YYYZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f4_17, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YYZZ)

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

                    mpol::compPrimitiveOctupolePG_Z_YYZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -24.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f4_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YZZZ)

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

                    mpol::compPrimitiveOctupolePG_Z_YZZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f4_2, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_ZZZZ)

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

                    mpol::compPrimitiveOctupolePG_Z_ZZZZ(buffer_xxx,
                                                         buffer_xxy,
                                                         buffer_xxz,
                                                         buffer_xyy,
                                                         buffer_xyz,
                                                         buffer_xzz,
                                                         buffer_yyy,
                                                         buffer_yyz,
                                                         buffer_yzz,
                                                         buffer_zzz,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 8.0, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace mpol
