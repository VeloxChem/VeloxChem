#include "OverlapGeom301RecDD.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveOverlapGeom301DD_XX_XX.hpp"
#include "PrimitiveOverlapGeom301DD_XX_XY.hpp"
#include "PrimitiveOverlapGeom301DD_XX_XZ.hpp"
#include "PrimitiveOverlapGeom301DD_XX_YY.hpp"
#include "PrimitiveOverlapGeom301DD_XX_YZ.hpp"
#include "PrimitiveOverlapGeom301DD_XX_ZZ.hpp"
#include "PrimitiveOverlapGeom301DD_XY_XX.hpp"
#include "PrimitiveOverlapGeom301DD_XY_XY.hpp"
#include "PrimitiveOverlapGeom301DD_XY_XZ.hpp"
#include "PrimitiveOverlapGeom301DD_XY_YY.hpp"
#include "PrimitiveOverlapGeom301DD_XY_YZ.hpp"
#include "PrimitiveOverlapGeom301DD_XY_ZZ.hpp"
#include "PrimitiveOverlapGeom301DD_XZ_XX.hpp"
#include "PrimitiveOverlapGeom301DD_XZ_XY.hpp"
#include "PrimitiveOverlapGeom301DD_XZ_XZ.hpp"
#include "PrimitiveOverlapGeom301DD_XZ_YY.hpp"
#include "PrimitiveOverlapGeom301DD_XZ_YZ.hpp"
#include "PrimitiveOverlapGeom301DD_XZ_ZZ.hpp"
#include "PrimitiveOverlapGeom301DD_YY_XX.hpp"
#include "PrimitiveOverlapGeom301DD_YY_XY.hpp"
#include "PrimitiveOverlapGeom301DD_YY_XZ.hpp"
#include "PrimitiveOverlapGeom301DD_YY_YY.hpp"
#include "PrimitiveOverlapGeom301DD_YY_YZ.hpp"
#include "PrimitiveOverlapGeom301DD_YY_ZZ.hpp"
#include "PrimitiveOverlapGeom301DD_YZ_XX.hpp"
#include "PrimitiveOverlapGeom301DD_YZ_XY.hpp"
#include "PrimitiveOverlapGeom301DD_YZ_XZ.hpp"
#include "PrimitiveOverlapGeom301DD_YZ_YY.hpp"
#include "PrimitiveOverlapGeom301DD_YZ_YZ.hpp"
#include "PrimitiveOverlapGeom301DD_YZ_ZZ.hpp"
#include "PrimitiveOverlapGeom301DD_ZZ_XX.hpp"
#include "PrimitiveOverlapGeom301DD_ZZ_XY.hpp"
#include "PrimitiveOverlapGeom301DD_ZZ_XZ.hpp"
#include "PrimitiveOverlapGeom301DD_ZZ_YY.hpp"
#include "PrimitiveOverlapGeom301DD_ZZ_YZ.hpp"
#include "PrimitiveOverlapGeom301DD_ZZ_ZZ.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGeom301DD(CSubMatrix*      matrix_xxx_x,
                     CSubMatrix*      matrix_xxx_y,
                     CSubMatrix*      matrix_xxx_z,
                     CSubMatrix*      matrix_xxy_x,
                     CSubMatrix*      matrix_xxy_y,
                     CSubMatrix*      matrix_xxy_z,
                     CSubMatrix*      matrix_xxz_x,
                     CSubMatrix*      matrix_xxz_y,
                     CSubMatrix*      matrix_xxz_z,
                     CSubMatrix*      matrix_xyy_x,
                     CSubMatrix*      matrix_xyy_y,
                     CSubMatrix*      matrix_xyy_z,
                     CSubMatrix*      matrix_xyz_x,
                     CSubMatrix*      matrix_xyz_y,
                     CSubMatrix*      matrix_xyz_z,
                     CSubMatrix*      matrix_xzz_x,
                     CSubMatrix*      matrix_xzz_y,
                     CSubMatrix*      matrix_xzz_z,
                     CSubMatrix*      matrix_yyy_x,
                     CSubMatrix*      matrix_yyy_y,
                     CSubMatrix*      matrix_yyy_z,
                     CSubMatrix*      matrix_yyz_x,
                     CSubMatrix*      matrix_yyz_y,
                     CSubMatrix*      matrix_yyz_z,
                     CSubMatrix*      matrix_yzz_x,
                     CSubMatrix*      matrix_yzz_y,
                     CSubMatrix*      matrix_yzz_z,
                     CSubMatrix*      matrix_zzz_x,
                     CSubMatrix*      matrix_zzz_y,
                     CSubMatrix*      matrix_zzz_z,
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

    alignas(64) TDoubleArray buffer_xxx_x;

    alignas(64) TDoubleArray buffer_xxx_y;

    alignas(64) TDoubleArray buffer_xxx_z;

    alignas(64) TDoubleArray buffer_xxy_x;

    alignas(64) TDoubleArray buffer_xxy_y;

    alignas(64) TDoubleArray buffer_xxy_z;

    alignas(64) TDoubleArray buffer_xxz_x;

    alignas(64) TDoubleArray buffer_xxz_y;

    alignas(64) TDoubleArray buffer_xxz_z;

    alignas(64) TDoubleArray buffer_xyy_x;

    alignas(64) TDoubleArray buffer_xyy_y;

    alignas(64) TDoubleArray buffer_xyy_z;

    alignas(64) TDoubleArray buffer_xyz_x;

    alignas(64) TDoubleArray buffer_xyz_y;

    alignas(64) TDoubleArray buffer_xyz_z;

    alignas(64) TDoubleArray buffer_xzz_x;

    alignas(64) TDoubleArray buffer_xzz_y;

    alignas(64) TDoubleArray buffer_xzz_z;

    alignas(64) TDoubleArray buffer_yyy_x;

    alignas(64) TDoubleArray buffer_yyy_y;

    alignas(64) TDoubleArray buffer_yyy_z;

    alignas(64) TDoubleArray buffer_yyz_x;

    alignas(64) TDoubleArray buffer_yyz_y;

    alignas(64) TDoubleArray buffer_yyz_z;

    alignas(64) TDoubleArray buffer_yzz_x;

    alignas(64) TDoubleArray buffer_yzz_y;

    alignas(64) TDoubleArray buffer_yzz_z;

    alignas(64) TDoubleArray buffer_zzz_x;

    alignas(64) TDoubleArray buffer_zzz_y;

    alignas(64) TDoubleArray buffer_zzz_z;

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

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XX_XX(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XX_XY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XX_XY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XX_XZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XX_XZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XX_YY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XX_YY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XX_YZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XX_YZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XX_ZZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XX_ZZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_XX)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XY_XX(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_XY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XY_XY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_XZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XY_XZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_YY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XY_YY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_YZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XY_YZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY_ZZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XY_ZZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_XX)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XZ_XX(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_XY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XZ_XY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_XZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XZ_XZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_YY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XZ_YY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_YZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XZ_YZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ_ZZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_XZ_ZZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_XX)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YY_XX(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_XY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YY_XY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_XZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YY_XZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_YY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YY_YY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_YZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YY_YZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY_ZZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YY_ZZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_XX)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YZ_XX(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_XY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YZ_XY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_XZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YZ_XZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_YY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YZ_YY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_YZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YZ_YZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ_ZZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_YZ_ZZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_XX)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_ZZ_XX(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_XY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_ZZ_XY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_XZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_ZZ_XZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_YY)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_ZZ_YY(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_x, buffer_xxx_x, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_y, buffer_xxx_y, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxx_z, buffer_xxx_z, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_x, buffer_xxy_x, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_y, buffer_xxy_y, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxy_z, buffer_xxy_z, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_x, buffer_xxz_x, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_y, buffer_xxz_y, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xxz_z, buffer_xxz_z, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_x, buffer_xyy_x, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_y, buffer_xyy_y, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyy_z, buffer_xyy_z, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_x, buffer_xyz_x, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_y, buffer_xyz_y, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xyz_z, buffer_xyz_z, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_x, buffer_xzz_x, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_y, buffer_xzz_y, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_xzz_z, buffer_xzz_z, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_x, buffer_yyy_x, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_y, buffer_yyy_y, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyy_z, buffer_yyy_z, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_x, buffer_yyz_x, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_y, buffer_yyz_y, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yyz_z, buffer_yyz_z, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_x, buffer_yzz_x, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_y, buffer_yzz_y, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_yzz_z, buffer_yzz_z, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_x, buffer_zzz_x, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_y, buffer_zzz_y, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_zzz_z, buffer_zzz_z, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_YZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_ZZ_YZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ_ZZ)

            simd::zero(buffer_xxx_x);

            simd::zero(buffer_xxx_y);

            simd::zero(buffer_xxx_z);

            simd::zero(buffer_xxy_x);

            simd::zero(buffer_xxy_y);

            simd::zero(buffer_xxy_z);

            simd::zero(buffer_xxz_x);

            simd::zero(buffer_xxz_y);

            simd::zero(buffer_xxz_z);

            simd::zero(buffer_xyy_x);

            simd::zero(buffer_xyy_y);

            simd::zero(buffer_xyy_z);

            simd::zero(buffer_xyz_x);

            simd::zero(buffer_xyz_y);

            simd::zero(buffer_xyz_z);

            simd::zero(buffer_xzz_x);

            simd::zero(buffer_xzz_y);

            simd::zero(buffer_xzz_z);

            simd::zero(buffer_yyy_x);

            simd::zero(buffer_yyy_y);

            simd::zero(buffer_yyy_z);

            simd::zero(buffer_yyz_x);

            simd::zero(buffer_yyz_y);

            simd::zero(buffer_yyz_z);

            simd::zero(buffer_yzz_x);

            simd::zero(buffer_yzz_y);

            simd::zero(buffer_yzz_z);

            simd::zero(buffer_zzz_x);

            simd::zero(buffer_zzz_y);

            simd::zero(buffer_zzz_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    ovlrec::compPrimitiveOverlapGeom301DD_ZZ_ZZ(buffer_xxx_x,
                                                                buffer_xxx_y,
                                                                buffer_xxx_z,
                                                                buffer_xxy_x,
                                                                buffer_xxy_y,
                                                                buffer_xxy_z,
                                                                buffer_xxz_x,
                                                                buffer_xxz_y,
                                                                buffer_xxz_z,
                                                                buffer_xyy_x,
                                                                buffer_xyy_y,
                                                                buffer_xyy_z,
                                                                buffer_xyz_x,
                                                                buffer_xyz_y,
                                                                buffer_xyz_z,
                                                                buffer_xzz_x,
                                                                buffer_xzz_y,
                                                                buffer_xzz_z,
                                                                buffer_yyy_x,
                                                                buffer_yyy_y,
                                                                buffer_yyy_z,
                                                                buffer_yyz_x,
                                                                buffer_yyz_y,
                                                                buffer_yyz_z,
                                                                buffer_yzz_x,
                                                                buffer_yzz_y,
                                                                buffer_yzz_z,
                                                                buffer_zzz_x,
                                                                buffer_zzz_y,
                                                                buffer_zzz_z,
                                                                bra_exp,
                                                                bra_norm,
                                                                bra_coord,
                                                                ket_exps,
                                                                ket_norms,
                                                                ket_coords_x,
                                                                ket_coords_y,
                                                                ket_coords_z,
                                                                ket_dim);
                }
            }

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);
        }
    }
}

}  // namespace ovlrec
