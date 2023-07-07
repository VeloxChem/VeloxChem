#include "OverlapGeom301RecSF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveOverlapGeom301SF_0_XXX.hpp"
#include "PrimitiveOverlapGeom301SF_0_XXY.hpp"
#include "PrimitiveOverlapGeom301SF_0_XXZ.hpp"
#include "PrimitiveOverlapGeom301SF_0_XYY.hpp"
#include "PrimitiveOverlapGeom301SF_0_XYZ.hpp"
#include "PrimitiveOverlapGeom301SF_0_XZZ.hpp"
#include "PrimitiveOverlapGeom301SF_0_YYY.hpp"
#include "PrimitiveOverlapGeom301SF_0_YYZ.hpp"
#include "PrimitiveOverlapGeom301SF_0_YZZ.hpp"
#include "PrimitiveOverlapGeom301SF_0_ZZZ.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compOverlapGeom301SF(CSubMatrix*      matrix_xxx_x,
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
                     const bool       ang_order,
                     const int64_t    bra_first,
                     const int64_t    bra_last) -> void
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

            // compute primitive integrals block (0_XXX)

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

                    ovlrec::compPrimitiveOverlapGeom301SF_0_XXX(buffer_xxx_x,
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

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_XXY)

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

                    ovlrec::compPrimitiveOverlapGeom301SF_0_XXY(buffer_xxx_x,
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

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_XXZ)

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

                    ovlrec::compPrimitiveOverlapGeom301SF_0_XXZ(buffer_xxx_x,
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

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_XYY)

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

                    ovlrec::compPrimitiveOverlapGeom301SF_0_XYY(buffer_xxx_x,
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

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_XYZ)

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

                    ovlrec::compPrimitiveOverlapGeom301SF_0_XYZ(buffer_xxx_x,
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

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_XZZ)

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

                    ovlrec::compPrimitiveOverlapGeom301SF_0_XZZ(buffer_xxx_x,
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

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_YYY)

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

                    ovlrec::compPrimitiveOverlapGeom301SF_0_YYY(buffer_xxx_x,
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

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_YYZ)

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

                    ovlrec::compPrimitiveOverlapGeom301SF_0_YYZ(buffer_xxx_x,
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

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_YZZ)

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

                    ovlrec::compPrimitiveOverlapGeom301SF_0_YZZ(buffer_xxx_x,
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

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (0_ZZZ)

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

                    ovlrec::compPrimitiveOverlapGeom301SF_0_ZZZ(buffer_xxx_x,
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

            t2cfunc::distribute(matrix_xxx_x, buffer_xxx_x, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_y, buffer_xxx_y, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx_z, buffer_xxx_z, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_x, buffer_xxy_x, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_y, buffer_xxy_y, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy_z, buffer_xxy_z, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_x, buffer_xxz_x, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_y, buffer_xxz_y, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz_z, buffer_xxz_z, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_x, buffer_xyy_x, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_y, buffer_xyy_y, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy_z, buffer_xyy_z, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_x, buffer_xyz_x, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_y, buffer_xyz_y, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz_z, buffer_xyz_z, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_x, buffer_xzz_x, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_y, buffer_xzz_y, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz_z, buffer_xzz_z, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_x, buffer_yyy_x, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_y, buffer_yyy_y, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy_z, buffer_yyy_z, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_x, buffer_yyz_x, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_y, buffer_yyz_y, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz_z, buffer_yyz_z, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_x, buffer_yzz_x, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_y, buffer_yzz_y, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz_z, buffer_yzz_z, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_x, buffer_zzz_x, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_y, buffer_zzz_y, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz_z, buffer_zzz_z, 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, ang_order);
        }
    }
}

}  // namespace ovlrec
