#include "OctupoleRecPF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "T2CDistributor.hpp"

#include "PrimitiveOctupolePF_X_XXX.hpp"
#include "PrimitiveOctupolePF_X_XXY.hpp"
#include "PrimitiveOctupolePF_X_XXZ.hpp"
#include "PrimitiveOctupolePF_X_XYY.hpp"
#include "PrimitiveOctupolePF_X_XYZ.hpp"
#include "PrimitiveOctupolePF_X_XZZ.hpp"
#include "PrimitiveOctupolePF_X_YYY.hpp"
#include "PrimitiveOctupolePF_X_YYZ.hpp"
#include "PrimitiveOctupolePF_X_YZZ.hpp"
#include "PrimitiveOctupolePF_X_ZZZ.hpp"
#include "PrimitiveOctupolePF_Y_XXX.hpp"
#include "PrimitiveOctupolePF_Y_XXY.hpp"
#include "PrimitiveOctupolePF_Y_XXZ.hpp"
#include "PrimitiveOctupolePF_Y_XYY.hpp"
#include "PrimitiveOctupolePF_Y_XYZ.hpp"
#include "PrimitiveOctupolePF_Y_XZZ.hpp"
#include "PrimitiveOctupolePF_Y_YYY.hpp"
#include "PrimitiveOctupolePF_Y_YYZ.hpp"
#include "PrimitiveOctupolePF_Y_YZZ.hpp"
#include "PrimitiveOctupolePF_Y_ZZZ.hpp"
#include "PrimitiveOctupolePF_Z_XXX.hpp"
#include "PrimitiveOctupolePF_Z_XXY.hpp"
#include "PrimitiveOctupolePF_Z_XXZ.hpp"
#include "PrimitiveOctupolePF_Z_XYY.hpp"
#include "PrimitiveOctupolePF_Z_XYZ.hpp"
#include "PrimitiveOctupolePF_Z_XZZ.hpp"
#include "PrimitiveOctupolePF_Z_YYY.hpp"
#include "PrimitiveOctupolePF_Z_YYZ.hpp"
#include "PrimitiveOctupolePF_Z_YZZ.hpp"
#include "PrimitiveOctupolePF_Z_ZZZ.hpp"

namespace octurec { // octurec namespace

auto
compOctupolePF(      CSubMatrix* matrix_xxx,
                     CSubMatrix* matrix_xxy,
                     CSubMatrix* matrix_xxz,
                     CSubMatrix* matrix_xyy,
                     CSubMatrix* matrix_xyz,
                     CSubMatrix* matrix_xzz,
                     CSubMatrix* matrix_yyy,
                     CSubMatrix* matrix_yyz,
                     CSubMatrix* matrix_yzz,
                     CSubMatrix* matrix_zzz,
               const TPoint3D& point,
               const CGtoBlock&  bra_gto_block,
               const CGtoBlock&  ket_gto_block,
               const bool        ang_order,
               const int64_t     bra_first,
               const int64_t     bra_last) -> void
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

        simd::loadCoordinates(ket_coords_x,
                              ket_coords_y,
                              ket_coords_z,
                              ket_gto_coords,
                              ket_first,
                              ket_last);

        for (int64_t j = bra_first; j < bra_last; j++) 
        {
            const auto bra_coord = bra_gto_coords[j];

            // compute primitive integrals block (X_XXX)

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

                    octurec::compPrimitiveOctupolePF_X_XXX(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXY)

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

                    octurec::compPrimitiveOctupolePF_X_XXY(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XXZ)

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

                    octurec::compPrimitiveOctupolePF_X_XXZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XYY)

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

                    octurec::compPrimitiveOctupolePF_X_XYY(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XYZ)

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

                    octurec::compPrimitiveOctupolePF_X_XYZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_XZZ)

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

                    octurec::compPrimitiveOctupolePF_X_XZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YYY)

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

                    octurec::compPrimitiveOctupolePF_X_YYY(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YYZ)

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

                    octurec::compPrimitiveOctupolePF_X_YYZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_YZZ)

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

                    octurec::compPrimitiveOctupolePF_X_YZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (X_ZZZ)

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

                    octurec::compPrimitiveOctupolePF_X_ZZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXX)

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

                    octurec::compPrimitiveOctupolePF_Y_XXX(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXY)

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

                    octurec::compPrimitiveOctupolePF_Y_XXY(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XXZ)

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

                    octurec::compPrimitiveOctupolePF_Y_XXZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XYY)

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

                    octurec::compPrimitiveOctupolePF_Y_XYY(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XYZ)

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

                    octurec::compPrimitiveOctupolePF_Y_XYZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_XZZ)

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

                    octurec::compPrimitiveOctupolePF_Y_XZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YYY)

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

                    octurec::compPrimitiveOctupolePF_Y_YYY(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YYZ)

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

                    octurec::compPrimitiveOctupolePF_Y_YYZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_YZZ)

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

                    octurec::compPrimitiveOctupolePF_Y_YZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Y_ZZZ)

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

                    octurec::compPrimitiveOctupolePF_Y_ZZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXX)

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

                    octurec::compPrimitiveOctupolePF_Z_XXX(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXY)

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

                    octurec::compPrimitiveOctupolePF_Z_XXY(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XXZ)

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

                    octurec::compPrimitiveOctupolePF_Z_XXZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XYY)

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

                    octurec::compPrimitiveOctupolePF_Z_XYY(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XYZ)

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

                    octurec::compPrimitiveOctupolePF_Z_XYZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_XZZ)

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

                    octurec::compPrimitiveOctupolePF_Z_XZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YYY)

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

                    octurec::compPrimitiveOctupolePF_Z_YYY(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YYZ)

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

                    octurec::compPrimitiveOctupolePF_Z_YYZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxx, buffer_xxx, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_YZZ)

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

                    octurec::compPrimitiveOctupolePF_Z_YZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (Z_ZZZ)

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

                    octurec::compPrimitiveOctupolePF_Z_ZZZ(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, ang_order);

        }
    }
}

} // octurec namespace

