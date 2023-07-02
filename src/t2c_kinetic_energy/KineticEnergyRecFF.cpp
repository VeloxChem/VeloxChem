#include "KineticEnergyRecFF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveKineticEnergyFF_XXX_T.hpp"
#include "PrimitiveKineticEnergyFF_XXY_T.hpp"
#include "PrimitiveKineticEnergyFF_XXZ_T.hpp"
#include "PrimitiveKineticEnergyFF_XYY_T.hpp"
#include "PrimitiveKineticEnergyFF_XYZ_T.hpp"
#include "PrimitiveKineticEnergyFF_XZZ_T.hpp"
#include "PrimitiveKineticEnergyFF_YYY_T.hpp"
#include "PrimitiveKineticEnergyFF_YYZ_T.hpp"
#include "PrimitiveKineticEnergyFF_YZZ_T.hpp"
#include "PrimitiveKineticEnergyFF_ZZZ_T.hpp"
#include "T2CDistributor.hpp"

namespace kinrec {  // kinrec namespace

auto
compKineticEnergyFF(CSubMatrix* matrix, const CGtoBlock& gto_block, const int64_t bra_first, const int64_t bra_last) -> void
{
    // spherical transformation factors

    const double f3_5 = std::sqrt(2.5);

    const double f3_15 = 2.0 * std::sqrt(15.0);

    const double f3_3 = std::sqrt(1.5);

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

    const auto nbatches = batch::getNumberOfBatches(ncgtos, simd_width);

    for (int64_t i = 0; i < nbatches; i++)
    {
        const auto [ket_first, ket_last] = batch::getBatchRange(i, ncgtos, simd_width);

        const auto ket_dim = ket_last - ket_first;

        simd::loadCoordinates(ket_coords_x, ket_coords_y, ket_coords_z, gto_coords, ket_first, ket_last);

        for (int64_t j = bra_first; j < bra_last; j++)
        {
            const auto bra_coord = gto_coords[j];

            // compute primitive integrals block (XXX)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    kinrec::compPrimitiveKineticEnergyFF_XXX_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, f3_5 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_5 * 3.0 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_5 * 3.0 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, f3_5 * f3_15, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, f3_5 * 4.0 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_5 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, f3_5 * 4.0 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, f3_5 * 2.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            // compute primitive integrals block (XXY)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    kinrec::compPrimitiveKineticEnergyFF_XXY_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f3_5 * f3_15, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, 3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, 3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * f3_5 * 2.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            // compute primitive integrals block (XXZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    kinrec::compPrimitiveKineticEnergyFF_XXZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f3_15 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, 0.5 * f3_15 * f3_15, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -3.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, 0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f3_15 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -3.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, 0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, 0.5 * f3_15 * 2.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            // compute primitive integrals block (XYY)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    kinrec::compPrimitiveKineticEnergyFF_XYY_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_5 * f3_15, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f3_5 * 2.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            // compute primitive integrals block (XYZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    kinrec::compPrimitiveKineticEnergyFF_XYZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -f3_15 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, f3_15 * f3_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_15 * 3.0 * f3_5, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_15 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_15 * 3.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_15 * 0.5 * f3_15, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_15 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_15 * 3.0 * f3_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, f3_15 * f3_15, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, f3_15 * 4.0 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_15 * f3_5, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_15 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_15 * 3.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_15 * 0.5 * f3_15, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, f3_15 * 4.0 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, f3_15 * 2.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            // compute primitive integrals block (XZZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    kinrec::compPrimitiveKineticEnergyFF_XZZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, 4.0 * f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            // compute primitive integrals block (YYY)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    kinrec::compPrimitiveKineticEnergyFF_YYY_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_5 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_5 * 3.0 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_5 * 3.0 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_5 * f3_15, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_5 * 4.0 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_5 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_5 * 4.0 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_5 * 2.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            // compute primitive integrals block (YYZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    kinrec::compPrimitiveKineticEnergyFF_YYZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f3_15 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -0.5 * f3_15 * f3_15, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -3.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f3_15 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -3.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -0.5 * f3_15 * 2.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            // compute primitive integrals block (YZZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    kinrec::compPrimitiveKineticEnergyFF_YZZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, 4.0 * f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZ)

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

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    kinrec::compPrimitiveKineticEnergyFF_ZZZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -2.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 2.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 2.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -2.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -2.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 2.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -2.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -2.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, 2.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, 2.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -2.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -2.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -2.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -2.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, 2.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, 2.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);
        }
    }
}

auto
compKineticEnergyFF(CSubMatrix*      matrix,
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

            // compute primitive integrals block (XXX)

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

                    kinrec::compPrimitiveKineticEnergyFF_XXX_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY)

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

                    kinrec::compPrimitiveKineticEnergyFF_XXY_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ)

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

                    kinrec::compPrimitiveKineticEnergyFF_XXZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY)

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

                    kinrec::compPrimitiveKineticEnergyFF_XYY_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ)

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

                    kinrec::compPrimitiveKineticEnergyFF_XYZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ)

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

                    kinrec::compPrimitiveKineticEnergyFF_XZZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY)

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

                    kinrec::compPrimitiveKineticEnergyFF_YYY_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ)

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

                    kinrec::compPrimitiveKineticEnergyFF_YYZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ)

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

                    kinrec::compPrimitiveKineticEnergyFF_YZZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xxz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xyy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yyz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ)

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

                    kinrec::compPrimitiveKineticEnergyFF_ZZZ_T(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);
        }
    }
}

}  // namespace kinrec
