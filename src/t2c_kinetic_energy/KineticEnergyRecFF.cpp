#include "KineticEnergyRecFF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "MathConst.hpp"
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

auto
compPrimitiveKineticEnergyFF_XXX_T(TDoubleArray&       buffer_xxx,
                                   TDoubleArray&       buffer_xxy,
                                   TDoubleArray&       buffer_xxz,
                                   TDoubleArray&       buffer_xyy,
                                   TDoubleArray&       buffer_xyz,
                                   TDoubleArray&       buffer_xzz,
                                   TDoubleArray&       buffer_yyy,
                                   TDoubleArray&       buffer_yyz,
                                   TDoubleArray&       buffer_yzz,
                                   TDoubleArray&       buffer_zzz,
                                   const double        bra_exp,
                                   const double        bra_norm,
                                   const TPoint3D&     bra_coord,
                                   const TDoubleArray& ket_exps,
                                   const TDoubleArray& ket_norms,
                                   const TDoubleArray& ket_coords_x,
                                   const TDoubleArray& ket_coords_y,
                                   const TDoubleArray& ket_coords_z,
                                   const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * fz_0 - (9.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * fz_0 -
                               (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 -
                               (9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * fz_0);

        fints_xxx[i] += fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * fz_0 + 15.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 +
                               45.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * fz_0 + 15.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * fz_0 -
                               (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_xxx[i] +=
            fss * (54.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 + 18.0 * fe_0 * fe_0 * rpa_x * rpa_x * fz_0 + 18.0 * fe_0 * fe_0 * rpb_x * rpb_x * fz_0 +
                   (45.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_x * fz_0);

        fints_xxx[i] += fss * 12.0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x +
                               (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x + (27.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_xxx[i] +=
            ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                               3.0 * fbe_0 * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * fz_0 +
                               15.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += fss * (30.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * fz_0 +
                               18.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 + 12.0 * fe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                               fke_0 * rpa_x * rpa_x * rpa_x * rpb_y * fz_0);

        fints_xxy[i] += fss * 12.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * fz_0;

        fints_xxy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x + 3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xxy[i] += ftt * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x;

        fints_xxz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                               3.0 * fbe_0 * rpa_x * rpb_z * rpb_x * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * fz_0 +
                               15.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += fss * (30.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 + 5.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * fz_0 +
                               18.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 + 12.0 * fe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                               fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * fz_0);

        fints_xxz[i] += fss * 12.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * fz_0;

        fints_xxz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x + 3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xxz[i] += ftt * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x;

        fints_xyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * fz_0);

        fints_xyy[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * fz_0 + 15.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 + 5.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_xyy[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * fz_0 + 6.0 * fe_0 * fe_0 * rpb_y * rpb_y * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_x * rpa_x * rpa_x * rpb_x * fz_0);

        fints_xyy[i] += fss * 12.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * fz_0;

        fints_xyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_xyy[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 + 15.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 +
                               6.0 * fe_0 * fe_0 * rpb_z * rpb_y * fz_0);

        fints_xyz[i] += fss * 12.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * fz_0;

        fints_xyz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y + rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * fz_0);

        fints_xzz[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * fz_0 + 15.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 + 5.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_xzz[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * fz_0 + 6.0 * fe_0 * fe_0 * rpb_z * rpb_z * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_x * rpa_x * rpa_x * rpb_x * fz_0);

        fints_xzz[i] += fss * 12.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * fz_0;

        fints_xzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_xzz[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_y * rpb_y * rpb_y * fz_0 -
                               (9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * fz_0 + 15.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * fz_0 +
                               15.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * fz_0);

        fints_yyy[i] += fss * (18.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 - 3.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_y * fz_0 +
                               12.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * fz_0 + 15.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * fz_0);

        fints_yyz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 - fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * fz_0 +
                               12.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * fz_0 + 15.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y * fz_0);

        fints_yzz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 - fke_0 * rpa_x * rpa_x * rpa_x * rpb_y * fz_0 +
                               12.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * fz_0 - 3.0 * fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 -
                               (9.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * fz_0 + 15.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 +
                               15.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * fz_0);

        fints_zzz[i] += fss * (18.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 - 3.0 * fke_0 * rpa_x * rpa_x * rpa_x * rpb_z * fz_0 +
                               12.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * fz_0);

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveKineticEnergyFF_XXY_T(TDoubleArray&       buffer_xxx,
                                   TDoubleArray&       buffer_xxy,
                                   TDoubleArray&       buffer_xxz,
                                   TDoubleArray&       buffer_xyy,
                                   TDoubleArray&       buffer_xyz,
                                   TDoubleArray&       buffer_xzz,
                                   TDoubleArray&       buffer_yyy,
                                   TDoubleArray&       buffer_yyz,
                                   TDoubleArray&       buffer_yzz,
                                   TDoubleArray&       buffer_zzz,
                                   const double        bra_exp,
                                   const double        bra_norm,
                                   const TPoint3D&     bra_coord,
                                   const TDoubleArray& ket_exps,
                                   const TDoubleArray& ket_norms,
                                   const TDoubleArray& ket_coords_x,
                                   const TDoubleArray& ket_coords_y,
                                   const TDoubleArray& ket_coords_z,
                                   const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * fz_0 - fbe_0 * rpa_y * rpb_x * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * fz_0 +
                               30.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += fss * (15.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * fz_0 +
                               12.0 * fe_0 * fe_0 * rpa_y * rpa_x * fz_0 + 18.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 -
                               3.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_x * fz_0);

        fints_xxx[i] += fss * 12.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxx[i] += ftt * (3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xxx[i] += ftt * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x;

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_y * rpb_y * rpb_x * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * fz_0);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * fz_0 + 20.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * fz_0 + 5.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xxy[i] +=
            fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0 + 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0 + 8.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 +
                   2.0 * fe_0 * fe_0 * rpa_x * rpa_x * fz_0 + 2.0 * fe_0 * fe_0 * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_y * rpa_x * rpa_x * rpb_y * fz_0 +
                               12.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += ftt * (2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xxy[i] += ftt * (fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * fz_0 + 20.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * fz_0);

        fints_xxz[i] += fss * (5.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 -
                               fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * fz_0 + 12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += ftt * (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z +
                               rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                               fbe_0 * rpa_y * rpb_y * rpb_y * rpb_x * fz_0 - fe_0 * fke_0 * rpa_y * rpa_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * fz_0);

        fints_xyy[i] += fss * (10.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * fz_0 + 5.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * fz_0 + 10.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_y * rpa_x * fz_0);

        fints_xyy[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 + 8.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 + 4.0 * fe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                   fke_0 * rpa_y * rpa_x * rpa_x * rpb_x * fz_0 + 12.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x + fe_0 * rpa_x * rpa_x * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x);

        fints_xyy[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + fe_0 * fe_0 * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x +
                               rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_x * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_y * rpb_x * fz_0 +
                               10.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 + 5.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * fz_0);

        fints_xyz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 + 2.0 * fe_0 * fe_0 * rpb_z * rpb_x * fz_0 +
                               12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xyz[i] += ftt * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x;

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_z * rpb_x * fz_0 -
                               fe_0 * fke_0 * rpa_y * rpa_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * fz_0 +
                               10.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * fz_0);

        fints_xzz[i] += fss * (5.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_y * rpa_x * fz_0 + 2.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 -
                               fke_0 * rpa_y * rpa_x * rpa_x * rpb_x * fz_0);

        fints_xzz[i] += fss * 12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * fz_0;

        fints_xzz[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xzz[i] += ftt * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x;

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_y * rpb_y * rpb_y * rpb_y * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * fz_0);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * fz_0 + 15.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * fz_0 + 15.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_yyy[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0 + 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * fz_0 + 6.0 * fe_0 * fe_0 * rpb_y * rpb_y * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_y * fz_0);

        fints_yyy[i] += fss * 12.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * fz_0;

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y +
                               (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_yyy[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_y * rpb_y * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * fz_0);

        fints_yyz[i] += fss * (5.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * fz_0 + 10.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 + 4.0 * fe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                               fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * fz_0);

        fints_yyz[i] += fss * 12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * fz_0;

        fints_yyz[i] +=
            ftt * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y +
                   fe_0 * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yyz[i] += ftt * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y;

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * fz_0);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * fz_0 + 5.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * fz_0 + 5.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_yzz[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0 + 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * fz_0 + 2.0 * fe_0 * fe_0 * rpb_z * rpb_z * fz_0 +
                   (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_y * rpa_x * rpa_x * rpb_y * fz_0);

        fints_yzz[i] += fss * 12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * fz_0;

        fints_yzz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_yzz[i] +=
            ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * fz_0 + 15.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * fz_0 +
                               5.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * fz_0);

        fints_zzz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 - 3.0 * fke_0 * rpa_y * rpa_x * rpa_x * rpb_z * fz_0 +
                               12.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * fz_0);

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveKineticEnergyFF_XXZ_T(TDoubleArray&       buffer_xxx,
                                   TDoubleArray&       buffer_xxy,
                                   TDoubleArray&       buffer_xxz,
                                   TDoubleArray&       buffer_xyy,
                                   TDoubleArray&       buffer_xyz,
                                   TDoubleArray&       buffer_xzz,
                                   TDoubleArray&       buffer_yyy,
                                   TDoubleArray&       buffer_yyz,
                                   TDoubleArray&       buffer_yzz,
                                   TDoubleArray&       buffer_zzz,
                                   const double        bra_exp,
                                   const double        bra_norm,
                                   const TPoint3D&     bra_coord,
                                   const TDoubleArray& ket_exps,
                                   const TDoubleArray& ket_norms,
                                   const TDoubleArray& ket_coords_x,
                                   const TDoubleArray& ket_coords_y,
                                   const TDoubleArray& ket_coords_z,
                                   const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * fz_0 - fbe_0 * rpa_z * rpb_x * rpb_x * rpb_x * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * fz_0 +
                               30.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += fss * (15.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * fz_0 +
                               12.0 * fe_0 * fe_0 * rpa_z * rpa_x * fz_0 + 18.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_x * fz_0);

        fints_xxx[i] += fss * 12.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxx[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xxx[i] += ftt * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x;

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * fz_0 - fbe_0 * rpa_z * rpb_y * rpb_x * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * fz_0 + 20.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * fz_0);

        fints_xxy[i] += fss * (5.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 -
                               fke_0 * rpa_z * rpa_x * rpa_x * rpb_y * fz_0 + 12.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y +
                               rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_z * rpb_z * rpb_x * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * fz_0);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * fz_0 + 20.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * fz_0 + 5.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x * fz_0);

        fints_xxz[i] +=
            fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 + 8.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 +
                   2.0 * fe_0 * fe_0 * rpa_x * rpa_x * fz_0 + 2.0 * fe_0 * fe_0 * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * fz_0 +
                               12.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_xxz[i] += ftt * (fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * fz_0 - fbe_0 * rpa_z * rpb_y * rpb_y * rpb_x * fz_0 -
                               fe_0 * fke_0 * rpa_z * rpa_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * fz_0);

        fints_xyy[i] += fss * (5.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_x * fz_0 + 2.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 -
                               fke_0 * rpa_z * rpa_x * rpa_x * rpb_x * fz_0);

        fints_xyy[i] += fss * 12.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * fz_0;

        fints_xyy[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xyy[i] += ftt * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x;

        fints_xyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_x * fz_0 - fbe_0 * rpa_z * rpb_z * rpb_y * rpb_x * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * fz_0 + 5.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 + 2.0 * fe_0 * fe_0 * rpb_y * rpb_x * fz_0 +
                               12.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xyz[i] += ftt * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x;

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * fz_0);

        fints_xzz[i] += fss * (10.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * fz_0 + 5.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * fz_0 + 10.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_x * fz_0);

        fints_xzz[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 + 8.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 + 4.0 * fe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                   fke_0 * rpa_z * rpa_x * rpa_x * rpb_x * fz_0 + 12.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xzz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x + fe_0 * rpa_x * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x);

        fints_xzz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + fe_0 * fe_0 * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x +
                               rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * fz_0 - fbe_0 * rpa_z * rpb_y * rpb_y * rpb_y * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * fz_0 + 15.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 - 3.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_y * fz_0 +
                               12.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_z * rpb_z * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * fz_0);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * fz_0 +
                               5.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * fz_0 + 5.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_yyz[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 + 2.0 * fe_0 * fe_0 * rpa_x * rpa_x * fz_0 + 2.0 * fe_0 * fe_0 * rpb_y * rpb_y * fz_0 +
                   (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * fz_0);

        fints_yyz[i] += fss * 12.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * fz_0;

        fints_yyz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_yyz[i] +=
            ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_y * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y * fz_0);

        fints_yzz[i] += fss * (5.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * fz_0 + 10.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 + 4.0 * fe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                               fke_0 * rpa_z * rpa_x * rpa_x * rpb_y * fz_0);

        fints_yzz[i] += fss * 12.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * fz_0;

        fints_yzz[i] +=
            ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y +
                   fe_0 * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yzz[i] += ftt * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y;

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * fz_0);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpa_x * fz_0 + 15.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * fz_0 +
                               5.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * fz_0 + 15.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_zzz[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 + 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * fz_0 + 6.0 * fe_0 * fe_0 * rpb_z * rpb_z * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fke_0 * rpa_z * rpa_x * rpa_x * rpb_z * fz_0);

        fints_zzz[i] += fss * 12.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * fz_0;

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z +
                               (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_zzz[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveKineticEnergyFF_XYY_T(TDoubleArray&       buffer_xxx,
                                   TDoubleArray&       buffer_xxy,
                                   TDoubleArray&       buffer_xxz,
                                   TDoubleArray&       buffer_xyy,
                                   TDoubleArray&       buffer_xyz,
                                   TDoubleArray&       buffer_xzz,
                                   TDoubleArray&       buffer_yyy,
                                   TDoubleArray&       buffer_yyz,
                                   TDoubleArray&       buffer_yzz,
                                   TDoubleArray&       buffer_zzz,
                                   const double        bra_exp,
                                   const double        bra_norm,
                                   const TPoint3D&     bra_coord,
                                   const TDoubleArray& ket_exps,
                                   const TDoubleArray& ket_norms,
                                   const TDoubleArray& ket_coords_x,
                                   const TDoubleArray& ket_coords_y,
                                   const TDoubleArray& ket_coords_z,
                                   const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * fz_0);

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * fz_0 + 15.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_xxx[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_y * rpa_y * fz_0 + 6.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * rpb_x * rpb_x * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_x * fz_0);

        fints_xxx[i] += fss * 12.0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xxx[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * fz_0 - fbe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                               fbe_0 * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_y * rpa_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * fz_0);

        fints_xxy[i] += fss * (10.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * fz_0 +
                               10.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_y * rpa_x * fz_0);

        fints_xxy[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 + 4.0 * fe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                   fke_0 * rpa_y * rpa_y * rpa_x * rpb_y * fz_0 + 12.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] +=
            ftt * (fe_0 * rpa_y * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y + fe_0 * rpa_y * rpa_y * rpb_y * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x);

        fints_xxy[i] += ftt * (fe_0 * fe_0 * rpa_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x +
                               rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_x * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * fz_0);

        fints_xxz[i] += fss * (10.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * fz_0 + 5.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 + 4.0 * fe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                               fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * fz_0);

        fints_xxz[i] += fss * 12.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * fz_0;

        fints_xxz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z + fe_0 * rpa_y * rpa_y * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xxz[i] += ftt * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x;

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * fz_0);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * fz_0 + 20.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] +=
            fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0 + 8.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0 + 2.0 * fe_0 * fe_0 * rpa_y * rpa_y * fz_0 +
                   6.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * rpb_y * rpb_y * fz_0);

        fints_xyy[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_y * rpa_y * rpa_x * rpb_x * fz_0 +
                               12.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += ftt * (2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x +
                               fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xyy[i] +=
            ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y +
                   (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_y * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 +
                               10.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 + 2.0 * fe_0 * fe_0 * rpb_z * rpb_y * fz_0 +
                               12.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_xyz[i] += ftt * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x;

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * fz_0);

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 + 5.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_xzz[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_y * rpa_y * fz_0 + 2.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * rpb_z * rpb_z * fz_0 +
                   (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_y * rpa_y * rpa_x * rpb_x * fz_0);

        fints_xzz[i] += fss * 12.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * fz_0;

        fints_xzz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xzz[i] +=
            ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * fz_0 - fbe_0 * rpa_x * rpb_y * rpb_y * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_y * rpa_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * fz_0 +
                               30.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += fss * (15.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * fz_0 + 5.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * fz_0 +
                               12.0 * fe_0 * fe_0 * rpa_y * rpa_x * fz_0 + 18.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 -
                               3.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_y * fz_0);

        fints_yyy[i] += fss * 12.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * fz_0;

        fints_yyy[i] += ftt * (3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y);

        fints_yyy[i] += ftt * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y;

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * fz_0 + 20.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * fz_0);

        fints_yyz[i] += fss * (5.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 -
                               fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * fz_0 + 12.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += ftt * (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z +
                               rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 -
                               fe_0 * fke_0 * rpa_y * rpa_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * fz_0 +
                               10.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * fz_0);

        fints_yzz[i] += fss * (5.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y * fz_0 + 5.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_y * rpa_x * fz_0 + 2.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 -
                               fke_0 * rpa_y * rpa_y * rpa_x * rpb_y * fz_0);

        fints_yzz[i] += fss * 12.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * fz_0;

        fints_yzz[i] += ftt * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y);

        fints_yzz[i] += ftt * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y;

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * fz_0 + 15.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * fz_0 +
                               5.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * fz_0);

        fints_zzz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 - 3.0 * fke_0 * rpa_y * rpa_y * rpa_x * rpb_z * fz_0 +
                               12.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * fz_0);

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveKineticEnergyFF_XYZ_T(TDoubleArray&       buffer_xxx,
                                   TDoubleArray&       buffer_xxy,
                                   TDoubleArray&       buffer_xxz,
                                   TDoubleArray&       buffer_xyy,
                                   TDoubleArray&       buffer_xyz,
                                   TDoubleArray&       buffer_xzz,
                                   TDoubleArray&       buffer_yyy,
                                   TDoubleArray&       buffer_yyz,
                                   TDoubleArray&       buffer_yzz,
                                   TDoubleArray&       buffer_zzz,
                                   const double        bra_exp,
                                   const double        bra_norm,
                                   const TPoint3D&     bra_coord,
                                   const TDoubleArray& ket_exps,
                                   const TDoubleArray& ket_norms,
                                   const TDoubleArray& ket_coords_x,
                                   const TDoubleArray& ket_coords_y,
                                   const TDoubleArray& ket_coords_z,
                                   const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_y * fz_0 + 15.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_x * fz_0);

        fints_xxx[i] += fss * 12.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y + rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpa_x * fz_0);

        fints_xxy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 - fke_0 * rpa_z * rpa_y * rpa_x * rpb_y * fz_0 +
                               12.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y + fe_0 * rpa_z * rpa_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xxy[i] += ftt * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x;

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_y * rpa_x * fz_0);

        fints_xxz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 - fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * fz_0 +
                               12.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z + fe_0 * rpa_z * rpa_y * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xxz[i] += ftt * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x;

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_y * fz_0 + 5.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * fz_0 + 10.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpa_y * fz_0);

        fints_xyy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 - fke_0 * rpa_z * rpa_y * rpa_x * rpb_x * fz_0 +
                               12.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] +=
            ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y +
                   fe_0 * rpa_z * rpa_x * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y);

        fints_xyy[i] += ftt * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x;

        fints_xyz[i] += fss * (5.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * fz_0 + 5.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0);

        fints_xyz[i] += fss * (2.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 + (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 +
                               12.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xyz[i] +=
            ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_y * fz_0 + 5.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * fz_0 + 10.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpa_y * fz_0);

        fints_xzz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 - fke_0 * rpa_z * rpa_y * rpa_x * rpb_x * fz_0 +
                               12.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xzz[i] +=
            ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z +
                   fe_0 * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_xzz[i] += ftt * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x;

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_x * fz_0 + 15.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * fz_0 +
                               15.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * fz_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpa_x * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_y * fz_0);

        fints_yyy[i] += fss * 12.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * fz_0;

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x + rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * fz_0 + 5.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_y * rpa_x * fz_0);

        fints_yyz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 - fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * fz_0 +
                               12.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z + fe_0 * rpa_z * rpa_x * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y);

        fints_yyz[i] += ftt * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y;

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * fz_0 + 10.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpa_x * fz_0);

        fints_yzz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 - fke_0 * rpa_z * rpa_y * rpa_x * rpb_y * fz_0 +
                               12.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yzz[i] +=
            ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z +
                   fe_0 * rpa_y * rpa_x * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z);

        fints_yzz[i] += ftt * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y;

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_x * fz_0 + 15.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * fz_0 +
                               15.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * fz_0 + 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_y * rpa_x * rpb_z * fz_0);

        fints_zzz[i] += fss * 12.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * fz_0;

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x + rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveKineticEnergyFF_XZZ_T(TDoubleArray&       buffer_xxx,
                                   TDoubleArray&       buffer_xxy,
                                   TDoubleArray&       buffer_xxz,
                                   TDoubleArray&       buffer_xyy,
                                   TDoubleArray&       buffer_xyz,
                                   TDoubleArray&       buffer_xzz,
                                   TDoubleArray&       buffer_yyy,
                                   TDoubleArray&       buffer_yyz,
                                   TDoubleArray&       buffer_yzz,
                                   TDoubleArray&       buffer_zzz,
                                   const double        bra_exp,
                                   const double        bra_norm,
                                   const TPoint3D&     bra_coord,
                                   const TDoubleArray& ket_exps,
                                   const TDoubleArray& ket_norms,
                                   const TDoubleArray& ket_coords_x,
                                   const TDoubleArray& ket_coords_y,
                                   const TDoubleArray& ket_coords_z,
                                   const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * fz_0);

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * fz_0 + 15.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_xxx[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_z * fz_0 + 6.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * rpb_x * rpb_x * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_x * fz_0);

        fints_xxx[i] += fss * 12.0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * fz_0;

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xxx[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * fz_0 - fbe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                               fbe_0 * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * fz_0);

        fints_xxy[i] += fss * (10.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 + 4.0 * fe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_x * rpb_y * fz_0);

        fints_xxy[i] += fss * 12.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * fz_0;

        fints_xxy[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y + fe_0 * rpa_z * rpa_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xxy[i] += ftt * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x;

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                               fbe_0 * rpa_x * rpb_z * rpb_x * rpb_x * fz_0 - fe_0 * fke_0 * rpa_z * rpa_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * fz_0);

        fints_xxz[i] += fss * (10.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * fz_0 + 5.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_x * fz_0);

        fints_xxz[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 + 4.0 * fe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                   fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * fz_0 + 12.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxz[i] +=
            ftt * (fe_0 * rpa_z * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z + fe_0 * rpa_z * rpa_z * rpb_z * rpb_x +
                   (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x);

        fints_xxz[i] += ftt * (fe_0 * fe_0 * rpa_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x +
                               rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * fz_0);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * fz_0 + 5.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_xyy[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * fz_0 + 2.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * rpb_y * rpb_y * fz_0 +
                   (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_z * rpa_z * rpa_x * rpb_x * fz_0);

        fints_xyy[i] += fss * 12.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * fz_0;

        fints_xyy[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xyy[i] +=
            ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_y * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_y * rpb_x * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 + 2.0 * fe_0 * fe_0 * rpb_z * rpb_y * fz_0 +
                               12.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_xyz[i] += ftt * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x;

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_x * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * fz_0);

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_x * fz_0 + 20.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * fz_0 +
                               5.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xzz[i] +=
            fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0 + 8.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 + 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * fz_0 +
                   6.0 * fe_0 * fe_0 * rpa_x * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * rpb_z * rpb_z * fz_0);

        fints_xzz[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_z * rpa_z * rpa_x * rpb_x * fz_0 +
                               12.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xzz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x +
                               fe_0 * fe_0 * rpa_z * rpb_z);

        fints_xzz[i] +=
            ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z +
                   (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * fz_0 - fbe_0 * rpa_x * rpb_y * rpb_y * rpb_y * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * fz_0 + 15.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += fss * (6.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 - 3.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_y * fz_0 +
                               12.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 -
                               fe_0 * fke_0 * rpa_z * rpa_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += fss * (5.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * fz_0 + 5.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_x * fz_0 + 2.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * fz_0);

        fints_yyz[i] += fss * 12.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * fz_0;

        fints_yyz[i] += ftt * (fe_0 * rpa_z * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z);

        fints_yyz[i] += ftt * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y;

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_y * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_y * fz_0 + 20.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * fz_0);

        fints_yzz[i] += fss * (5.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y * fz_0 + 6.0 * fe_0 * fe_0 * rpa_x * rpb_y * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_x * rpb_y * fz_0 + 12.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yzz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y +
                               rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_x * rpb_z * fz_0 - fbe_0 * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_x * rpb_z * fz_0 +
                               30.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * fz_0);

        fints_zzz[i] += fss * (15.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * fz_0 + 5.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z * fz_0 +
                               12.0 * fe_0 * fe_0 * rpa_z * rpa_x * fz_0 + 18.0 * fe_0 * fe_0 * rpa_x * rpb_z * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_z * rpa_x * rpb_z * fz_0);

        fints_zzz[i] += fss * 12.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * fz_0;

        fints_zzz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z);

        fints_zzz[i] += ftt * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z;
    }
}

auto
compPrimitiveKineticEnergyFF_YYY_T(TDoubleArray&       buffer_xxx,
                                   TDoubleArray&       buffer_xxy,
                                   TDoubleArray&       buffer_xxz,
                                   TDoubleArray&       buffer_xyy,
                                   TDoubleArray&       buffer_xyz,
                                   TDoubleArray&       buffer_xzz,
                                   TDoubleArray&       buffer_yyy,
                                   TDoubleArray&       buffer_yyz,
                                   TDoubleArray&       buffer_yzz,
                                   TDoubleArray&       buffer_zzz,
                                   const double        bra_exp,
                                   const double        bra_norm,
                                   const TPoint3D&     bra_coord,
                                   const TDoubleArray& ket_exps,
                                   const TDoubleArray& ket_norms,
                                   const TDoubleArray& ket_coords_x,
                                   const TDoubleArray& ket_coords_y,
                                   const TDoubleArray& ket_coords_z,
                                   const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_x * rpb_x * rpb_x * fz_0 -
                               (9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * fz_0 + 15.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * fz_0);

        fints_xxx[i] += fss * (18.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 - 3.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_x * fz_0 +
                               12.0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_y * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * fz_0);

        fints_xxy[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * fz_0 + 15.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_xxy[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0 + 6.0 * fe_0 * fe_0 * rpa_y * rpa_y * fz_0 + 6.0 * fe_0 * fe_0 * rpb_x * rpb_x * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_y * rpa_y * rpa_y * rpb_y * fz_0);

        fints_xxy[i] += fss * 12.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * fz_0;

        fints_xxy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_xxy[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * fz_0 + 15.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * fz_0);

        fints_xxz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 - fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * fz_0 +
                               12.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                               3.0 * fbe_0 * rpa_y * rpb_y * rpb_y * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += fss * (30.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * fz_0 +
                               18.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 + 12.0 * fe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                               fke_0 * rpa_y * rpa_y * rpa_y * rpb_x * fz_0);

        fints_xyy[i] += fss * 12.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * fz_0;

        fints_xyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xyy[i] += ftt * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x;

        fints_xyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_y * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * fz_0 + 15.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * fz_0 +
                               6.0 * fe_0 * fe_0 * rpb_z * rpb_x * fz_0);

        fints_xyz[i] += fss * 12.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * fz_0;

        fints_xyz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x + rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_z * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * fz_0 + 15.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x * fz_0);

        fints_xzz[i] += fss * (6.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 - fke_0 * rpa_y * rpa_y * rpa_y * rpb_x * fz_0 +
                               12.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * fz_0 - (9.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * fz_0 -
                               (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_y * rpb_y * rpb_y * fz_0 -
                               (9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * fz_0);

        fints_yyy[i] += fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * fz_0 + 15.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * fz_0 +
                               45.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * fz_0 + 15.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * fz_0 -
                               (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_yyy[i] +=
            fss * (54.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0 + 18.0 * fe_0 * fe_0 * rpa_y * rpa_y * fz_0 + 18.0 * fe_0 * fe_0 * rpb_y * rpb_y * fz_0 +
                   (45.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_y * fz_0);

        fints_yyy[i] += fss * 12.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * fz_0;

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y +
                               (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y + (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_yyy[i] +=
            ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                               3.0 * fbe_0 * rpa_y * rpb_z * rpb_y * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * fz_0 +
                               15.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += fss * (30.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 + 5.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * fz_0 +
                               18.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 + 12.0 * fe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                               fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * fz_0);

        fints_yyz[i] += fss * 12.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * fz_0;

        fints_yyz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yyz[i] += ftt * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y;

        fints_yzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_z * rpb_y * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * fz_0);

        fints_yzz[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * fz_0 + 15.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * fz_0 +
                               15.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 + 5.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_yzz[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0 + 6.0 * fe_0 * fe_0 * rpa_y * rpa_y * fz_0 + 6.0 * fe_0 * fe_0 * rpb_z * rpb_z * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_y * rpa_y * rpa_y * rpb_y * fz_0);

        fints_yzz[i] += fss * 12.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * fz_0;

        fints_yzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_yzz[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * fz_0 - 3.0 * fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * fz_0 -
                               (9.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * fz_0 + 15.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * fz_0 +
                               15.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * fz_0);

        fints_zzz[i] += fss * (18.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 - 3.0 * fke_0 * rpa_y * rpa_y * rpa_y * rpb_z * fz_0 +
                               12.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * fz_0);

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveKineticEnergyFF_YYZ_T(TDoubleArray&       buffer_xxx,
                                   TDoubleArray&       buffer_xxy,
                                   TDoubleArray&       buffer_xxz,
                                   TDoubleArray&       buffer_xyy,
                                   TDoubleArray&       buffer_xyz,
                                   TDoubleArray&       buffer_xzz,
                                   TDoubleArray&       buffer_yyy,
                                   TDoubleArray&       buffer_yyz,
                                   TDoubleArray&       buffer_yzz,
                                   TDoubleArray&       buffer_zzz,
                                   const double        bra_exp,
                                   const double        bra_norm,
                                   const TPoint3D&     bra_coord,
                                   const TDoubleArray& ket_exps,
                                   const TDoubleArray& ket_norms,
                                   const TDoubleArray& ket_coords_x,
                                   const TDoubleArray& ket_coords_y,
                                   const TDoubleArray& ket_coords_z,
                                   const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * fz_0 - fbe_0 * rpa_z * rpb_x * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * fz_0 + 15.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 - 3.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_x * fz_0 +
                               12.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * fz_0 - fbe_0 * rpa_z * rpb_y * rpb_x * rpb_x * fz_0 -
                               fe_0 * fke_0 * rpa_z * rpa_y * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += fss * (5.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * fz_0 + 5.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_y * fz_0 + 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 -
                               fke_0 * rpa_z * rpa_y * rpa_y * rpb_y * fz_0);

        fints_xxy[i] += fss * 12.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * fz_0;

        fints_xxy[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y);

        fints_xxy[i] += ftt * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x;

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_z * rpb_z * rpb_x * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * fz_0);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * fz_0 + 5.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * fz_0 +
                               5.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_xxz[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 + 2.0 * fe_0 * fe_0 * rpa_y * rpa_y * fz_0 + 2.0 * fe_0 * fe_0 * rpb_x * rpb_x * fz_0 +
                   (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * fz_0);

        fints_xxz[i] += fss * 12.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * fz_0;

        fints_xxz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_xxz[i] +=
            ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * fz_0 - fbe_0 * rpa_z * rpb_y * rpb_y * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * fz_0 + 20.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * fz_0);

        fints_xyy[i] += fss * (5.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 -
                               fke_0 * rpa_z * rpa_y * rpa_y * rpb_x * fz_0 + 12.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x +
                               rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_x * fz_0 - fbe_0 * rpa_z * rpb_z * rpb_y * rpb_x * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * fz_0 + 5.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * rpb_y * rpb_x * fz_0 +
                               12.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xyz[i] += ftt * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x;

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x * fz_0);

        fints_xzz[i] += fss * (5.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * fz_0 + 10.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 + 4.0 * fe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                               fke_0 * rpa_z * rpa_y * rpa_y * rpb_x * fz_0);

        fints_xzz[i] += fss * 12.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * fz_0;

        fints_xzz[i] +=
            ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x +
                   fe_0 * rpa_y * rpa_y * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xzz[i] += ftt * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x;

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * fz_0 - fbe_0 * rpa_z * rpb_y * rpb_y * rpb_y * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * fz_0 +
                               30.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += fss * (15.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * fz_0 + 5.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * fz_0 +
                               12.0 * fe_0 * fe_0 * rpa_z * rpa_y * fz_0 + 18.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_y * fz_0);

        fints_yyy[i] += fss * 12.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * fz_0;

        fints_yyy[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y);

        fints_yyy[i] += ftt * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y;

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_z * rpb_z * rpb_y * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * fz_0);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * fz_0 + 20.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * fz_0 + 5.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * fz_0);

        fints_yyz[i] +=
            fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 + 8.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0 +
                   2.0 * fe_0 * fe_0 * rpa_y * rpa_y * fz_0 + 2.0 * fe_0 * fe_0 * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * fz_0 +
                               12.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * fz_0);

        fints_yyz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_yyz[i] += ftt * (fe_0 * fe_0 * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y +
                               (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                               fbe_0 * rpa_z * rpb_z * rpb_z * rpb_y * fz_0 - fe_0 * fke_0 * rpa_z * rpa_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * fz_0);

        fints_yzz[i] += fss * (10.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * fz_0 + 5.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * fz_0 + 10.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_y * fz_0);

        fints_yzz[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 + 8.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 + 4.0 * fe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                   fke_0 * rpa_z * rpa_y * rpa_y * rpb_y * fz_0 + 12.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yzz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y + fe_0 * rpa_y * rpa_y * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y);

        fints_yzz[i] += ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + fe_0 * fe_0 * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y +
                               rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * fz_0);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpa_y * fz_0 + 15.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * fz_0 +
                               5.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * fz_0 + 15.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_zzz[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 + 6.0 * fe_0 * fe_0 * rpa_y * rpa_y * fz_0 + 6.0 * fe_0 * fe_0 * rpb_z * rpb_z * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fke_0 * rpa_z * rpa_y * rpa_y * rpb_z * fz_0);

        fints_zzz[i] += fss * 12.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * fz_0;

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z +
                               (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_zzz[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z);
    }
}

auto
compPrimitiveKineticEnergyFF_YZZ_T(TDoubleArray&       buffer_xxx,
                                   TDoubleArray&       buffer_xxy,
                                   TDoubleArray&       buffer_xxz,
                                   TDoubleArray&       buffer_xyy,
                                   TDoubleArray&       buffer_xyz,
                                   TDoubleArray&       buffer_xzz,
                                   TDoubleArray&       buffer_yyy,
                                   TDoubleArray&       buffer_yyz,
                                   TDoubleArray&       buffer_yzz,
                                   TDoubleArray&       buffer_zzz,
                                   const double        bra_exp,
                                   const double        bra_norm,
                                   const TPoint3D&     bra_coord,
                                   const TDoubleArray& ket_exps,
                                   const TDoubleArray& ket_norms,
                                   const TDoubleArray& ket_coords_x,
                                   const TDoubleArray& ket_coords_y,
                                   const TDoubleArray& ket_coords_z,
                                   const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * fz_0 - fbe_0 * rpa_y * rpb_x * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * fz_0 + 15.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += fss * (6.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 - 3.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_x * fz_0 +
                               12.0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_y * rpb_y * rpb_x * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * fz_0);

        fints_xxy[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x * fz_0 -
                               (1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_xxy[i] +=
            fss * (2.0 * fe_0 * fe_0 * rpa_z * rpa_z * fz_0 + 2.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0 + 2.0 * fe_0 * fe_0 * rpb_x * rpb_x * fz_0 +
                   (3.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_z * rpa_z * rpa_y * rpb_y * fz_0);

        fints_xxy[i] += fss * 12.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * fz_0;

        fints_xxy[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xxy[i] +=
            ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 -
                               fe_0 * fke_0 * rpa_z * rpa_y * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x * fz_0);

        fints_xxz[i] += fss * (5.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * fz_0 + 5.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_y * fz_0 + 2.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * fz_0);

        fints_xxz[i] += fss * 12.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * fz_0;

        fints_xxz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_xxz[i] += ftt * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x;

        fints_xyy[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * fz_0 - fbe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                               fbe_0 * rpa_y * rpb_y * rpb_y * rpb_x * fz_0 - (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * fz_0);

        fints_xyy[i] += fss * (10.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x * fz_0 +
                               2.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 + 4.0 * fe_0 * fe_0 * rpb_y * rpb_x * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_y * rpb_x * fz_0);

        fints_xyy[i] += fss * 12.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * fz_0;

        fints_xyy[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x + fe_0 * rpa_z * rpa_z * rpb_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xyy[i] += ftt * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x;

        fints_xyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_x * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_y * rpb_x * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 + 2.0 * fe_0 * fe_0 * rpb_z * rpb_x * fz_0 +
                               12.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += ftt * (fe_0 * rpa_z * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x +
                               (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xyz[i] += ftt * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x;

        fints_xzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_x * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_z * rpb_x * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_x * fz_0 + 20.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * fz_0);

        fints_xzz[i] += fss * (5.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x * fz_0 + 6.0 * fe_0 * fe_0 * rpa_y * rpb_x * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_y * rpb_x * fz_0 + 12.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xzz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x +
                               rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_y * rpb_y * rpb_y * rpb_y * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * fz_0);

        fints_yyy[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * fz_0 + 15.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * fz_0 +
                               15.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * fz_0 + 5.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_yyy[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_z * rpa_z * fz_0 + 6.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0 + 6.0 * fe_0 * fe_0 * rpb_y * rpb_y * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_y * fz_0);

        fints_yyy[i] += fss * 12.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * fz_0;

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_yyy[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * fz_0 - fbe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                               fbe_0 * rpa_y * rpb_z * rpb_y * rpb_y * fz_0 - fe_0 * fke_0 * rpa_z * rpa_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * fz_0);

        fints_yyz[i] += fss * (10.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * fz_0 +
                               10.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * fz_0 + 5.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y * fz_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpa_y * fz_0);

        fints_yyz[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 + 2.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 + 4.0 * fe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                   fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * fz_0 + 12.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * fz_0);

        fints_yyz[i] +=
            ftt * (fe_0 * rpa_z * rpa_y * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z + fe_0 * rpa_z * rpa_z * rpb_z * rpb_y +
                   (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y);

        fints_yyz[i] += ftt * (fe_0 * fe_0 * rpa_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y +
                               rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_y * fz_0 - (1.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * fz_0 -
                               (1.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_z * rpb_y * fz_0 -
                               (1.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * fz_0);

        fints_yzz[i] += fss * (-(1.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_y * fz_0 + 20.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * fz_0 +
                               5.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yzz[i] +=
            fss * (-(1.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0 + 8.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 + 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * fz_0 +
                   6.0 * fe_0 * fe_0 * rpa_y * rpb_y * fz_0 + 2.0 * fe_0 * fe_0 * rpb_z * rpb_z * fz_0);

        fints_yzz[i] += fss * ((9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_z * rpa_z * rpa_y * rpb_y * fz_0 +
                               12.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yzz[i] += ftt * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y +
                               fe_0 * fe_0 * rpa_z * rpb_z);

        fints_yzz[i] +=
            ftt * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z +
                   (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_y * rpb_z * fz_0 - fbe_0 * rpa_y * rpb_z * rpb_z * rpb_z * fz_0 -
                               3.0 * fe_0 * fke_0 * rpa_z * rpa_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_y * rpb_z * fz_0 +
                               30.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * fz_0);

        fints_zzz[i] += fss * (15.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * fz_0 + 5.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z * fz_0 +
                               12.0 * fe_0 * fe_0 * rpa_z * rpa_y * fz_0 + 18.0 * fe_0 * fe_0 * rpa_y * rpb_z * fz_0 -
                               3.0 * fke_0 * rpa_z * rpa_z * rpa_y * rpb_z * fz_0);

        fints_zzz[i] += fss * 12.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * fz_0;

        fints_zzz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z +
                               (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_zzz[i] += ftt * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z;
    }
}

auto
compPrimitiveKineticEnergyFF_ZZZ_T(TDoubleArray&       buffer_xxx,
                                   TDoubleArray&       buffer_xxy,
                                   TDoubleArray&       buffer_xxz,
                                   TDoubleArray&       buffer_xyy,
                                   TDoubleArray&       buffer_xyz,
                                   TDoubleArray&       buffer_xzz,
                                   TDoubleArray&       buffer_yyy,
                                   TDoubleArray&       buffer_yyz,
                                   TDoubleArray&       buffer_yzz,
                                   TDoubleArray&       buffer_zzz,
                                   const double        bra_exp,
                                   const double        bra_norm,
                                   const TPoint3D&     bra_coord,
                                   const TDoubleArray& ket_exps,
                                   const TDoubleArray& ket_norms,
                                   const TDoubleArray& ket_coords_x,
                                   const TDoubleArray& ket_coords_y,
                                   const TDoubleArray& ket_coords_z,
                                   const int64_t       ket_dim) -> void
{
    // set up math constants

    const auto fpi = mathconst::getPiValue();

    // set up coordinates for bra side

    const auto bra_rx = bra_coord[0];

    const auto bra_ry = bra_coord[1];

    const auto bra_rz = bra_coord[2];

    // set up coordinates for ket side

    auto ket_rx = ket_coords_x.data();

    auto ket_ry = ket_coords_y.data();

    auto ket_rz = ket_coords_z.data();

    // set exponents and normalization factors on ket side

    auto ket_fe = ket_exps.data();

    auto ket_fn = ket_norms.data();

    // set up pointer to integrals buffer(s)

    auto fints_xxx = buffer_xxx.data();

    auto fints_xxy = buffer_xxy.data();

    auto fints_xxz = buffer_xxz.data();

    auto fints_xyy = buffer_xyy.data();

    auto fints_xyz = buffer_xyz.data();

    auto fints_xzz = buffer_xzz.data();

    auto fints_yyy = buffer_yyy.data();

    auto fints_yyz = buffer_yyz.data();

    auto fints_yzz = buffer_yzz.data();

    auto fints_zzz = buffer_zzz.data();

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
                             ket_fe,    \
                             ket_fn,    \
                             ket_rx,    \
                             ket_ry,    \
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto fke_0 = 1.0 / ket_fe[i];

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_x * rpb_x * rpb_x * fz_0 -
                               (9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * fz_0 + 15.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * fz_0);

        fints_xxx[i] += fss * (18.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 - 3.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_x * fz_0 +
                               12.0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x * fz_0);

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_y * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * fz_0 + 15.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * fz_0);

        fints_xxy[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 - fke_0 * rpa_z * rpa_z * rpa_z * rpb_y * fz_0 +
                               12.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * fz_0);

        fints_xxy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_x * rpb_x * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_z * rpb_x * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * fz_0);

        fints_xxz[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * fz_0 + 15.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_xxz[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * fz_0 + 6.0 * fe_0 * fe_0 * rpb_x * rpb_x * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * fz_0);

        fints_xxz[i] += fss * 12.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x * fz_0;

        fints_xxz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z);

        fints_xxz[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_y * rpb_y * rpb_x * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * fz_0 + 15.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x * fz_0 +
                               5.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * fz_0);

        fints_xyy[i] += fss * (6.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 - fke_0 * rpa_z * rpa_z * rpa_z * rpb_x * fz_0 +
                               12.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * fz_0);

        fints_xyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_x * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_z * rpb_y * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x * fz_0 + 15.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * fz_0 +
                               6.0 * fe_0 * fe_0 * rpb_y * rpb_x * fz_0);

        fints_xyz[i] += fss * 12.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * fz_0;

        fints_xyz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x + rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x);

        fints_xzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_x * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                               3.0 * fbe_0 * rpa_z * rpb_z * rpb_z * rpb_x * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_x * fz_0 +
                               15.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x * fz_0);

        fints_xzz[i] += fss * (30.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * fz_0 +
                               18.0 * fe_0 * fe_0 * rpa_z * rpb_x * fz_0 + 12.0 * fe_0 * fe_0 * rpb_z * rpb_x * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_z * rpb_x * fz_0);

        fints_xzz[i] += fss * 12.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x * fz_0;

        fints_xzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x + 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xzz[i] += ftt * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x;

        fints_yyy[i] += fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_y * rpb_y * rpb_y * fz_0 -
                               (9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * fz_0 + 15.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y * fz_0 +
                               15.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * fz_0);

        fints_yyy[i] += fss * (18.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 - 3.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_y * fz_0 +
                               12.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * fz_0);

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * fz_0 - (3.0 / 2.0) * fbe_0 * fe_0 * rpb_y * rpb_y * fz_0 -
                               (3.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_z * rpb_y * rpb_y * fz_0 -
                               (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * fz_0);

        fints_yyz[i] += fss * (-(3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * fz_0 + 15.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y * fz_0 +
                               15.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * fz_0 -
                               (3.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_yyz[i] +=
            fss * (6.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 + 6.0 * fe_0 * fe_0 * rpa_z * rpa_z * fz_0 + 6.0 * fe_0 * fe_0 * rpb_y * rpb_y * fz_0 +
                   (9.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * fz_0);

        fints_yyz[i] += fss * 12.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * fz_0;

        fints_yyz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z +
                               (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z);

        fints_yyz[i] +=
            ftt * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * (-(3.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_y * fz_0 - 3.0 * fbe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                               3.0 * fbe_0 * rpa_z * rpb_z * rpb_z * rpb_y * fz_0 - (3.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_y * fz_0 +
                               15.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y * fz_0);

        fints_yzz[i] += fss * (30.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * fz_0 + 5.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * fz_0 +
                               18.0 * fe_0 * fe_0 * rpa_z * rpb_y * fz_0 + 12.0 * fe_0 * fe_0 * rpb_z * rpb_y * fz_0 -
                               fke_0 * rpa_z * rpa_z * rpa_z * rpb_y * fz_0);

        fints_yzz[i] += fss * 12.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * fz_0;

        fints_yzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y + 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y +
                               (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y +
                               (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yzz[i] += ftt * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y;

        fints_zzz[i] += fss * (-(9.0 / 2.0) * fbe_0 * fe_0 * rpa_z * rpb_z * fz_0 - (9.0 / 2.0) * fbe_0 * fe_0 * rpb_z * rpb_z * fz_0 -
                               (9.0 / 4.0) * fbe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fbe_0 * rpa_z * rpb_z * rpb_z * rpb_z * fz_0 -
                               (9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpb_z * fz_0);

        fints_zzz[i] += fss * (-(9.0 / 2.0) * fe_0 * fke_0 * rpa_z * rpa_z * fz_0 + 15.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z * fz_0 +
                               45.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * fz_0 + 15.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * fz_0 -
                               (9.0 / 4.0) * fe_0 * fe_0 * fke_0 * fz_0);

        fints_zzz[i] +=
            fss * (54.0 * fe_0 * fe_0 * rpa_z * rpb_z * fz_0 + 18.0 * fe_0 * fe_0 * rpa_z * rpa_z * fz_0 + 18.0 * fe_0 * fe_0 * rpb_z * rpb_z * fz_0 +
                   (45.0 / 4.0) * fe_0 * fe_0 * fe_0 * fz_0 - 3.0 * fke_0 * rpa_z * rpa_z * rpa_z * rpb_z * fz_0);

        fints_zzz[i] += fss * 12.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * fz_0;

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z +
                               (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z + (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z +
                               (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z);

        fints_zzz[i] +=
            ftt * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z);
    }
}

}  // namespace kinrec
