#include "NuclearPotentialRecFF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "MathConst.hpp"
#include "BoysFunc.hpp"
#include "T2CDistributor.hpp"

namespace npotrec { // npotrec namespace

auto
compNuclearPotentialFF(      CSubMatrix* matrix,
                       const double charge,
                       const TPoint3D& point,
                       const CGtoBlock&  gto_block,
                       const int64_t     bra_first,
                       const int64_t     bra_last) -> void
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

        simd::loadCoordinates(ket_coords_x,
                              ket_coords_y,
                              ket_coords_z,
                              gto_coords,
                              ket_first,
                              ket_last);

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

                    npotrec::compPrimitiveNuclearPotentialFF_XXX_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, gto_indexes,
                                4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_5 * f3_3, gto_indexes,
                                6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, gto_indexes,
                                4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, f3_5 * f3_5, gto_indexes,
                                6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, gto_indexes,
                                4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_5 * 3.0 * f3_5, gto_indexes,
                                6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, gto_indexes,
                                4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_5 * f3_3, gto_indexes,
                                6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, gto_indexes,
                                4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_5 * 3.0, gto_indexes,
                                6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, gto_indexes,
                                4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_5 * 0.5 * f3_15, gto_indexes,
                                6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, gto_indexes,
                                4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_5 * f3_3, gto_indexes,
                                6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, gto_indexes,
                                4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_5 * 3.0 * f3_5, gto_indexes,
                                6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, gto_indexes,
                                4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, f3_5 * f3_15, gto_indexes,
                                6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, gto_indexes,
                                4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, f3_5 * 4.0 * f3_3, gto_indexes,
                                6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, gto_indexes,
                                4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_5 * f3_5, gto_indexes,
                                6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, gto_indexes,
                                4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_5 * f3_3, gto_indexes,
                                6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, gto_indexes,
                                4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_5 * 3.0, gto_indexes,
                                6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, gto_indexes,
                                4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_5 * 0.5 * f3_15, gto_indexes,
                                6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, gto_indexes,
                                4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, f3_5 * 4.0 * f3_3, gto_indexes,
                                6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, gto_indexes,
                                4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, f3_5 * 2.0, gto_indexes,
                                6, 3, j, ket_first, ket_last);

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

                    npotrec::compPrimitiveNuclearPotentialFF_XXY_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5 * f3_3, gto_indexes,
                                0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, gto_indexes,
                                2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5 * f3_5, gto_indexes,
                                0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, gto_indexes,
                                2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_5 * 3.0 * f3_5, gto_indexes,
                                0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, gto_indexes,
                                2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f3_5 * f3_3, gto_indexes,
                                0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, gto_indexes,
                                2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * f3_5 * 3.0, gto_indexes,
                                0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, gto_indexes,
                                2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f3_5 * 0.5 * f3_15, gto_indexes,
                                0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, gto_indexes,
                                2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_5 * f3_3, gto_indexes,
                                0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, gto_indexes,
                                2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_5 * 3.0 * f3_5, gto_indexes,
                                0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, gto_indexes,
                                2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f3_5 * f3_15, gto_indexes,
                                0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, gto_indexes,
                                2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, 3.0 * f3_5 * 4.0 * f3_3, gto_indexes,
                                0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, gto_indexes,
                                2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5 * f3_5, gto_indexes,
                                0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, gto_indexes,
                                2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5 * f3_3, gto_indexes,
                                0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, gto_indexes,
                                2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * f3_5 * 3.0, gto_indexes,
                                0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, gto_indexes,
                                2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * f3_5 * 0.5 * f3_15, gto_indexes,
                                0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, gto_indexes,
                                2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, 3.0 * f3_5 * 4.0 * f3_3, gto_indexes,
                                0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, gto_indexes,
                                2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * f3_5 * 2.0, gto_indexes,
                                0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, gto_indexes,
                                2, 3, j, ket_first, ket_last);

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

                    npotrec::compPrimitiveNuclearPotentialFF_XXZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_3, gto_indexes,
                                3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f3_15 * f3_3, gto_indexes,
                                5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5, gto_indexes,
                                3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f3_15 * f3_5, gto_indexes,
                                5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * 3.0 * f3_5, gto_indexes,
                                3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 0.5 * f3_15 * 3.0 * f3_5, gto_indexes,
                                5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_3, gto_indexes,
                                3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -0.5 * f3_15 * f3_3, gto_indexes,
                                5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 3.0, gto_indexes,
                                3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -0.5 * f3_15 * 3.0, gto_indexes,
                                5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 0.5 * f3_15, gto_indexes,
                                3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 0.5 * f3_15 * 0.5 * f3_15, gto_indexes,
                                5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_3, gto_indexes,
                                3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -0.5 * f3_15 * f3_3, gto_indexes,
                                5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * 3.0 * f3_5, gto_indexes,
                                3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -0.5 * f3_15 * 3.0 * f3_5, gto_indexes,
                                5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_15, gto_indexes,
                                3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, 0.5 * f3_15 * f3_15, gto_indexes,
                                5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -3.0 * 4.0 * f3_3, gto_indexes,
                                3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, 0.5 * f3_15 * 4.0 * f3_3, gto_indexes,
                                5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5, gto_indexes,
                                3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f3_15 * f3_5, gto_indexes,
                                5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_3, gto_indexes,
                                3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f3_15 * f3_3, gto_indexes,
                                5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 3.0, gto_indexes,
                                3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -0.5 * f3_15 * 3.0, gto_indexes,
                                5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 0.5 * f3_15, gto_indexes,
                                3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -0.5 * f3_15 * 0.5 * f3_15, gto_indexes,
                                5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -3.0 * 4.0 * f3_3, gto_indexes,
                                3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, 0.5 * f3_15 * 4.0 * f3_3, gto_indexes,
                                5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * 2.0, gto_indexes,
                                3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, 0.5 * f3_15 * 2.0, gto_indexes,
                                5, 3, j, ket_first, ket_last);

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

                    npotrec::compPrimitiveNuclearPotentialFF_XYY_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, gto_indexes,
                                4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5 * f3_3, gto_indexes,
                                6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, gto_indexes,
                                4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5 * f3_5, gto_indexes,
                                6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, gto_indexes,
                                4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f3_5 * 3.0 * f3_5, gto_indexes,
                                6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, gto_indexes,
                                4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_5 * f3_3, gto_indexes,
                                6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, gto_indexes,
                                4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f3_5 * 3.0, gto_indexes,
                                6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, gto_indexes,
                                4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * f3_5 * 0.5 * f3_15, gto_indexes,
                                6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, gto_indexes,
                                4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_5 * f3_3, gto_indexes,
                                6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, gto_indexes,
                                4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_5 * 3.0 * f3_5, gto_indexes,
                                6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, gto_indexes,
                                4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_5 * f3_15, gto_indexes,
                                6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, gto_indexes,
                                4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -3.0 * f3_5 * 4.0 * f3_3, gto_indexes,
                                6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, gto_indexes,
                                4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5 * f3_5, gto_indexes,
                                6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, gto_indexes,
                                4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5 * f3_3, gto_indexes,
                                6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, gto_indexes,
                                4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f3_5 * 3.0, gto_indexes,
                                6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, gto_indexes,
                                4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f3_5 * 0.5 * f3_15, gto_indexes,
                                6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, gto_indexes,
                                4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -3.0 * f3_5 * 4.0 * f3_3, gto_indexes,
                                6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, gto_indexes,
                                4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f3_5 * 2.0, gto_indexes,
                                6, 3, j, ket_first, ket_last);

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

                    npotrec::compPrimitiveNuclearPotentialFF_XYZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, -f3_15 * f3_3, gto_indexes,
                                1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, f3_15 * f3_5, gto_indexes,
                                1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_15 * 3.0 * f3_5, gto_indexes,
                                1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_15 * f3_3, gto_indexes,
                                1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_15 * 3.0, gto_indexes,
                                1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_15 * 0.5 * f3_15, gto_indexes,
                                1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_15 * f3_3, gto_indexes,
                                1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_15 * 3.0 * f3_5, gto_indexes,
                                1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, f3_15 * f3_15, gto_indexes,
                                1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, f3_15 * 4.0 * f3_3, gto_indexes,
                                1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_15 * f3_5, gto_indexes,
                                1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_15 * f3_3, gto_indexes,
                                1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_15 * 3.0, gto_indexes,
                                1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_15 * 0.5 * f3_15, gto_indexes,
                                1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, f3_15 * 4.0 * f3_3, gto_indexes,
                                1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, f3_15 * 2.0, gto_indexes,
                                1, 3, j, ket_first, ket_last);

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

                    npotrec::compPrimitiveNuclearPotentialFF_XZZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f3_3 * f3_3, gto_indexes,
                                4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f3_3 * f3_5, gto_indexes,
                                4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 4.0 * f3_3 * 3.0 * f3_5, gto_indexes,
                                4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f3_3 * f3_3, gto_indexes,
                                4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f3_3 * 3.0, gto_indexes,
                                4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 4.0 * f3_3 * 0.5 * f3_15, gto_indexes,
                                4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * f3_3, gto_indexes,
                                4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * 3.0 * f3_5, gto_indexes,
                                4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f3_3 * f3_15, gto_indexes,
                                4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes,
                                4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_5, gto_indexes,
                                4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_3, gto_indexes,
                                4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 3.0, gto_indexes,
                                4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 0.5 * f3_15, gto_indexes,
                                4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes,
                                4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, 4.0 * f3_3 * 2.0, gto_indexes,
                                4, 3, j, ket_first, ket_last);

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

                    npotrec::compPrimitiveNuclearPotentialFF_YYY_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_5 * f3_3, gto_indexes,
                                0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, gto_indexes,
                                2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_5 * f3_5, gto_indexes,
                                0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, gto_indexes,
                                2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_5 * 3.0 * f3_5, gto_indexes,
                                0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, gto_indexes,
                                2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_5 * f3_3, gto_indexes,
                                0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, gto_indexes,
                                2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_5 * 3.0, gto_indexes,
                                0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, gto_indexes,
                                2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_5 * 0.5 * f3_15, gto_indexes,
                                0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, gto_indexes,
                                2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_5 * f3_3, gto_indexes,
                                0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, gto_indexes,
                                2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_5 * 3.0 * f3_5, gto_indexes,
                                0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, gto_indexes,
                                2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_5 * f3_15, gto_indexes,
                                0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, gto_indexes,
                                2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_5 * 4.0 * f3_3, gto_indexes,
                                0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, gto_indexes,
                                2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_5 * f3_5, gto_indexes,
                                0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, gto_indexes,
                                2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_5 * f3_3, gto_indexes,
                                0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, gto_indexes,
                                2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_5 * 3.0, gto_indexes,
                                0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, gto_indexes,
                                2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_5 * 0.5 * f3_15, gto_indexes,
                                0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, gto_indexes,
                                2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_5 * 4.0 * f3_3, gto_indexes,
                                0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, gto_indexes,
                                2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_5 * 2.0, gto_indexes,
                                0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, gto_indexes,
                                2, 3, j, ket_first, ket_last);

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

                    npotrec::compPrimitiveNuclearPotentialFF_YYZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_3, gto_indexes,
                                3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f3_15 * f3_3, gto_indexes,
                                5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5, gto_indexes,
                                3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f3_15 * f3_5, gto_indexes,
                                5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * 3.0 * f3_5, gto_indexes,
                                3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -0.5 * f3_15 * 3.0 * f3_5, gto_indexes,
                                5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_3, gto_indexes,
                                3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 0.5 * f3_15 * f3_3, gto_indexes,
                                5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 3.0, gto_indexes,
                                3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 0.5 * f3_15 * 3.0, gto_indexes,
                                5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 0.5 * f3_15, gto_indexes,
                                3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -0.5 * f3_15 * 0.5 * f3_15, gto_indexes,
                                5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_3, gto_indexes,
                                3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 0.5 * f3_15 * f3_3, gto_indexes,
                                5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * 3.0 * f3_5, gto_indexes,
                                3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, 0.5 * f3_15 * 3.0 * f3_5, gto_indexes,
                                5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_15, gto_indexes,
                                3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, -0.5 * f3_15 * f3_15, gto_indexes,
                                5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -3.0 * 4.0 * f3_3, gto_indexes,
                                3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, -0.5 * f3_15 * 4.0 * f3_3, gto_indexes,
                                5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5, gto_indexes,
                                3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f3_15 * f3_5, gto_indexes,
                                5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_3, gto_indexes,
                                3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f3_15 * f3_3, gto_indexes,
                                5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 3.0, gto_indexes,
                                3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 0.5 * f3_15 * 3.0, gto_indexes,
                                5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 0.5 * f3_15, gto_indexes,
                                3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, 0.5 * f3_15 * 0.5 * f3_15, gto_indexes,
                                5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -3.0 * 4.0 * f3_3, gto_indexes,
                                3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, -0.5 * f3_15 * 4.0 * f3_3, gto_indexes,
                                5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * 2.0, gto_indexes,
                                3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, -0.5 * f3_15 * 2.0, gto_indexes,
                                5, 3, j, ket_first, ket_last);

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

                    npotrec::compPrimitiveNuclearPotentialFF_YZZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f3_3 * f3_3, gto_indexes,
                                2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f3_3 * f3_5, gto_indexes,
                                2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 4.0 * f3_3 * 3.0 * f3_5, gto_indexes,
                                2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f3_3 * f3_3, gto_indexes,
                                2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f3_3 * 3.0, gto_indexes,
                                2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 4.0 * f3_3 * 0.5 * f3_15, gto_indexes,
                                2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * f3_3, gto_indexes,
                                2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * 3.0 * f3_5, gto_indexes,
                                2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f3_3 * f3_15, gto_indexes,
                                2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes,
                                2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_5, gto_indexes,
                                2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_3, gto_indexes,
                                2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 3.0, gto_indexes,
                                2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 0.5 * f3_15, gto_indexes,
                                2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes,
                                2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, 4.0 * f3_3 * 2.0, gto_indexes,
                                2, 3, j, ket_first, ket_last);

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

                    npotrec::compPrimitiveNuclearPotentialFF_ZZZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, -2.0 * f3_3, gto_indexes,
                                3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxx, 2.0 * f3_5, gto_indexes,
                                3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, 2.0 * 3.0 * f3_5, gto_indexes,
                                3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxy, -2.0 * f3_3, gto_indexes,
                                3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, -2.0 * 3.0, gto_indexes,
                                3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xxz, 2.0 * 0.5 * f3_15, gto_indexes,
                                3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -2.0 * f3_3, gto_indexes,
                                3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyy, -2.0 * 3.0 * f3_5, gto_indexes,
                                3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xyz, 2.0 * f3_15, gto_indexes,
                                3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xzz, 2.0 * 4.0 * f3_3, gto_indexes,
                                3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -2.0 * f3_5, gto_indexes,
                                3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyy, -2.0 * f3_3, gto_indexes,
                                3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -2.0 * 3.0, gto_indexes,
                                3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yyz, -2.0 * 0.5 * f3_15, gto_indexes,
                                3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yzz, 2.0 * 4.0 * f3_3, gto_indexes,
                                3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zzz, 2.0 * 2.0, gto_indexes,
                                3, 3, j, ket_first, ket_last);

        }
    }
}

auto
compNuclearPotentialFF(      CSubMatrix* matrix,
                       const double charge,
                       const TPoint3D& point,
                       const CGtoBlock&  bra_gto_block,
                       const CGtoBlock&  ket_gto_block,
                       const int64_t     bra_first,
                       const int64_t     bra_last,
                       const mat_t       mat_type) -> void
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

                    npotrec::compPrimitiveNuclearPotentialFF_XXX_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, mat_type);

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

                    npotrec::compPrimitiveNuclearPotentialFF_XXY_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, mat_type);

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

                    npotrec::compPrimitiveNuclearPotentialFF_XXZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, mat_type);

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

                    npotrec::compPrimitiveNuclearPotentialFF_XYY_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                6, 3, j, ket_first, ket_last, mat_type);

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

                    npotrec::compPrimitiveNuclearPotentialFF_XYZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                1, 3, j, ket_first, ket_last, mat_type);

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

                    npotrec::compPrimitiveNuclearPotentialFF_XZZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                4, 3, j, ket_first, ket_last, mat_type);

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

                    npotrec::compPrimitiveNuclearPotentialFF_YYY_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, mat_type);

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

                    npotrec::compPrimitiveNuclearPotentialFF_YYZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                5, 3, j, ket_first, ket_last, mat_type);

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

                    npotrec::compPrimitiveNuclearPotentialFF_YZZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                2, 3, j, ket_first, ket_last, mat_type);

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

                    npotrec::compPrimitiveNuclearPotentialFF_ZZZ_T(buffer_xxx,
                                                                   buffer_xxy,
                                                                   buffer_xxz,
                                                                   buffer_xyy,
                                                                   buffer_xyz,
                                                                   buffer_xzz,
                                                                   buffer_yyy,
                                                                   buffer_yyz,
                                                                   buffer_yzz,
                                                                   buffer_zzz,
                                                                   charge,
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

            t2cfunc::distribute(matrix, buffer_xxx, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxx, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xxz, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyy, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xyz, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyy, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yyz, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yzz, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zzz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 3, j, ket_first, ket_last, mat_type);

        }
    }
}

auto
compPrimitiveNuclearPotentialFF_XXX_T(      TDoubleArray& buffer_xxx,
                                            TDoubleArray& buffer_xxy,
                                            TDoubleArray& buffer_xxz,
                                            TDoubleArray& buffer_xyy,
                                            TDoubleArray& buffer_xyz,
                                            TDoubleArray& buffer_xzz,
                                            TDoubleArray& buffer_yyy,
                                            TDoubleArray& buffer_yyz,
                                            TDoubleArray& buffer_yzz,
                                            TDoubleArray& buffer_zzz,
                       const double charge,
                       const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

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

    // set up Boys function variables

    const CBoysFunc<6> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<7> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<7>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_xxx,\
                             fints_xxy,\
                             fints_xxz,\
                             fints_xyy,\
                             fints_xyz,\
                             fints_xzz,\
                             fints_yyy,\
                             fints_yyz,\
                             fints_yzz,\
                             fints_zzz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxx[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x + (27.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_xxx[i] += fss * b0_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(27.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x - (27.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_x - (27.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x - (45.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(45.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x - (45.0 / 8.0) * fe_0 * fe_0 * fe_0 - 3.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_x - 3.0 * rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (27.0 * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x + (27.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_x + (27.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_x + 9.0 * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_x + (27.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (45.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_xxx[i] += fss * b2_vals[i] * ((45.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (45.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 + 3.0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (9.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-27.0 * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x - 15.0 * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x - 9.0 * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x - 15.0 * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x - 9.0 * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(45.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (45.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (45.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (15.0 / 8.0) * fe_0 * fe_0 * fe_0 - 9.0 * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-9.0 * rpa_x * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x - rpa_x * rpa_x * rpa_x * rpc_x * rpc_x * rpc_x - rpb_x * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (15.0 * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x + 15.0 * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x + (15.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x + (45.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + 9.0 * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (3.0 * rpa_x * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] * (-(15.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b6_vals[i] * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x + 3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xxy[i] += fss * b0_vals[i] * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x;

        fints_xxy[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y - 3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x);

        fints_xxy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y);

        fints_xxy[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - 3.0 * fe_0 * fe_0 * rpb_y * rpb_x - (15.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y - 3.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_x);

        fints_xxy[i] += fss * b1_vals[i] * (-2.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpc_x - rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_x + 9.0 * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x + 9.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y + 6.0 * fe_0 * rpb_y * rpb_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x + (15.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (3.0 * fe_0 * fe_0 * rpb_x * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 3.0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_x * rpc_x + 6.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (rpa_x * rpa_x * rpa_x * rpb_y * rpc_x * rpc_x + 2.0 * rpa_x * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x - 9.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x - 9.0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x - 6.0 * fe_0 * rpb_y * rpb_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-5.0 * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y - (15.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 6.0 * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x - 3.0 * rpa_x * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-6.0 * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x - rpa_x * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x - rpb_y * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (9.0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x + 5.0 * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + 5.0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (15.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (3.0 * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x + 6.0 * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x + 2.0 * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x + rpb_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpb_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b6_vals[i] * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x + 3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xxz[i] += fss * b0_vals[i] * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x;

        fints_xxz[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z - 3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x);

        fints_xxz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z);

        fints_xxz[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - 3.0 * fe_0 * fe_0 * rpb_z * rpb_x - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z - 3.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_x);

        fints_xxz[i] += fss * b1_vals[i] * (-2.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpc_x - rpa_x * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_x + 9.0 * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x + 9.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z + 6.0 * fe_0 * rpb_z * rpb_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (3.0 * fe_0 * fe_0 * rpb_x * rpc_z + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 3.0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_x * rpc_x + 6.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (rpa_x * rpa_x * rpa_x * rpb_z * rpc_x * rpc_x + 2.0 * rpa_x * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x - 9.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x - 9.0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x - 6.0 * fe_0 * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-5.0 * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x - 6.0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 6.0 * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x - 3.0 * rpa_x * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-6.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x - rpb_z * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (9.0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + 5.0 * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + 5.0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (3.0 * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x + 6.0 * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x + 2.0 * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - rpb_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b6_vals[i] * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_xyy[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x - (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_x - 3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y);

        fints_xyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 3.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_x - 2.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y - rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y + 9.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_x + 3.0 * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 3.0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_x * rpc_x + 6.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x + 3.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x + rpa_x * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x - 6.0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y - 3.0 * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x - 3.0 * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-6.0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x - 3.0 * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x - rpa_x * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x - rpb_y * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x + 6.0 * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (3.0 * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + rpb_y * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpb_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xyz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y + rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x - (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_y);

        fints_xyz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_xyz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - 3.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_x - rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpc_x - rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y - rpa_x * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y + (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z + (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpb_z * rpb_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 3.0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (3.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x + rpa_x * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x + rpa_x * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_z * rpb_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 3.0 * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-3.0 * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - 3.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x - 3.0 * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x - 3.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x - 3.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_x * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x - rpb_z * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 3.0 * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpb_z * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x + rpb_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + rpb_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpb_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_xzz[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x - (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x - 3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z);

        fints_xzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_x - 2.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z - rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z + 9.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_x + 3.0 * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + 3.0 * fe_0 * fe_0 * rpb_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 3.0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_x * rpc_x + 6.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x + 3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x + rpa_x * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x - 6.0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x - 3.0 * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - 3.0 * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x - 3.0 * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-6.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 3.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x - rpb_z * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x + 6.0 * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (3.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + rpb_z * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - 2.0 * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_yyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - 3.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_x);

        fints_yyy[i] += fss * b1_vals[i] * (-3.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y + (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 3.0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x + 9.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x + 3.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x);

        fints_yyy[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_yyy[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 9.0 * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x - 9.0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x - rpa_x * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 9.0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_yyy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 3.0 * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_yyy[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_yyz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_yyz[i] += fss * b1_vals[i] * (-3.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_x - 2.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y - rpa_x * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z + 3.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (3.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_x * rpc_x + 6.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x + 3.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x + rpa_x * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y + 2.0 * rpa_x * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 6.0 * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-3.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x - 3.0 * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x - 6.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y - rpb_z * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 3.0 * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + 6.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * (3.0 * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x + rpb_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_yyz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_yyz[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_yzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_yzz[i] += fss * b1_vals[i] * (-3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_x - 2.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z - rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z + 3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (3.0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_x * rpc_x + 6.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x + 3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y + rpa_x * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x - 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 6.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-3.0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - 6.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - 3.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 6.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + 3.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * (3.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x + rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_yzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - 2.0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_yzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_zzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - 3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x);

        fints_zzz[i] += fss * b1_vals[i] * (-3.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpa_x * rpc_z + (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 3.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x + 9.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x + 3.0 * rpa_x * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x);

        fints_zzz[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_zzz[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 9.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x - 9.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x - rpa_x * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 9.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + 3.0 * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x + 3.0 * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_zzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - 3.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - 3.0 * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_zzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialFF_XXY_T(      TDoubleArray& buffer_xxx,
                                            TDoubleArray& buffer_xxy,
                                            TDoubleArray& buffer_xxz,
                                            TDoubleArray& buffer_xyy,
                                            TDoubleArray& buffer_xyz,
                                            TDoubleArray& buffer_xzz,
                                            TDoubleArray& buffer_yyy,
                                            TDoubleArray& buffer_yyz,
                                            TDoubleArray& buffer_yzz,
                                            TDoubleArray& buffer_zzz,
                       const double charge,
                       const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

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

    // set up Boys function variables

    const CBoysFunc<6> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<7> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<7>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_xxx,\
                             fints_xxy,\
                             fints_xxz,\
                             fints_xyy,\
                             fints_xyz,\
                             fints_xzz,\
                             fints_yyy,\
                             fints_yyz,\
                             fints_yzz,\
                             fints_zzz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxx[i] += fss * b0_vals[i] * (3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xxx[i] += fss * b0_vals[i] * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x;

        fints_xxx[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x - 3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_y - 3.0 * fe_0 * fe_0 * rpa_y * rpa_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - 2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b1_vals[i] * (-3.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_x - rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y);

        fints_xxx[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x + 6.0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x + 9.0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpc_y);

        fints_xxx[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 6.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x + 3.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_x * rpc_x + rpa_y * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x + 3.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x - 9.0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x - 5.0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x - 9.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x - 6.0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x - 9.0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y);

        fints_xxx[i] += fss * b3_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - (15.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 6.0 * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-6.0 * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x - 3.0 * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x - rpb_x * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (5.0 * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x + 9.0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + 5.0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (15.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x + 6.0 * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + rpa_x * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b6_vals[i] * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxy[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xxy[i] += fss * b0_vals[i] * (fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxy[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x - 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y - 2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y - fe_0 * rpa_x * rpb_x * rpb_x * rpc_x);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y - fe_0 * rpa_x * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xxy[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - 2.0 * fe_0 * fe_0 * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_x - 2.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_x);

        fints_xxy[i] += fss * b1_vals[i] * (-rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y - rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x + 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y + 3.0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y + 2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y + 3.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y + 2.0 * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x + fe_0 * rpa_x * rpb_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y + fe_0 * rpa_x * rpa_x * rpb_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x + 3.0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y + fe_0 * fe_0 * rpa_x * rpb_x + 3.0 * fe_0 * fe_0 * rpa_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y + 3.0 * fe_0 * fe_0 * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 4.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x + rpa_y * rpa_x * rpa_x * rpb_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x + rpa_y * rpb_y * rpb_x * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x + rpa_x * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-2.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y - 2.0 * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - fe_0 * rpa_x * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x - fe_0 * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y - 3.0 * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-4.0 * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x - 2.0 * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x - rpa_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x - 4.0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-2.0 * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x - rpa_x * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x - rpb_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x + fe_0 * rpa_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (fe_0 * rpb_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x + rpa_y * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x + 2.0 * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x + 4.0 * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (rpa_x * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + rpb_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x - rpa_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpb_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-2.0 * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxz[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_x - 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z);

        fints_xxz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z - 2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - 2.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_x);

        fints_xxz[i] += fss * b1_vals[i] * (-2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_x - rpa_y * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z - rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x + 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z + 3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z + 3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z + 2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y + 3.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + 3.0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * (4.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x + 2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpb_z * rpc_x * rpc_x + 2.0 * rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x + rpa_y * rpb_z * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x + rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x - 4.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x - rpa_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-4.0 * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - rpb_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x + rpa_y * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x + 2.0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (4.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xyy[i] += fss * b0_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x + fe_0 * rpa_x * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x);

        fints_xyy[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + fe_0 * fe_0 * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x + rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_y - fe_0 * rpa_y * rpa_x * rpb_y * rpb_y - fe_0 * rpa_y * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x - 2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_x - fe_0 * rpa_x * rpb_y * rpb_y * rpc_y);

        fints_xyy[i] += fss * b1_vals[i] * (-fe_0 * rpa_x * rpa_x * rpb_y * rpb_x - fe_0 * rpa_x * rpa_x * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_y - fe_0 * fe_0 * rpa_y * rpa_x);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - 2.0 * fe_0 * fe_0 * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y - fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - 2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_x - 2.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y - rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_y + fe_0 * rpa_y * rpa_x * rpb_x * rpc_x + fe_0 * rpa_y * rpa_x * rpc_y * rpc_y + fe_0 * rpa_y * rpa_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_y * rpb_x * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_x + 2.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y + 2.0 * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x + fe_0 * rpa_x * rpb_y * rpb_y * rpc_y + 3.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpa_x * rpa_x * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + fe_0 * rpb_y * rpb_x * rpc_y * rpc_y + fe_0 * rpb_y * rpb_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * fe_0 * rpa_x * rpb_y + 3.0 * fe_0 * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x + 3.0 * fe_0 * fe_0 * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 4.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x + rpa_y * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (rpa_y * rpb_y * rpb_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x + 2.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y + rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpa_x * rpc_y * rpc_y - fe_0 * rpa_y * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x - 2.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y - 2.0 * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x - fe_0 * rpa_x * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x - fe_0 * rpb_y * rpb_x * rpc_y * rpc_y - fe_0 * rpb_y * rpb_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x - fe_0 * rpb_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-4.0 * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_y * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x - rpa_y * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x - rpa_y * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-4.0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x - rpa_x * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x + fe_0 * rpa_x * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (fe_0 * rpb_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x + rpa_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + 4.0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (2.0 * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + rpb_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - rpa_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xyz[i] += fss * b0_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xyz[i] += fss * b0_vals[i] * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x;

        fints_xyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpa_x * rpb_z * rpb_y - fe_0 * rpa_y * rpa_x * rpb_z * rpc_y - fe_0 * rpa_y * rpa_x * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z - fe_0 * rpa_x * rpb_z * rpb_y * rpc_y - fe_0 * rpa_x * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_y - fe_0 * fe_0 * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - 2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_x - rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y - rpa_y * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z - rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_z * rpc_y + fe_0 * rpa_y * rpa_x * rpb_y * rpc_z + fe_0 * rpa_y * rpa_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y + fe_0 * rpa_x * rpb_z * rpb_y * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_x * rpb_z * rpb_x * rpc_x + fe_0 * rpa_x * rpb_z * rpc_y * rpc_y + fe_0 * rpa_x * rpb_z * rpc_x * rpc_x + fe_0 * rpa_x * rpb_y * rpc_z * rpc_y + fe_0 * rpa_x * rpb_x * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z + fe_0 * fe_0 * rpa_x * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x + rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * (rpa_y * rpb_z * rpb_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x + rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x + rpa_x * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y + rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_x * rpb_z * rpc_y * rpc_y - fe_0 * rpa_x * rpb_z * rpc_x * rpc_x - fe_0 * rpa_x * rpb_y * rpc_z * rpc_y - fe_0 * rpa_x * rpb_x * rpc_z * rpc_x - fe_0 * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_x * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x - rpa_y * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_y * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - rpa_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x - rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y - rpb_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 2.0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + rpa_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x + rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (2.0 * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + rpb_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpb_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xzz[i] += fss * b0_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xzz[i] += fss * b0_vals[i] * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x;

        fints_xzz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z - fe_0 * rpa_y * rpa_x * rpb_z * rpb_z - fe_0 * rpa_y * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x - fe_0 * rpa_x * rpb_z * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_y - fe_0 * fe_0 * rpa_y * rpa_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - 2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_x - 2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z - rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x - rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z + fe_0 * rpa_y * rpa_x * rpb_x * rpc_x + fe_0 * rpa_y * rpa_x * rpc_z * rpc_z + fe_0 * rpa_y * rpa_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_z * rpb_x * rpc_z + 3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + fe_0 * rpa_x * rpb_z * rpb_z * rpc_y + fe_0 * rpa_x * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x + fe_0 * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 4.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x + 2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z + rpa_y * rpb_z * rpb_z * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y + rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpa_x * rpc_z * rpc_z - fe_0 * rpa_y * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - fe_0 * rpa_x * rpb_x * rpc_y * rpc_x - fe_0 * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x - fe_0 * rpb_z * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 4.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-rpa_y * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x - 4.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - 2.0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-rpb_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x + fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + fe_0 * rpa_x * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + rpa_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + 4.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpb_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - 2.0 * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_yyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_yyy[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_x - 3.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y - rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y + 9.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpb_y * rpb_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y + 3.0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (3.0 * fe_0 * fe_0 * rpa_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + (9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 6.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x + 3.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y + rpa_y * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x + 3.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x - 6.0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * fe_0 * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x - rpa_y * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-3.0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x - 6.0 * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y + 6.0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + 2.0 * rpa_y * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * (3.0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + 6.0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyy[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y - rpa_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyy[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yyz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y + fe_0 * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yyz[i] += fss * b0_vals[i] * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y;

        fints_yyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpa_x * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z - fe_0 * rpa_y * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z - 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x - fe_0 * rpa_x * rpa_x * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y - fe_0 * rpa_x * rpa_x * rpb_y * rpc_z);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - fe_0 * fe_0 * rpb_z * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z - 2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_x - 2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y - rpa_y * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z - rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_z * rpc_x + fe_0 * rpa_y * rpa_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z + fe_0 * rpa_y * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x + fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z + 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x + 3.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y + fe_0 * rpa_x * rpa_x * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y + fe_0 * fe_0 * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 4.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x + rpa_y * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y + 2.0 * rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y + rpa_y * rpb_z * rpb_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y + rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpa_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x - fe_0 * rpa_y * rpb_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x - 2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-fe_0 * rpb_z * rpb_y * rpc_y * rpc_y - fe_0 * rpb_z * rpb_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y - fe_0 * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-fe_0 * rpb_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 2.0 * rpa_y * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x - 4.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-rpa_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x - 4.0 * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y - 2.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-rpb_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + fe_0 * rpb_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + rpa_y * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + 4.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * (rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpb_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyz[i] += fss * b5_vals[i] * (-2.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyz[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yzz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_yzz[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpa_x * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y - fe_0 * rpa_y * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y - fe_0 * rpa_x * rpb_z * rpb_z * rpc_x - fe_0 * rpa_x * rpa_x * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_x);

        fints_yzz[i] += fss * b1_vals[i] * (-2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z - rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y - rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_y * rpc_x + fe_0 * rpa_y * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_y + fe_0 * rpa_y * rpb_z * rpb_y * rpc_z + fe_0 * rpa_y * rpb_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x + 2.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x + fe_0 * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * rpa_x * rpb_y * rpc_y * rpc_x + fe_0 * rpa_x * rpa_x * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y + fe_0 * fe_0 * rpa_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + fe_0 * fe_0 * rpb_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_yzz[i] += fss * b2_vals[i] * (4.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x + 2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y + rpa_y * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z + rpa_y * rpb_z * rpb_z * rpb_y * rpc_x * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x + 2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y + rpa_x * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpa_x * rpc_y * rpc_x - fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x - fe_0 * rpa_x * rpb_y * rpc_y * rpc_x - fe_0 * rpa_x * rpc_z * rpc_z * rpc_x - fe_0 * rpa_x * rpc_y * rpc_y * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y - fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - fe_0 * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (1.0 / 8.0) * fe_0 * fe_0 * fe_0 - 4.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x - rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x - rpa_y * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - 4.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x + fe_0 * rpa_x * rpc_z * rpc_z * rpc_x + fe_0 * rpa_x * rpc_y * rpc_y * rpc_x + fe_0 * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (fe_0 * rpb_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + 2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + 4.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + rpb_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yzz[i] += fss * b5_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yzz[i] += fss * b5_vals[i] * (-2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_zzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_zzz[i] += fss * b1_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x - 3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z - rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z + 3.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] * (6.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x + 3.0 * rpa_y * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z + rpa_y * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x + 3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_zzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 6.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_zzz[i] += fss * b3_vals[i] * (-rpa_y * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z - 3.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x - 6.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - 3.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z + 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x + 3.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * (6.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_zzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - 3.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_zzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialFF_XXZ_T(      TDoubleArray& buffer_xxx,
                                            TDoubleArray& buffer_xxy,
                                            TDoubleArray& buffer_xxz,
                                            TDoubleArray& buffer_xyy,
                                            TDoubleArray& buffer_xyz,
                                            TDoubleArray& buffer_xzz,
                                            TDoubleArray& buffer_yyy,
                                            TDoubleArray& buffer_yyz,
                                            TDoubleArray& buffer_yzz,
                                            TDoubleArray& buffer_zzz,
                       const double charge,
                       const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

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

    // set up Boys function variables

    const CBoysFunc<6> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<7> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<7>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_xxx,\
                             fints_xxy,\
                             fints_xxz,\
                             fints_xyy,\
                             fints_xyz,\
                             fints_xzz,\
                             fints_yyy,\
                             fints_yyz,\
                             fints_yzz,\
                             fints_zzz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxx[i] += fss * b0_vals[i] * (3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xxx[i] += fss * b0_vals[i] * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpb_x;

        fints_xxx[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x - 3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z - 3.0 * fe_0 * fe_0 * rpa_z * rpa_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - 2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b1_vals[i] * (-3.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_x - rpa_x * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxx[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x + 6.0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x + 9.0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x + 3.0 * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x + 3.0 * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_xxx[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 6.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_x * rpc_x + 3.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_x * rpc_x + rpa_z * rpb_x * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x + 3.0 * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x - 9.0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x - 5.0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x - 9.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x - 6.0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x - 9.0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_xxx[i] += fss * b3_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 6.0 * rpa_z * rpa_x * rpb_x * rpc_x * rpc_x * rpc_x - rpa_z * rpa_x * rpa_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_z * rpb_x * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-6.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x - 3.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x - rpb_x * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (5.0 * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x + 6.0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + 9.0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + 5.0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_x * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * rpa_z * rpb_x * rpc_x * rpc_x * rpc_x * rpc_x + 6.0 * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + rpa_x * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - rpa_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b6_vals[i] * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxy[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxy[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpb_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x - 2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y - 2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpc_x);

        fints_xxy[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_x - rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_y - rpa_x * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y + 3.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y + 3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y + 2.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z + 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + 3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (4.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x + rpa_z * rpa_x * rpa_x * rpb_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_x + rpa_z * rpb_y * rpb_x * rpb_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y - 3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_y * rpc_x * rpc_x * rpc_x - 4.0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x - rpa_z * rpa_x * rpa_x * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpb_y * rpb_x * rpc_x * rpc_x * rpc_x - rpa_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-4.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - rpb_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x + rpa_z * rpb_y * rpc_x * rpc_x * rpc_x * rpc_x + 2.0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (4.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - rpa_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpb_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxz[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_xxz[i] += fss * b0_vals[i] * (fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x - 2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z);

        fints_xxz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z - 2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z - fe_0 * rpa_x * rpb_x * rpb_x * rpc_x);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z - fe_0 * rpa_x * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_xxz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - 2.0 * fe_0 * fe_0 * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_xxz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpc_x - 2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_x);

        fints_xxz[i] += fss * b1_vals[i] * (-rpa_z * rpa_x * rpa_x * rpb_x * rpb_x * rpc_z - rpa_x * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x + 2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z + 3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z + 3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z + 2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z + 3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z + 2.0 * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x + fe_0 * rpa_x * rpb_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z + fe_0 * rpa_x * rpa_x * rpb_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x + 3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z + fe_0 * fe_0 * rpa_x * rpb_x + 3.0 * fe_0 * fe_0 * rpa_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 4.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + rpa_z * rpa_x * rpa_x * rpb_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_x + rpa_z * rpb_z * rpb_x * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-2.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z - 2.0 * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x - fe_0 * rpa_x * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x - fe_0 * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - 3.0 * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-4.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_z * rpb_z * rpb_x * rpc_x * rpc_x * rpc_x - rpa_z * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x - 4.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-2.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x - rpa_x * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x - rpb_z * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x + fe_0 * rpa_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (fe_0 * rpb_x * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x + rpa_z * rpb_z * rpc_x * rpc_x * rpc_x * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + 4.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (rpa_x * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x - rpa_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - rpb_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-2.0 * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xyy[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xyy[i] += fss * b0_vals[i] * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_x;

        fints_xyy[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y - fe_0 * rpa_z * rpa_x * rpb_y * rpb_y - fe_0 * rpa_z * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x - fe_0 * rpa_x * rpb_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z - fe_0 * fe_0 * rpa_z * rpa_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - 2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpc_x - 2.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_y - rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_x - rpa_x * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y + fe_0 * rpa_z * rpa_x * rpb_x * rpc_x + fe_0 * rpa_z * rpa_x * rpc_y * rpc_y + fe_0 * rpa_z * rpa_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpb_y * rpb_x * rpc_y + 3.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y + fe_0 * rpa_x * rpb_y * rpb_y * rpc_z + fe_0 * rpa_x * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpb_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x + fe_0 * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 4.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_x + rpa_z * rpa_x * rpa_x * rpb_x * rpc_y * rpc_y + rpa_z * rpb_y * rpb_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y + rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_x * rpc_y * rpc_y - fe_0 * rpa_z * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x - 2.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y - fe_0 * rpa_x * rpb_x * rpc_z * rpc_x - fe_0 * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-fe_0 * rpa_x * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x - fe_0 * rpb_y * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 4.0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-rpa_z * rpb_y * rpb_y * rpc_x * rpc_x * rpc_x - 4.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-rpb_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x + fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x + rpa_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + 4.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpb_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - rpa_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xyz[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xyz[i] += fss * b0_vals[i] * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_x;

        fints_xyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_x * rpb_z * rpb_y - fe_0 * rpa_z * rpa_x * rpb_z * rpc_y - fe_0 * rpa_z * rpa_x * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z - fe_0 * rpa_x * rpb_z * rpb_y * rpc_z - fe_0 * rpa_x * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z - fe_0 * fe_0 * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpc_x - rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_y - rpa_z * rpa_x * rpa_x * rpb_y * rpb_x * rpc_z - rpa_x * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_z * rpc_y + fe_0 * rpa_z * rpa_x * rpb_y * rpc_z + fe_0 * rpa_z * rpa_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y + fe_0 * rpa_x * rpb_z * rpb_y * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + fe_0 * rpa_x * rpb_y * rpb_x * rpc_x + fe_0 * rpa_x * rpb_y * rpc_z * rpc_z + fe_0 * rpa_x * rpb_y * rpc_x * rpc_x + fe_0 * rpa_x * rpb_x * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y + fe_0 * fe_0 * rpa_x * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x + rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_x + rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_x + rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * (rpa_z * rpb_z * rpb_y * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x + rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x + rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y + rpa_x * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - fe_0 * rpa_x * rpb_y * rpc_z * rpc_z - fe_0 * rpa_x * rpb_y * rpc_x * rpc_x - fe_0 * rpa_x * rpb_x * rpc_y * rpc_x - fe_0 * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_x - rpa_z * rpb_z * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - rpa_z * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x - rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y - rpb_z * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + fe_0 * rpa_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + rpa_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x + rpa_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (2.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + rpb_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpb_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xzz[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x + fe_0 * rpa_x * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x);

        fints_xzz[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + fe_0 * fe_0 * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x + rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z - fe_0 * rpa_z * rpa_x * rpb_z * rpb_z - fe_0 * rpa_z * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_x - fe_0 * rpa_x * rpb_z * rpb_z * rpc_z);

        fints_xzz[i] += fss * b1_vals[i] * (-fe_0 * rpa_x * rpa_x * rpb_z * rpb_x - fe_0 * rpa_x * rpa_x * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z - fe_0 * fe_0 * rpa_z * rpa_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - 2.0 * fe_0 * fe_0 * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z - fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - 2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpc_x - 2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z - rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-rpa_x * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z + fe_0 * rpa_z * rpa_x * rpb_x * rpc_x + fe_0 * rpa_z * rpa_x * rpc_z * rpc_z + fe_0 * rpa_z * rpa_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpb_z * rpb_x * rpc_z + 3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_x + 2.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z + 2.0 * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x + fe_0 * rpa_x * rpb_z * rpb_z * rpc_z + 3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * rpa_x * rpa_x * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x + fe_0 * rpb_z * rpb_x * rpc_z * rpc_z + fe_0 * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * fe_0 * rpa_x * rpb_z + 3.0 * fe_0 * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x + 3.0 * fe_0 * fe_0 * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 4.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_x + rpa_z * rpa_x * rpa_x * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * (rpa_z * rpb_z * rpb_z * rpb_x * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x + 2.0 * rpa_x * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z + rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_x * rpc_z * rpc_z - fe_0 * rpa_z * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z - 2.0 * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x - fe_0 * rpa_x * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_x - fe_0 * rpb_z * rpb_x * rpc_z * rpc_z - fe_0 * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - fe_0 * rpb_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-4.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x - rpa_z * rpb_z * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-4.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x - rpa_x * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + fe_0 * rpa_x * rpc_z * rpc_z * rpc_z + 3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (fe_0 * rpb_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + rpa_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + 4.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (2.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x + 2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + rpb_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - rpa_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - 2.0 * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-rpb_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x;

        fints_yyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z);

        fints_yyy[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpc_x - 3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_y - rpa_x * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z);

        fints_yyy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y + 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (6.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x + 3.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_y * rpc_y + rpa_z * rpb_y * rpb_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x + 3.0 * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 6.0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b3_vals[i] * (-rpa_z * rpa_x * rpa_x * rpc_y * rpc_y * rpc_y - 3.0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x - 6.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - 3.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y - rpb_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * (6.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - rpa_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_yyy[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yyz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_yyz[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_x * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z - fe_0 * rpa_z * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z - fe_0 * rpa_x * rpb_y * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z - fe_0 * rpa_x * rpa_x * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpc_x);

        fints_yyz[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_y - rpa_z * rpa_x * rpa_x * rpb_y * rpb_y * rpc_z - rpa_x * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_z * rpc_x + fe_0 * rpa_z * rpa_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z + fe_0 * rpa_z * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x + fe_0 * rpa_z * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z + fe_0 * rpa_x * rpb_z * rpc_z * rpc_x + 2.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (fe_0 * rpa_x * rpb_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z + fe_0 * rpa_x * rpa_x * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z + fe_0 * fe_0 * rpa_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z + fe_0 * fe_0 * rpb_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_yyz[i] += fss * b2_vals[i] * (4.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x + rpa_z * rpa_x * rpa_x * rpb_z * rpc_y * rpc_y + 2.0 * rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_y + rpa_z * rpb_z * rpb_y * rpb_y * rpc_x * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x + 2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y + rpa_x * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x - fe_0 * rpa_z * rpb_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x - fe_0 * rpa_x * rpb_z * rpc_z * rpc_x - 2.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x - fe_0 * rpa_x * rpc_z * rpc_z * rpc_x - fe_0 * rpa_x * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_y * rpc_y - fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (1.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x - 4.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x - rpa_z * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x - 4.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x - rpa_x * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-2.0 * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y - rpb_z * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + fe_0 * rpa_x * rpc_z * rpc_z * rpc_x + fe_0 * rpa_x * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x + fe_0 * rpb_y * rpc_z * rpc_z * rpc_y + fe_0 * rpb_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * (rpa_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + 4.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + rpb_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yyz[i] += fss * b5_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b5_vals[i] * (-rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yzz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y + fe_0 * rpa_x * rpa_x * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yzz[i] += fss * b0_vals[i] * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_y;

        fints_yzz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_x * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y - fe_0 * rpa_z * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y - 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x - fe_0 * rpa_x * rpa_x * rpb_z * rpb_y - fe_0 * rpa_x * rpa_x * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - fe_0 * fe_0 * rpb_z * rpb_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpc_x - 2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z - rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_y - rpa_x * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_y * rpc_x + fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_y + fe_0 * rpa_z * rpb_z * rpb_y * rpc_z + fe_0 * rpa_z * rpb_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x + 2.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x + 2.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x + fe_0 * rpa_x * rpa_x * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y + fe_0 * rpb_z * rpb_y * rpc_z * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y + fe_0 * fe_0 * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 4.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_y + rpa_z * rpa_x * rpa_x * rpb_y * rpc_z * rpc_z + rpa_z * rpb_z * rpb_z * rpb_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z + rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_x * rpc_y * rpc_x - fe_0 * rpa_z * rpb_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - 2.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-fe_0 * rpb_z * rpb_y * rpc_z * rpc_z - fe_0 * rpb_z * rpb_y * rpc_x * rpc_x - fe_0 * rpb_z * rpc_z * rpc_z * rpc_y - fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 4.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-rpa_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - 4.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y - rpa_x * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-rpb_z * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x + fe_0 * rpb_z * rpc_z * rpc_z * rpc_y + fe_0 * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + 4.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * (rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - 2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_yzz[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_yzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_zzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x);

        fints_zzz[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpa_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpc_x - 3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z - rpa_x * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z + 9.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x + 3.0 * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x + (9.0 / 2.0) * fe_0 * rpa_x * rpa_x * rpb_z * rpc_z + 3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * fe_0 * rpa_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpa_x + (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 6.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x + 3.0 * rpa_z * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z + rpa_z * rpb_z * rpb_z * rpb_z * rpc_x * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * (2.0 * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x + 3.0 * rpa_x * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x - 6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_x * rpa_x * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * fe_0 * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x - rpa_z * rpa_x * rpa_x * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x - 6.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x - 3.0 * rpa_x * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z + 6.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * (3.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + 6.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + rpa_x * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z + 3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_zzz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z - rpa_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - 3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_zzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialFF_XYY_T(      TDoubleArray& buffer_xxx,
                                            TDoubleArray& buffer_xxy,
                                            TDoubleArray& buffer_xxz,
                                            TDoubleArray& buffer_xyy,
                                            TDoubleArray& buffer_xyz,
                                            TDoubleArray& buffer_xzz,
                                            TDoubleArray& buffer_yyy,
                                            TDoubleArray& buffer_yyz,
                                            TDoubleArray& buffer_yzz,
                                            TDoubleArray& buffer_zzz,
                       const double charge,
                       const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

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

    // set up Boys function variables

    const CBoysFunc<6> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<7> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<7>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_xxx,\
                             fints_xxy,\
                             fints_xxz,\
                             fints_xyy,\
                             fints_xyz,\
                             fints_xzz,\
                             fints_yyy,\
                             fints_yyz,\
                             fints_yzz,\
                             fints_zzz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxx[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xxx[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y - 3.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_xxx[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y - 3.0 * rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x - rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y + 3.0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x + 9.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_x + 3.0 * fe_0 * rpa_y * rpa_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_x + 3.0 * fe_0 * fe_0 * rpa_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 6.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x + 3.0 * rpa_y * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (3.0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x + rpa_x * rpb_x * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxx[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x - 9.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x - 6.0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpa_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y - 3.0 * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x - 6.0 * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-rpa_y * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_y * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x - rpb_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (6.0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x + 2.0 * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (6.0 * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + rpa_y * rpa_y * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpa_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxy[i] += fss * b0_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y + fe_0 * rpa_y * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x);

        fints_xxy[i] += fss * b0_vals[i] * (fe_0 * fe_0 * rpa_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x + rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x);

        fints_xxy[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpa_x * rpb_y * rpc_y - 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x - fe_0 * rpa_y * rpa_x * rpb_x * rpb_x - 2.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_y - fe_0 * rpa_y * rpb_x * rpb_x * rpc_x);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_y - fe_0 * rpa_y * rpa_y * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_x - fe_0 * rpa_y * rpa_y * rpb_x * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-fe_0 * rpa_x * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_x - fe_0 * fe_0 * rpa_y * rpa_x);

        fints_xxy[i] += fss * b1_vals[i] * (-2.0 * fe_0 * fe_0 * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xxy[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y - 2.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y - 2.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x - rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_y * rpc_y + 2.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_x + fe_0 * rpa_y * rpa_x * rpc_y * rpc_y + fe_0 * rpa_y * rpa_x * rpc_x * rpc_x + 2.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_x + 2.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y + 2.0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x + fe_0 * rpa_y * rpb_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_x + fe_0 * rpa_y * rpa_y * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x + fe_0 * rpa_x * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y + fe_0 * rpb_y * rpb_x * rpc_y * rpc_y + fe_0 * rpb_y * rpb_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x + fe_0 * fe_0 * rpa_y * rpb_x + 3.0 * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x + 3.0 * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 4.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x + rpa_y * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x + rpa_y * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x + rpa_x * rpb_y * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpa_x * rpc_y * rpc_y - fe_0 * rpa_y * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_x - 2.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y - 2.0 * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x - fe_0 * rpa_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y - fe_0 * rpb_y * rpb_x * rpc_y * rpc_y - fe_0 * rpb_y * rpb_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x - fe_0 * rpb_x * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x - 4.0 * rpa_y * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x - 4.0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x - 2.0 * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x - rpa_y * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-rpa_y * rpa_y * rpb_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_y * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x - rpa_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y - rpb_y * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x + fe_0 * rpa_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + fe_0 * rpb_x * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x + 4.0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + rpa_y * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (rpa_x * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + rpb_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpa_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - rpb_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-2.0 * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xxz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z + fe_0 * rpa_y * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xxz[i] += fss * b0_vals[i] * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x;

        fints_xxz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpa_x * rpb_z * rpc_y - 2.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_z - fe_0 * rpa_y * rpa_y * rpb_z * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_x - fe_0 * rpa_y * rpa_y * rpb_x * rpc_z - fe_0 * rpa_x * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - fe_0 * fe_0 * rpb_z * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z - 2.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y - 2.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x - rpa_y * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z - rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_z * rpc_y + fe_0 * rpa_y * rpa_x * rpc_z * rpc_y + 2.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x + 2.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_x + fe_0 * rpa_y * rpa_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x + fe_0 * rpa_x * rpb_z * rpb_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x + fe_0 * rpa_x * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z + fe_0 * rpb_z * rpb_x * rpc_y * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x + fe_0 * fe_0 * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 4.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x + rpa_y * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x + 2.0 * rpa_y * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x + rpa_x * rpb_z * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpa_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x - 2.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x - fe_0 * rpa_x * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-fe_0 * rpb_z * rpb_x * rpc_y * rpc_y - fe_0 * rpb_z * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x - fe_0 * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-fe_0 * rpb_x * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 2.0 * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x - 4.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - 4.0 * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - 2.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-rpa_y * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x - rpa_y * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x - 2.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x - rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-rpb_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (fe_0 * rpb_x * rpc_z * rpc_y * rpc_y + fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + 4.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x + rpa_x * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (2.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - 2.0 * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-2.0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xyy[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x + fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xyy[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x - 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y - 2.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_x - fe_0 * rpa_y * rpb_y * rpb_y * rpc_y);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_x - fe_0 * rpa_y * rpa_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_x - 2.0 * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y - 2.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y);

        fints_xyy[i] += fss * b1_vals[i] * (-rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x - rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x + 3.0 * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y + 3.0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x + 2.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_x + 2.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x + fe_0 * rpa_y * rpb_y * rpb_y * rpc_y + 3.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_x + fe_0 * rpa_y * rpa_y * rpb_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y + 3.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_x + 3.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y + 3.0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x + fe_0 * fe_0 * rpa_y * rpb_y + 3.0 * fe_0 * fe_0 * rpa_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 4.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x + 2.0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x + rpa_y * rpa_y * rpa_x * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x + rpa_y * rpa_y * rpb_y * rpb_y * rpc_x * rpc_x + rpa_x * rpb_y * rpb_y * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x - 2.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y - 2.0 * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - fe_0 * rpb_y * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - 3.0 * fe_0 * fe_0 * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 4.0 * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y - 4.0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x - rpa_y * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-rpa_y * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y - rpa_x * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x - rpb_y * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x + fe_0 * rpa_y * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x + fe_0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x + 4.0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + rpa_y * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (rpa_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + rpb_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y - 2.0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - rpa_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-rpb_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xyz[i] += fss * b0_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_xyz[i] += fss * b0_vals[i] * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x;

        fints_xyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpa_x * rpb_z * rpb_x - fe_0 * rpa_y * rpa_x * rpb_z * rpc_x - fe_0 * rpa_y * rpa_x * rpb_x * rpc_z - fe_0 * rpa_y * rpb_z * rpb_y * rpc_y - fe_0 * rpa_y * rpb_z * rpb_x * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_x - fe_0 * fe_0 * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - 2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y - rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * rpc_y - rpa_y * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z - rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_z * rpc_x + fe_0 * rpa_y * rpa_x * rpb_x * rpc_z + fe_0 * rpa_y * rpa_x * rpc_z * rpc_x + fe_0 * rpa_y * rpb_z * rpb_y * rpc_y + fe_0 * rpa_y * rpb_z * rpb_x * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_z * rpc_y * rpc_y + fe_0 * rpa_y * rpb_z * rpc_x * rpc_x + fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + fe_0 * rpa_y * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z + fe_0 * fe_0 * rpa_y * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x + rpa_y * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x + rpa_y * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (rpa_y * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpb_z * rpb_y * rpc_x * rpc_x + rpa_y * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x + rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x + rpa_x * rpb_z * rpb_y * rpb_x * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpa_x * rpc_z * rpc_x - fe_0 * rpa_y * rpb_z * rpc_y * rpc_y - fe_0 * rpa_y * rpb_z * rpc_x * rpc_x - fe_0 * rpa_y * rpb_y * rpc_z * rpc_y - fe_0 * rpa_y * rpb_x * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - rpa_y * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x - rpa_y * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x - rpa_y * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x - rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_x * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x - rpa_x * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y - rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - rpb_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpa_y * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x + rpa_x * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + rpb_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpb_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - 2.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpb_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xzz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xzz[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpa_x * rpb_x * rpc_y - fe_0 * rpa_y * rpb_z * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_x - fe_0 * rpa_y * rpa_y * rpb_z * rpc_z);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_x - fe_0 * rpa_x * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y);

        fints_xzz[i] += fss * b1_vals[i] * (-2.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z - rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x - rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_x * rpc_y + fe_0 * rpa_y * rpa_x * rpc_y * rpc_x + 2.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + fe_0 * rpa_y * rpb_z * rpb_z * rpc_y + fe_0 * rpa_y * rpb_x * rpc_y * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_x + fe_0 * rpa_y * rpa_y * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * rpa_x * rpb_z * rpb_x * rpc_z + fe_0 * rpa_x * rpb_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x + fe_0 * fe_0 * rpa_y * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x + fe_0 * fe_0 * rpb_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_xzz[i] += fss * b2_vals[i] * (4.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y + 2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x + rpa_y * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x + rpa_y * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x + rpa_x * rpb_z * rpb_z * rpb_x * rpc_y * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpa_x * rpc_y * rpc_x - 2.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - fe_0 * rpa_y * rpb_x * rpc_y * rpc_x - fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - fe_0 * rpa_y * rpc_y * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x * rpc_x - fe_0 * rpa_x * rpb_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - fe_0 * rpb_z * rpb_x * rpc_z * rpc_x - fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - fe_0 * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (1.0 / 8.0) * fe_0 * fe_0 * fe_0 - 4.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y - 4.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - rpa_y * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x - rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-rpa_x * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x - rpb_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + fe_0 * rpa_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x + fe_0 * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] * (fe_0 * rpb_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + 2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (4.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpa_y * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] * (2.0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpb_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_yyy[i] += fss * b0_vals[i] * (3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y);

        fints_yyy[i] += fss * b0_vals[i] * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y;

        fints_yyy[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_y - 3.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y - 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_x - 3.0 * fe_0 * fe_0 * rpa_y * rpa_x);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - 2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y);

        fints_yyy[i] += fss * b1_vals[i] * (-3.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y - rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_y + 6.0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x + 9.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y + (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x + 3.0 * fe_0 * fe_0 * rpa_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x + (15.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 6.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x + 3.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x + rpa_x * rpb_y * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpa_y * rpa_x * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_x - 6.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x - 9.0 * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-5.0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y - 9.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (15.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 6.0 * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y - 6.0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x - rpa_y * rpa_y * rpa_x * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-3.0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * (6.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x + 5.0 * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y + 9.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x + 5.0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x + (15.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + rpa_y * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - rpa_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_yyz[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y - 3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_y - 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z - 2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_x - 3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - 2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_yyz[i] += fss * b1_vals[i] * (-2.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y - rpa_y * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z - rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_y + 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z + 3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y + 2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x + 3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y + 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z + 3.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (4.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y + 2.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x + rpa_y * rpa_y * rpa_x * rpb_z * rpc_y * rpc_y + 2.0 * rpa_y * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x + rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x + rpa_x * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x - 2.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y - 3.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y - 4.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y - 4.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - rpa_y * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-rpa_y * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 2.0 * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + 4.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * (rpa_x * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_yzz[i] += fss * b0_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y);

        fints_yzz[i] += fss * b0_vals[i] * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y;

        fints_yzz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z - fe_0 * rpa_y * rpa_x * rpb_z * rpb_z - fe_0 * rpa_y * rpa_x * rpb_y * rpc_y - fe_0 * rpa_y * rpb_z * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_x - fe_0 * rpa_x * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_x - fe_0 * fe_0 * rpa_y * rpa_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - 2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z - rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpc_y - rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z + fe_0 * rpa_y * rpa_x * rpb_y * rpc_y + fe_0 * rpa_y * rpa_x * rpc_z * rpc_z + fe_0 * rpa_y * rpa_x * rpc_y * rpc_y + 2.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_z * rpb_z * rpc_x + fe_0 * rpa_y * rpb_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * rpa_x * rpb_z * rpb_y * rpc_z + 3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x + fe_0 * fe_0 * rpa_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 4.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y + 2.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x + 2.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y + rpa_y * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z + 2.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x + rpa_x * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpa_x * rpc_z * rpc_z - fe_0 * rpa_y * rpa_x * rpc_y * rpc_y - 2.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x - fe_0 * rpa_y * rpb_y * rpc_y * rpc_x - fe_0 * rpa_y * rpc_z * rpc_z * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y - fe_0 * rpb_z * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 4.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y - 4.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-rpa_y * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x - rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - rpa_x * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * (fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + fe_0 * rpa_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y + 4.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yzz[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_zzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_zzz[i] += fss * b1_vals[i] * (-2.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y - 3.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z - rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpc_y + 3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * (6.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x + 3.0 * rpa_y * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z + 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x + rpa_x * rpb_z * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_zzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 6.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b3_vals[i] * (-6.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - rpa_y * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z - 3.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x - 3.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y - rpb_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 2.0 * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + 6.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * (rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x + 3.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + 3.0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_zzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - 2.0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - 3.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_zzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialFF_XYZ_T(      TDoubleArray& buffer_xxx,
                                            TDoubleArray& buffer_xxy,
                                            TDoubleArray& buffer_xxz,
                                            TDoubleArray& buffer_xyy,
                                            TDoubleArray& buffer_xyz,
                                            TDoubleArray& buffer_xzz,
                                            TDoubleArray& buffer_yyy,
                                            TDoubleArray& buffer_yyz,
                                            TDoubleArray& buffer_yzz,
                                            TDoubleArray& buffer_zzz,
                       const double charge,
                       const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

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

    // set up Boys function variables

    const CBoysFunc<6> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<7> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<7>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_xxx,\
                             fints_xxy,\
                             fints_xxz,\
                             fints_xyy,\
                             fints_xyz,\
                             fints_xzz,\
                             fints_yyy,\
                             fints_yyz,\
                             fints_yzz,\
                             fints_zzz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxx[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y + rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y);

        fints_xxx[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xxx[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - 3.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_x - rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_x - rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpc_y - rpa_y * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 3.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_x * rpc_x + 3.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x * rpc_x + 3.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (rpa_z * rpb_x * rpb_x * rpb_x * rpc_y * rpc_x + 3.0 * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x + rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxx[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - rpa_z * rpa_y * rpa_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x * rpc_x - 3.0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x * rpc_x - 3.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x - 3.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-3.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + rpa_z * rpa_y * rpc_x * rpc_x * rpc_x * rpc_x + rpa_z * rpa_x * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x * rpc_x + rpa_y * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (3.0 * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - rpa_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - rpa_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxy[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y + fe_0 * rpa_z * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xxy[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpb_x;

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y - fe_0 * rpa_z * rpa_y * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x - fe_0 * rpa_z * rpa_y * rpb_x * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y - fe_0 * rpa_z * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpb_x - fe_0 * rpa_z * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z - fe_0 * rpa_y * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x - fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xxy[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z - 2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_x - rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_x - rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpc_y - rpa_y * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x + fe_0 * rpa_z * rpa_y * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x + fe_0 * rpa_z * rpb_y * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpb_x * rpc_y * rpc_y + fe_0 * rpa_z * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x + fe_0 * rpa_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y + fe_0 * rpa_x * rpb_x * rpc_z * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z + fe_0 * rpb_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xxy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z + fe_0 * fe_0 * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + rpa_z * rpa_y * rpa_x * rpb_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x * rpc_x + rpa_z * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_x + rpa_z * rpa_x * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (rpa_z * rpb_y * rpb_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x + rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y + rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x + rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x - fe_0 * rpa_z * rpb_x * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y - fe_0 * rpa_x * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x - fe_0 * rpb_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - fe_0 * rpb_x * rpc_z * rpc_y * rpc_y - fe_0 * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-rpa_z * rpa_y * rpa_x * rpc_y * rpc_x * rpc_x - rpa_z * rpa_y * rpb_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x - rpa_z * rpa_x * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x * rpc_x - rpa_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x - rpa_y * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - rpa_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y - rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + fe_0 * rpb_x * rpc_z * rpc_y * rpc_y + fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + rpa_z * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x + rpa_z * rpa_x * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpb_y * rpc_y * rpc_x * rpc_x * rpc_x + 2.0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (rpa_y * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x + 2.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_x * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (2.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - rpa_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xxz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z + fe_0 * rpa_z * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xxz[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpb_x;

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z - fe_0 * rpa_z * rpa_y * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x - fe_0 * rpa_z * rpa_y * rpb_x * rpc_z);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y - fe_0 * rpa_z * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z - fe_0 * rpa_y * rpa_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x - fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y - 2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_x - rpa_z * rpa_y * rpa_x * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b1_vals[i] * (-rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_x - rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpc_y - rpa_y * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x + fe_0 * rpa_z * rpa_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y + fe_0 * rpa_z * rpb_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x + fe_0 * rpa_z * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x + fe_0 * rpa_y * rpb_z * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_x * rpc_z * rpc_z + fe_0 * rpa_y * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + fe_0 * rpa_x * rpb_x * rpc_y * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y + fe_0 * rpb_z * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xxz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y + fe_0 * fe_0 * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + rpa_z * rpa_y * rpa_x * rpb_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x * rpc_x + rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_x + rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * (rpa_z * rpb_z * rpb_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x + rpa_y * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z + rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x + rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x - fe_0 * rpa_z * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x - fe_0 * rpa_y * rpb_x * rpc_z * rpc_z - fe_0 * rpa_y * rpb_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - fe_0 * rpa_x * rpb_x * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x - fe_0 * rpb_z * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - fe_0 * rpb_x * rpc_z * rpc_z * rpc_y - fe_0 * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-rpa_z * rpa_y * rpa_x * rpc_z * rpc_x * rpc_x - rpa_z * rpa_y * rpb_z * rpc_x * rpc_x * rpc_x - 2.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x - rpa_z * rpa_x * rpb_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x * rpc_x - rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_y * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y - rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x + fe_0 * rpb_x * rpc_z * rpc_z * rpc_y + fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + rpa_z * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x + rpa_z * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpb_z * rpc_y * rpc_x * rpc_x * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (rpa_y * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + rpa_y * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + 2.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + rpa_x * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (2.0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-rpb_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xyy[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y + fe_0 * rpa_z * rpa_x * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y);

        fints_xyy[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_x;

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x - fe_0 * rpa_z * rpa_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_x * rpb_y * rpb_x - fe_0 * rpa_z * rpa_x * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y - fe_0 * rpa_z * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z - fe_0 * rpa_x * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y - fe_0 * fe_0 * rpa_z * rpb_y);

        fints_xyy[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_y - rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_x - rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpc_y - rpa_y * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x + fe_0 * rpa_z * rpa_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + fe_0 * rpa_z * rpb_y * rpb_x * rpc_x + fe_0 * rpa_z * rpb_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpb_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x + fe_0 * rpa_x * rpb_y * rpb_x * rpc_z + fe_0 * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y + fe_0 * rpb_y * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y);

        fints_xyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z + fe_0 * fe_0 * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (rpa_z * rpa_y * rpa_x * rpb_x * rpc_y * rpc_y + 2.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x + rpa_z * rpa_y * rpb_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y * rpc_y + rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (rpa_z * rpb_y * rpb_y * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y + rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x + rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x + rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x - fe_0 * rpa_z * rpb_y * rpc_y * rpc_y - fe_0 * rpa_z * rpb_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x - fe_0 * rpa_y * rpb_y * rpc_z * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - fe_0 * rpa_x * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - fe_0 * rpb_y * rpb_x * rpc_z * rpc_x - fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - fe_0 * rpb_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x - rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_x - rpa_z * rpa_x * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x - rpa_z * rpb_y * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - rpa_y * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-rpa_y * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + fe_0 * rpb_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + rpa_z * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b4_vals[i] * (2.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpb_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - rpa_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-2.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xyz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xyz[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_xyz[i] += fss * b1_vals[i] * (-rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_x - rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_y - rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * rpc_z - rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_x - rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpc_y);

        fints_xyz[i] += fss * b1_vals[i] * (-rpa_y * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_x + rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_y + rpa_z * rpa_y * rpb_z * rpb_y * rpc_x * rpc_x + rpa_z * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x + rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x + rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (rpa_z * rpa_x * rpb_z * rpb_x * rpc_y * rpc_y + rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y + rpa_z * rpb_z * rpb_y * rpb_x * rpc_y * rpc_x + rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x + rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * (rpa_y * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z + rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x + rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (1.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_x - rpa_z * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x - rpa_z * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x - rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x - rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y - rpa_z * rpb_z * rpb_y * rpc_y * rpc_x * rpc_x - rpa_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x - rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x - rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y - rpa_y * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x - rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x - rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - rpa_x * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y - rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y - rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + rpa_z * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpa_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpa_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] * (rpb_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xzz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z + fe_0 * rpa_y * rpa_x * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_xzz[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_x;

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x - fe_0 * rpa_z * rpa_y * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y - fe_0 * rpa_y * rpa_x * rpb_z * rpb_x - fe_0 * rpa_y * rpa_x * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z);

        fints_xzz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z - fe_0 * rpa_x * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xzz[i] += fss * b1_vals[i] * (-fe_0 * fe_0 * rpa_y * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y - 2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z - rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_x - rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpc_y - rpa_y * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_x + fe_0 * rpa_z * rpa_y * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + fe_0 * rpa_z * rpb_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x + fe_0 * rpa_y * rpb_z * rpb_x * rpc_x + fe_0 * rpa_y * rpb_z * rpc_z * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x + fe_0 * rpa_x * rpb_z * rpb_x * rpc_y + fe_0 * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y + fe_0 * rpb_z * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z + fe_0 * fe_0 * rpb_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (rpa_z * rpa_y * rpa_x * rpb_x * rpc_z * rpc_z + 2.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x + rpa_z * rpa_y * rpb_z * rpb_z * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y + rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (rpa_z * rpb_z * rpb_z * rpb_x * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z + rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x + rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x + rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_x - fe_0 * rpa_z * rpb_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_x - fe_0 * rpa_y * rpb_z * rpc_z * rpc_z - fe_0 * rpa_y * rpb_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z - fe_0 * rpa_x * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - fe_0 * rpb_z * rpb_x * rpc_y * rpc_x - fe_0 * rpb_z * rpc_z * rpc_z * rpc_y - fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x - rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - rpa_z * rpb_z * rpb_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x - rpa_y * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z - 2.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-rpa_y * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (fe_0 * rpb_z * rpc_z * rpc_z * rpc_y + fe_0 * rpb_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + rpa_z * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x + rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x + 2.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x + 2.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] * (2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpb_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_yyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x + rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - 3.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_y - rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_x - rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpc_y - rpa_y * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y + 3.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 3.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_y * rpc_y + 3.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x + 3.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (rpa_z * rpb_y * rpb_y * rpb_y * rpc_y * rpc_x + 3.0 * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y + rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x + rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - rpa_z * rpa_y * rpa_x * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yyy[i] += fss * b3_vals[i] * (-3.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x + rpa_z * rpa_x * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] * (3.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + 3.0 * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpa_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_yyz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z + fe_0 * rpa_z * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y);

        fints_yyz[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpb_y;

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x - fe_0 * rpa_z * rpa_x * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y);

        fints_yyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_x * rpb_y * rpc_z - fe_0 * rpa_z * rpb_z * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z - fe_0 * rpa_y * rpa_x * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x - fe_0 * rpa_x * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_yyz[i] += fss * b1_vals[i] * (-fe_0 * fe_0 * rpa_x * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x - 2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_y - rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b1_vals[i] * (-rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_x - rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpc_y - rpa_y * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y + fe_0 * rpa_z * rpa_x * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y + fe_0 * rpa_z * rpb_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x + fe_0 * rpa_z * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x + fe_0 * rpa_y * rpb_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x + fe_0 * rpa_x * rpb_z * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + fe_0 * rpa_x * rpb_y * rpc_z * rpc_z + fe_0 * rpa_x * rpb_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y + fe_0 * rpb_z * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y + fe_0 * fe_0 * rpb_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + rpa_z * rpa_y * rpa_x * rpb_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x + rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y * rpc_y + rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (rpa_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y + rpa_y * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z + rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x + rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x - fe_0 * rpa_z * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x - fe_0 * rpa_y * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - fe_0 * rpa_x * rpb_y * rpc_z * rpc_z - fe_0 * rpa_x * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y - fe_0 * rpb_z * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - fe_0 * rpb_y * rpc_z * rpc_z * rpc_x - fe_0 * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-rpa_z * rpa_y * rpa_x * rpc_z * rpc_y * rpc_y - rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x - rpa_z * rpa_x * rpb_z * rpc_y * rpc_y * rpc_y - 2.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x - rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - rpa_y * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x + fe_0 * rpb_y * rpc_z * rpc_z * rpc_x + fe_0 * rpb_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y + rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * (rpa_y * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y + rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yyz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b5_vals[i] * (-rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_yzz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z + fe_0 * rpa_y * rpa_x * rpb_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z);

        fints_yzz[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_y;

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x - fe_0 * rpa_z * rpa_x * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x - fe_0 * rpa_y * rpa_x * rpb_z * rpb_y - fe_0 * rpa_y * rpa_x * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z);

        fints_yzz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpb_z * rpb_y * rpc_x - fe_0 * rpa_x * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_yzz[i] += fss * b1_vals[i] * (-fe_0 * fe_0 * rpa_x * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x - 2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z - rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_y);

        fints_yzz[i] += fss * b1_vals[i] * (-rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_x - rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpc_y - rpa_y * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x + fe_0 * rpa_z * rpa_x * rpb_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y + fe_0 * rpa_z * rpb_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_x * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y + fe_0 * rpa_y * rpb_z * rpb_y * rpc_x + fe_0 * rpa_y * rpb_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x + fe_0 * rpa_x * rpb_z * rpb_y * rpc_y + fe_0 * rpa_x * rpb_z * rpc_z * rpc_z + fe_0 * rpa_x * rpb_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y + fe_0 * rpb_z * rpb_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z + fe_0 * fe_0 * rpb_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (rpa_z * rpa_y * rpa_x * rpb_y * rpc_z * rpc_z + 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x + rpa_z * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y + rpa_z * rpa_x * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (rpa_z * rpb_z * rpb_z * rpb_y * rpc_y * rpc_x + 2.0 * rpa_y * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z + rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y + rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x + rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_y * rpc_y - fe_0 * rpa_z * rpb_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpc_z * rpc_y - fe_0 * rpa_y * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - fe_0 * rpa_x * rpb_z * rpc_z * rpc_z - fe_0 * rpa_x * rpb_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z - fe_0 * rpb_z * rpb_y * rpc_y * rpc_x - fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - fe_0 * rpb_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x - rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y - rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - rpa_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y - rpa_y * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z - 2.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - rpa_x * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y - rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b4_vals[i] * (fe_0 * rpb_z * rpc_z * rpc_z * rpc_x + fe_0 * rpb_z * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * (rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + 2.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x + 2.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b5_vals[i] * (-2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_zzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x + rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - 3.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z - rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_x - rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpc_y - rpa_y * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpa_x * rpb_z * rpc_z + 3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 3.0 * rpa_z * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z + 3.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x + 3.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] * (rpa_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_x + 3.0 * rpa_y * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z + rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x + rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_y * rpa_x * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_zzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - rpa_z * rpa_y * rpa_x * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x - 3.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y - 3.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - 3.0 * rpa_y * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z - 3.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x + rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_y * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b4_vals[i] * (3.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + 3.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_zzz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y - 3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_zzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialFF_XZZ_T(      TDoubleArray& buffer_xxx,
                                            TDoubleArray& buffer_xxy,
                                            TDoubleArray& buffer_xxz,
                                            TDoubleArray& buffer_xyy,
                                            TDoubleArray& buffer_xyz,
                                            TDoubleArray& buffer_xzz,
                                            TDoubleArray& buffer_yyy,
                                            TDoubleArray& buffer_yyz,
                                            TDoubleArray& buffer_yzz,
                                            TDoubleArray& buffer_zzz,
                       const double charge,
                       const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

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

    // set up Boys function variables

    const CBoysFunc<6> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<7> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<7>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_xxx,\
                             fints_xxy,\
                             fints_xxz,\
                             fints_xyy,\
                             fints_xyz,\
                             fints_xzz,\
                             fints_yyy,\
                             fints_yyz,\
                             fints_yzz,\
                             fints_zzz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxx[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xxx[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z - 3.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_xxx[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpb_x * rpc_z - 3.0 * rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * rpc_x - rpa_z * rpa_z * rpb_x * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z + 3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x + 9.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_x + 3.0 * fe_0 * rpa_z * rpa_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_x + 3.0 * fe_0 * fe_0 * rpa_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x + (9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 6.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_x + 3.0 * rpa_z * rpa_z * rpa_x * rpb_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (3.0 * rpa_z * rpa_z * rpb_x * rpb_x * rpc_x * rpc_x + rpa_x * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxx[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x - 9.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x - 6.0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpa_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - 3.0 * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_x * rpc_x - 6.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-rpa_z * rpa_z * rpa_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_z * rpa_z * rpb_x * rpc_x * rpc_x * rpc_x - 3.0 * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x - rpb_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (6.0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (6.0 * rpa_z * rpb_x * rpc_z * rpc_x * rpc_x * rpc_x + rpa_z * rpa_z * rpc_x * rpc_x * rpc_x * rpc_x + 3.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_x * rpc_x * rpc_x * rpc_x - 2.0 * rpa_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x * rpc_x;

        fints_xxy[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y + fe_0 * rpa_z * rpa_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xxy[i] += fss * b0_vals[i] * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x;

        fints_xxy[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_x * rpb_y * rpc_z - 2.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_y - fe_0 * rpa_z * rpa_z * rpb_y * rpb_x);

        fints_xxy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_x - fe_0 * rpa_z * rpa_z * rpb_x * rpc_y - fe_0 * rpa_x * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - fe_0 * fe_0 * rpb_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y - 2.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * rpc_x - rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * rpc_y - rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_y * rpc_z + fe_0 * rpa_z * rpa_x * rpc_z * rpc_y + 2.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z + 3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x + 2.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_x + fe_0 * rpa_z * rpa_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x + fe_0 * rpa_x * rpb_y * rpb_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x + fe_0 * rpa_x * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_y + fe_0 * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpb_y * rpb_x * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x + fe_0 * fe_0 * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 4.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_x + rpa_z * rpa_z * rpa_x * rpb_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_z * rpa_x * rpb_x * rpc_y * rpc_x + 2.0 * rpa_z * rpa_z * rpb_y * rpb_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (rpa_z * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x + rpa_x * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxy[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x - 2.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_x * rpc_x - fe_0 * rpa_x * rpb_x * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-fe_0 * rpb_y * rpb_x * rpc_z * rpc_z - fe_0 * rpb_y * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x - fe_0 * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-fe_0 * rpb_x * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 2.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_x * rpc_x - 4.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_x - 4.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-rpa_z * rpa_z * rpa_x * rpc_y * rpc_x * rpc_x - rpa_z * rpa_z * rpb_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_z * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x - 2.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x - rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-rpb_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (fe_0 * rpb_x * rpc_z * rpc_z * rpc_y + fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpb_y * rpc_z * rpc_x * rpc_x * rpc_x + 4.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_z * rpa_z * rpc_y * rpc_x * rpc_x * rpc_x + rpa_x * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (2.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpb_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-2.0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xxz[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z + fe_0 * rpa_z * rpa_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x);

        fints_xxz[i] += fss * b0_vals[i] * (fe_0 * fe_0 * rpa_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x + rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_x * rpb_z * rpc_z - 2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x - fe_0 * rpa_z * rpa_x * rpb_x * rpb_x - 2.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z - fe_0 * rpa_z * rpb_x * rpb_x * rpc_x);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_z - fe_0 * rpa_z * rpa_z * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_x - fe_0 * rpa_z * rpa_z * rpb_x * rpc_z);

        fints_xxz[i] += fss * b1_vals[i] * (-fe_0 * rpa_x * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_x - fe_0 * fe_0 * rpa_z * rpa_x);

        fints_xxz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * fe_0 * rpa_z * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z - 2.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * rpc_x - rpa_z * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b1_vals[i] * (-rpa_z * rpa_z * rpb_z * rpb_x * rpb_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_z * rpc_z + 2.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_x + fe_0 * rpa_z * rpa_x * rpc_z * rpc_z + fe_0 * rpa_z * rpa_x * rpc_x * rpc_x + 2.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x + 2.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z + 2.0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x + fe_0 * rpa_z * rpb_x * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_x + fe_0 * rpa_z * rpa_z * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x + fe_0 * rpa_x * rpb_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpb_x * rpc_z + fe_0 * rpb_z * rpb_x * rpc_z * rpc_z + fe_0 * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x + fe_0 * fe_0 * rpa_z * rpb_x + 3.0 * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x + 3.0 * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 4.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpb_x * rpb_x * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_x + rpa_z * rpa_z * rpa_x * rpb_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_z * rpa_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpa_z * rpb_z * rpb_x * rpc_x * rpc_x + rpa_z * rpa_z * rpb_x * rpb_x * rpc_z * rpc_x + rpa_x * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_x * rpc_z * rpc_z - fe_0 * rpa_z * rpa_x * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x - 2.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z - 2.0 * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x - fe_0 * rpa_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z - fe_0 * rpb_z * rpb_x * rpc_z * rpc_z - fe_0 * rpb_z * rpb_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x - fe_0 * rpb_x * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_x * rpc_x - 4.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_x - 4.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x - rpa_z * rpa_z * rpa_x * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-rpa_z * rpa_z * rpb_z * rpc_x * rpc_x * rpc_x - 2.0 * rpa_z * rpa_z * rpb_x * rpc_z * rpc_x * rpc_x - 2.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x - rpa_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z - rpb_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x + fe_0 * rpa_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + fe_0 * rpb_x * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_z * rpb_z * rpc_z * rpc_x * rpc_x * rpc_x + 4.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + rpa_z * rpa_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (rpa_x * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x + 2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - rpb_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-2.0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x;

        fints_xyy[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xyy[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_x * rpb_x * rpc_z - fe_0 * rpa_z * rpb_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_x - fe_0 * rpa_z * rpa_z * rpb_y * rpc_y);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_x - fe_0 * rpa_x * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyy[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * rpc_y - rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpc_x - rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_x * rpc_z + fe_0 * rpa_z * rpa_x * rpc_z * rpc_x + 2.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y + fe_0 * rpa_z * rpb_y * rpb_y * rpc_z + fe_0 * rpa_z * rpb_x * rpc_z * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_x + fe_0 * rpa_z * rpa_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpa_x * rpb_y * rpb_x * rpc_y + fe_0 * rpa_x * rpb_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpb_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x + fe_0 * fe_0 * rpa_z * rpc_z);

        fints_xyy[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x + fe_0 * fe_0 * rpb_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y);

        fints_xyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_xyy[i] += fss * b2_vals[i] * (4.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_x + 2.0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpc_y * rpc_x + rpa_z * rpa_z * rpa_x * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x + rpa_z * rpa_z * rpb_y * rpb_y * rpc_x * rpc_x + rpa_x * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyy[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_x * rpc_z * rpc_x - 2.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y - fe_0 * rpa_z * rpb_x * rpc_z * rpc_x - fe_0 * rpa_z * rpc_z * rpc_y * rpc_y - fe_0 * rpa_z * rpc_z * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x * rpc_x - fe_0 * rpa_x * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x - fe_0 * rpb_y * rpb_x * rpc_y * rpc_x - fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - fe_0 * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_xyy[i] += fss * b3_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (1.0 / 8.0) * fe_0 * fe_0 * fe_0 - 4.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y * rpc_y - 4.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_x * rpc_x - rpa_z * rpa_z * rpa_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x - rpa_z * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x - rpb_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_x + fe_0 * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_xyy[i] += fss * b4_vals[i] * (fe_0 * rpb_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (4.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpa_z * rpa_z * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyy[i] += fss * b4_vals[i] * (2.0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpb_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xyz[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_xyz[i] += fss * b0_vals[i] * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x;

        fints_xyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_x * rpb_y * rpb_x - fe_0 * rpa_z * rpa_x * rpb_y * rpc_x - fe_0 * rpa_z * rpa_x * rpb_x * rpc_y - fe_0 * rpa_z * rpb_z * rpb_y * rpc_z - fe_0 * rpa_z * rpb_y * rpb_x * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_x - fe_0 * fe_0 * rpa_z * rpb_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_x * rpc_z - rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * rpc_y - rpa_z * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z - rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_y * rpc_x + fe_0 * rpa_z * rpa_x * rpb_x * rpc_y + fe_0 * rpa_z * rpa_x * rpc_y * rpc_x + fe_0 * rpa_z * rpb_z * rpb_y * rpc_z + fe_0 * rpa_z * rpb_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpb_y * rpb_x * rpc_x + fe_0 * rpa_z * rpb_y * rpc_z * rpc_z + fe_0 * rpa_z * rpb_y * rpc_x * rpc_x + fe_0 * rpa_z * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y + fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpa_x * rpb_y * rpb_x * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_x + rpa_z * rpa_z * rpa_x * rpb_z * rpc_y * rpc_x + rpa_z * rpa_z * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (rpa_z * rpa_z * rpa_x * rpb_x * rpc_z * rpc_y + rpa_z * rpa_z * rpb_z * rpb_y * rpc_x * rpc_x + rpa_z * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x + rpa_z * rpa_z * rpb_y * rpb_x * rpc_z * rpc_x + rpa_x * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_x * rpc_y * rpc_x - fe_0 * rpa_z * rpb_z * rpc_z * rpc_y - fe_0 * rpa_z * rpb_y * rpc_z * rpc_z - fe_0 * rpa_z * rpb_y * rpc_x * rpc_x - fe_0 * rpa_z * rpb_x * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpc_z * rpc_z * rpc_y - fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_xyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x - rpa_z * rpa_z * rpa_x * rpc_z * rpc_y * rpc_x - rpa_z * rpa_z * rpb_z * rpc_y * rpc_x * rpc_x - rpa_z * rpa_z * rpb_y * rpc_z * rpc_x * rpc_x - rpa_z * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x - rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - rpa_x * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z - rpb_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (fe_0 * rpa_z * rpc_z * rpc_z * rpc_y + fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpa_z * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x + rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x + rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y + rpb_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_xzz[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x + fe_0 * fe_0 * rpa_z * rpb_z);

        fints_xzz[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_x - 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x - 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z - 2.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x - fe_0 * rpa_z * rpb_z * rpb_z * rpc_z);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_x - fe_0 * rpa_z * rpa_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_x - 2.0 * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_xzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z);

        fints_xzz[i] += fss * b1_vals[i] * (-rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpc_x - rpa_z * rpa_z * rpb_z * rpb_z * rpb_x * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_x + 3.0 * fe_0 * rpa_z * rpa_x * rpb_x * rpc_z + 3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x + 2.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x + 2.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x + fe_0 * rpa_z * rpb_z * rpb_z * rpc_z + 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_x + fe_0 * rpa_z * rpa_z * rpb_z * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_x * rpb_z * rpb_x * rpc_z + 3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_x + 3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z + 3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x + fe_0 * fe_0 * rpa_z * rpb_z + 3.0 * fe_0 * fe_0 * rpa_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_x);

        fints_xzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_x + 3.0 * fe_0 * fe_0 * rpb_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 4.0 * rpa_z * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z + 2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_x + 2.0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpc_z * rpc_x + rpa_z * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z + 2.0 * rpa_z * rpa_z * rpb_z * rpb_x * rpc_z * rpc_x + rpa_z * rpa_z * rpb_z * rpb_z * rpc_x * rpc_x + rpa_x * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_x - 2.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z - 2.0 * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpc_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_x * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x - fe_0 * rpb_z * rpc_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_x * rpc_x - 3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_x - 3.0 * fe_0 * fe_0 * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 4.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_x * rpc_z * rpc_z * rpc_z - 4.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_x * rpc_x - rpa_z * rpa_z * rpa_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-rpa_z * rpa_z * rpb_x * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z - rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x - rpb_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + fe_0 * rpa_z * rpc_z * rpc_z * rpc_z + 3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x + fe_0 * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_xzz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_x + 4.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x + rpa_z * rpa_z * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (rpa_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z + 2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x + rpb_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z - 2.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - 2.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x;

        fints_yyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_yyy[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpb_y * rpc_z - 3.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpc_y - rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z + 3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y + 3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_yyy[i] += fss * b2_vals[i] * (6.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_x + 3.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpc_y * rpc_y + 3.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x + rpa_x * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z);

        fints_yyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 6.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-6.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - rpa_z * rpa_z * rpa_x * rpc_y * rpc_y * rpc_y - 3.0 * rpa_z * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - rpb_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y * rpc_y + 6.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b4_vals[i] * (rpa_z * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yyy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpa_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyy[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_yyz[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z);

        fints_yyz[i] += fss * b0_vals[i] * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y;

        fints_yyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_x * rpb_z * rpc_z - 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y - fe_0 * rpa_z * rpa_x * rpb_y * rpb_y - fe_0 * rpa_z * rpb_y * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_x - fe_0 * rpa_x * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_x - fe_0 * fe_0 * rpa_z * rpa_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - 2.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpb_y * rpc_z - 2.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpc_y - rpa_z * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z - rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_x * rpb_z * rpc_z + 2.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_y + fe_0 * rpa_z * rpa_x * rpc_z * rpc_z + fe_0 * rpa_z * rpa_x * rpc_y * rpc_y + fe_0 * rpa_z * rpb_z * rpc_z * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x + fe_0 * rpa_z * rpb_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (fe_0 * rpa_x * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y + 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_x * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x + fe_0 * fe_0 * rpa_z * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 4.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_x * rpb_y * rpb_y * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_x + rpa_z * rpa_z * rpa_x * rpb_z * rpc_y * rpc_y + 2.0 * rpa_z * rpa_z * rpa_x * rpb_y * rpc_z * rpc_y + 2.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b2_vals[i] * (rpa_z * rpa_z * rpb_y * rpb_y * rpc_z * rpc_x + rpa_x * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_x * rpc_z * rpc_z - fe_0 * rpa_z * rpa_x * rpc_y * rpc_y - fe_0 * rpa_z * rpb_z * rpc_z * rpc_x - 2.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x - fe_0 * rpa_z * rpc_z * rpc_z * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z - fe_0 * rpb_z * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 2.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y * rpc_y - 4.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_y - 4.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_yyz[i] += fss * b3_vals[i] * (-rpa_z * rpa_z * rpa_x * rpc_z * rpc_y * rpc_y - rpa_z * rpa_z * rpb_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - rpa_x * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-rpb_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * (fe_0 * rpa_z * rpc_z * rpc_z * rpc_x + fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_yyz[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + 4.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpa_z * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x + rpa_x * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (2.0 * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y + 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_yyz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_yyz[i] += fss * b5_vals[i] * (-2.0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_yzz[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_y - 2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y - 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z - 2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_x - 3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - 2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_yzz[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z - rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpc_y - rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_y + 3.0 * fe_0 * rpa_z * rpa_x * rpb_y * rpc_z + 3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y + 2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x + 2.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpb_z * rpb_y * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_y + 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_yzz[i] += fss * b2_vals[i] * (4.0 * rpa_z * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z + 2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_y + 2.0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_x + 2.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpc_z * rpc_y + rpa_z * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpc_z * rpc_x + rpa_z * rpa_z * rpb_z * rpb_z * rpc_y * rpc_x + rpa_x * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_y - 2.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_x * rpb_y * rpc_z * rpc_z - 3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x - 3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_yzz[i] += fss * b3_vals[i] * (-4.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpa_x * rpb_y * rpc_z * rpc_z * rpc_z - 4.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - rpa_z * rpa_z * rpa_x * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x - rpa_z * rpa_z * rpb_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_x * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z - rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_y + 4.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x + rpa_z * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b4_vals[i] * (2.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + rpa_x * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z + 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x + rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_yzz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y - 2.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_yzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_zzz[i] += fss * b0_vals[i] * (3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z + (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z);

        fints_zzz[i] += fss * b0_vals[i] * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z;

        fints_zzz[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z - 3.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z - 3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_x - 3.0 * fe_0 * fe_0 * rpa_z * rpa_x);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x - (9.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpb_z - (15.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - 2.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_zzz[i] += fss * b1_vals[i] * (-3.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z - rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpc_z + 6.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z + 9.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x + 3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_x * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x + 9.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_x * rpb_z * rpb_z * rpc_z + (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_zzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_x + 3.0 * fe_0 * fe_0 * rpa_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpb_z + (15.0 / 2.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 6.0 * rpa_z * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_x + 3.0 * rpa_z * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * (3.0 * rpa_z * rpa_z * rpb_z * rpb_z * rpc_z * rpc_x + rpa_x * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpa_z * rpa_x * rpc_z * rpc_z - 9.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x - 6.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x - 9.0 * fe_0 * rpa_x * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-5.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z - 9.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpa_x * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 6.0 * rpa_z * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z - 6.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x - rpa_z * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_z * rpc_x - 3.0 * rpa_x * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * (6.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x + 5.0 * fe_0 * rpa_x * rpc_z * rpc_z * rpc_z + 9.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x + 5.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_zzz[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_x * rpc_z * rpc_z * rpc_z * rpc_z + 6.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + rpa_z * rpa_z * rpc_z * rpc_z * rpc_z * rpc_x + 3.0 * rpa_x * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z + 3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_zzz[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - rpa_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z - 3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_zzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x;

    }
}

auto
compPrimitiveNuclearPotentialFF_YYY_T(      TDoubleArray& buffer_xxx,
                                            TDoubleArray& buffer_xxy,
                                            TDoubleArray& buffer_xxz,
                                            TDoubleArray& buffer_xyy,
                                            TDoubleArray& buffer_xyz,
                                            TDoubleArray& buffer_xzz,
                                            TDoubleArray& buffer_yyy,
                                            TDoubleArray& buffer_yyz,
                                            TDoubleArray& buffer_yzz,
                                            TDoubleArray& buffer_zzz,
                       const double charge,
                       const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

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

    // set up Boys function variables

    const CBoysFunc<6> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<7> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<7>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_xxx,\
                             fints_xxy,\
                             fints_xxz,\
                             fints_xyy,\
                             fints_xyz,\
                             fints_xzz,\
                             fints_yyy,\
                             fints_yyz,\
                             fints_yzz,\
                             fints_zzz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxx[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - 3.0 * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * rpc_y);

        fints_xxx[i] += fss * b1_vals[i] * (-3.0 * rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 3.0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_y * rpc_y + 9.0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x + 3.0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xxx[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 9.0 * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x - 9.0 * rpa_y * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x - rpa_y * rpa_y * rpa_y * rpc_x * rpc_x * rpc_x - rpb_x * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xxx[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 9.0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x + 3.0 * rpa_y * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xxy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_xxy[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xxy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y - 3.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpc_x);

        fints_xxy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xxy[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x);

        fints_xxy[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 3.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpc_x - rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x + 9.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_y + 3.0 * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y + 3.0 * fe_0 * fe_0 * rpb_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 3.0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x + 3.0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (rpa_y * rpa_y * rpa_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x - 9.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x * rpc_x - 3.0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x - 3.0 * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * fe_0 * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x - 3.0 * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-3.0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x - 6.0 * rpa_y * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x - rpa_y * rpa_y * rpa_y * rpc_y * rpc_x * rpc_x - rpb_y * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xxy[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y + 6.0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + 3.0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (6.0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + rpb_x * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xxy[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - rpb_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xxz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_xxz[i] += fss * b1_vals[i] * (-3.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpc_x - rpa_y * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z + 3.0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * (3.0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x + 3.0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpa_y * rpb_z * rpc_x * rpc_x + 2.0 * rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 6.0 * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-3.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x - 6.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x - rpa_y * rpa_y * rpa_y * rpc_z * rpc_x * rpc_x - rpb_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xxz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 3.0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + 6.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (3.0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xxz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xxz[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xyy[i] += fss * b0_vals[i] * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x;

        fints_xyy[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x - 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - 3.0 * fe_0 * fe_0 * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - 3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xyy[i] += fss * b1_vals[i] * (-2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpc_y - rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_y + 9.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x + 9.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_x + 6.0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x + 3.0 * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * ((15.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 3.0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y * rpc_x + rpa_y * rpa_y * rpa_y * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_x - 9.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x - 6.0 * fe_0 * rpb_y * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x - 5.0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(15.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - (15.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 6.0 * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x - 6.0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-3.0 * rpa_y * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y - rpa_y * rpa_y * rpa_y * rpc_y * rpc_y * rpc_x - rpb_y * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b4_vals[i] * (9.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x + 6.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x + 5.0 * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y + 5.0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x + (15.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (6.0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y + rpb_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - rpb_x * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_xyz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x + rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - 3.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * rpc_y - rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_x - rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpc_y - rpa_y * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x + 3.0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 3.0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x + 3.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * (3.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpa_y * rpb_z * rpc_y * rpc_x + rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_x + rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpb_z * rpb_x * rpc_y * rpc_y - 3.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 3.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-3.0 * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x - 3.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x - 3.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_y * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x - rpb_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 3.0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + 3.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + 3.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + rpb_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y + rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - 3.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_xzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xzz[i] += fss * b1_vals[i] * (-3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z - rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z + 3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_x + 3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (3.0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y + 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_x + rpa_y * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x - 3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 6.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-3.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x - 6.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x - 3.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y - rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_x - rpb_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 6.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + 3.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] * (3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - 3.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_yyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y + (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_yyy[i] += fss * b0_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(27.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y - (27.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_y - (27.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y - (45.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(45.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y - (45.0 / 8.0) * fe_0 * fe_0 * fe_0 - 3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpc_y - 3.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (27.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y + (27.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y + (27.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y + 9.0 * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_y + (27.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (45.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_yyy[i] += fss * b2_vals[i] * ((45.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (45.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 + 3.0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (9.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-27.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y - 15.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y - 15.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y - 9.0 * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-(45.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (45.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (45.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y - (15.0 / 8.0) * fe_0 * fe_0 * fe_0 - 9.0 * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-9.0 * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y - rpa_y * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] * (15.0 * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y + 15.0 * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y + (15.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y + (45.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + 9.0 * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] * (3.0 * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b5_vals[i] * (-(15.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b6_vals[i] * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yyz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y + 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yyz[i] += fss * b0_vals[i] * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y;

        fints_yyz[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z - 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y);

        fints_yyz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_yyz[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - 3.0 * fe_0 * fe_0 * rpb_z * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z - 3.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_yyz[i] += fss * b1_vals[i] * (-2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y - rpa_y * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y + 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z + 6.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (3.0 * fe_0 * fe_0 * rpb_y * rpc_z + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 3.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (rpa_y * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y + 2.0 * rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y - 9.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y - 6.0 * fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-5.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y - 6.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (15.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 6.0 * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-6.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - rpa_y * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y - rpb_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (9.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + 5.0 * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y + 6.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + 5.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (3.0 * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpb_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - 2.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_yzz[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y - (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z);

        fints_yzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - 2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z + 9.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y + 3.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + 3.0 * fe_0 * fe_0 * rpb_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 3.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 6.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + rpa_y * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y - 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - 6.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - 3.0 * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - 3.0 * fe_0 * fe_0 * rpc_y * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-6.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y + 6.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + 6.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (3.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + rpb_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 2.0 * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_zzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_zzz[i] += fss * b1_vals[i] * (-3.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpa_y * rpc_z + (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 3.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_y * rpc_y + 9.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y + 3.0 * rpa_y * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_zzz[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_zzz[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 9.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y - 9.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y - rpa_y * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 9.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + 3.0 * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - 3.0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - 3.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_zzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialFF_YYZ_T(      TDoubleArray& buffer_xxx,
                                            TDoubleArray& buffer_xxy,
                                            TDoubleArray& buffer_xxz,
                                            TDoubleArray& buffer_xyy,
                                            TDoubleArray& buffer_xyz,
                                            TDoubleArray& buffer_xzz,
                                            TDoubleArray& buffer_yyy,
                                            TDoubleArray& buffer_yyz,
                                            TDoubleArray& buffer_yzz,
                                            TDoubleArray& buffer_zzz,
                       const double charge,
                       const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

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

    // set up Boys function variables

    const CBoysFunc<6> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<7> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<7>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_xxx,\
                             fints_xxy,\
                             fints_xxz,\
                             fints_xyy,\
                             fints_xyz,\
                             fints_xzz,\
                             fints_yyy,\
                             fints_yyz,\
                             fints_yzz,\
                             fints_zzz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxx[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xxx[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_y - 3.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpc_x - rpa_y * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxx[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y + 3.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (6.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y * rpc_x + 3.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_x * rpc_x + rpa_z * rpb_x * rpb_x * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y + 3.0 * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 6.0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-rpa_z * rpa_y * rpa_y * rpc_x * rpc_x * rpc_x - 3.0 * rpa_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_x - 6.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - 3.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x - rpb_x * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xxx[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_y * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (6.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - rpa_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xxy[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y);

        fints_xxy[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpb_x;

        fints_xxy[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_y * rpb_y * rpc_y - 2.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x - fe_0 * rpa_z * rpa_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y - fe_0 * rpa_y * rpb_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z - fe_0 * fe_0 * rpa_z * rpa_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_y - 2.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpc_x - rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpc_y - rpa_y * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_y * rpc_y + 2.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x + fe_0 * rpa_z * rpa_y * rpc_y * rpc_y + fe_0 * rpa_z * rpa_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + 2.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x + fe_0 * rpa_y * rpb_x * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpb_y * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y);

        fints_xxy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y + fe_0 * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 4.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y * rpc_y + rpa_z * rpa_y * rpa_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_y * rpc_x + rpa_z * rpb_y * rpb_x * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x + rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_y * rpc_y * rpc_y - fe_0 * rpa_z * rpa_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y - fe_0 * rpa_y * rpb_y * rpc_z * rpc_y - 2.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x - fe_0 * rpa_y * rpc_z * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y - fe_0 * rpb_y * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 2.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x * rpc_x - 4.0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_x - rpa_z * rpa_y * rpa_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-rpa_z * rpb_x * rpb_x * rpc_y * rpc_y * rpc_y - 4.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y - rpa_y * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x - 2.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-rpb_y * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xxy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y + fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_xxy[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xxy[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x * rpc_x + rpa_z * rpb_y * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 4.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (rpa_y * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpb_x * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xxy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - rpa_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpb_y * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-2.0 * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xxz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_xxz[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_y * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z - fe_0 * rpa_z * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z - fe_0 * rpa_y * rpb_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z - fe_0 * rpa_y * rpa_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_y);

        fints_xxz[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpc_x - rpa_z * rpa_y * rpa_y * rpb_x * rpb_x * rpc_z - rpa_y * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_z * rpc_y + fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z + fe_0 * rpa_z * rpb_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x + fe_0 * rpa_z * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z + fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + 2.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_x * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z + fe_0 * rpa_y * rpa_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_x * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z + fe_0 * fe_0 * rpa_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z + fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_xxz[i] += fss * b2_vals[i] * (4.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y + rpa_z * rpa_y * rpa_y * rpb_z * rpc_x * rpc_x + 2.0 * rpa_z * rpa_y * rpa_y * rpb_x * rpc_z * rpc_x + rpa_z * rpb_z * rpb_x * rpb_x * rpc_y * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y + 2.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x + rpa_y * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_y * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x - fe_0 * rpa_z * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x - fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - 2.0 * fe_0 * rpa_y * rpb_x * rpc_y * rpc_x - fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - fe_0 * rpa_y * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_x * rpc_x - fe_0 * rpb_z * rpb_x * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-fe_0 * rpb_x * rpc_z * rpc_z * rpc_x - fe_0 * rpb_x * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (1.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_x * rpc_x - 4.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x - rpa_z * rpa_y * rpa_y * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_x - rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y - 4.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y - rpa_y * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x - rpb_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xxz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + fe_0 * rpa_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_xxz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x + fe_0 * rpb_x * rpc_z * rpc_z * rpc_x + fe_0 * rpb_x * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + 2.0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (rpa_z * rpb_z * rpc_y * rpc_y * rpc_x * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + 4.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpa_y * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (2.0 * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xxz[i] += fss * b5_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-rpb_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xyy[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x - 2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x - 2.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - 2.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xyy[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpc_y - rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpc_x - rpa_y * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_y + 3.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_x + 3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y + 2.0 * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z + 2.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x + 3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (4.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y * rpc_y + 2.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_y * rpc_x + rpa_z * rpa_y * rpa_y * rpb_x * rpc_y * rpc_y + rpa_z * rpb_y * rpb_y * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y + 2.0 * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - 2.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x - 3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-4.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y * rpc_y - rpa_z * rpa_y * rpa_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y * rpc_y - rpa_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-4.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x - rpa_y * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y - rpb_y * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_x + rpa_z * rpb_x * rpc_y * rpc_y * rpc_y * rpc_y + 4.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpb_x * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_xyz[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xyz[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_x;

        fints_xyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_y * rpb_z * rpb_x - fe_0 * rpa_z * rpa_y * rpb_z * rpc_x - fe_0 * rpa_z * rpa_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z - fe_0 * rpa_y * rpb_z * rpb_x * rpc_z - fe_0 * rpa_y * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z - fe_0 * fe_0 * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_y - rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpc_y - rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * rpc_z - rpa_y * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_z * rpc_x + fe_0 * rpa_z * rpa_y * rpb_x * rpc_z + fe_0 * rpa_z * rpa_y * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y + fe_0 * rpa_y * rpb_z * rpb_x * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_z * rpc_z * rpc_x + fe_0 * rpa_y * rpb_y * rpb_x * rpc_y + fe_0 * rpa_y * rpb_y * rpc_y * rpc_x + fe_0 * rpa_y * rpb_x * rpc_z * rpc_z + fe_0 * rpa_y * rpb_x * rpc_y * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x + fe_0 * fe_0 * rpa_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y * rpc_y + 2.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y + rpa_z * rpa_y * rpa_y * rpb_z * rpc_y * rpc_x + rpa_z * rpa_y * rpa_y * rpb_y * rpc_z * rpc_x + rpa_z * rpa_y * rpa_y * rpb_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * (rpa_z * rpb_z * rpb_y * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x + rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y + rpa_y * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpb_z * rpc_z * rpc_x - fe_0 * rpa_y * rpb_y * rpc_y * rpc_x - fe_0 * rpa_y * rpb_x * rpc_z * rpc_z - fe_0 * rpa_y * rpb_x * rpc_y * rpc_y - fe_0 * rpa_y * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y - rpa_z * rpa_y * rpa_y * rpc_z * rpc_y * rpc_x - rpa_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_z * rpb_z * rpb_x * rpc_y * rpc_y * rpc_y - rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x - rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x - rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y - rpb_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + fe_0 * rpa_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_x + rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y + rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpb_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_xzz[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x + fe_0 * rpa_y * rpa_y * rpb_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xzz[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_x;

        fints_xzz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_x - fe_0 * rpa_z * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x - 2.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y - fe_0 * rpa_y * rpa_y * rpb_z * rpb_x - fe_0 * rpa_y * rpa_y * rpb_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - fe_0 * fe_0 * rpb_z * rpb_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - 2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_y - 2.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z - rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpc_x - rpa_y * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_x * rpc_y + fe_0 * rpa_z * rpa_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_x + fe_0 * rpa_z * rpb_z * rpb_x * rpc_z + fe_0 * rpa_z * rpb_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y + 2.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y + 2.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y + fe_0 * rpa_y * rpa_y * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x + fe_0 * rpb_z * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x + fe_0 * fe_0 * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 4.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_z * rpc_x + rpa_z * rpa_y * rpa_y * rpb_x * rpc_z * rpc_z + rpa_z * rpb_z * rpb_z * rpb_x * rpc_y * rpc_y + 2.0 * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z + rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_y * rpc_y * rpc_x - fe_0 * rpa_z * rpb_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - 2.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-fe_0 * rpb_z * rpb_x * rpc_z * rpc_z - fe_0 * rpb_z * rpb_x * rpc_y * rpc_y - fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - fe_0 * rpb_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 4.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y - rpa_z * rpa_y * rpa_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xzz[i] += fss * b3_vals[i] * (-rpa_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_x - 4.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x - rpa_y * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-rpb_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + fe_0 * rpb_z * rpc_z * rpc_z * rpc_x + fe_0 * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y + 4.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] * (rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x + 2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y + rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - 2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_yyy[i] += fss * b0_vals[i] * (3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y);

        fints_yyy[i] += fss * b0_vals[i] * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpb_y;

        fints_yyy[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y - 3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y - 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z - 3.0 * fe_0 * fe_0 * rpa_z * rpa_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_y);

        fints_yyy[i] += fss * b1_vals[i] * (-3.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpc_y - rpa_y * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z);

        fints_yyy[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y + 6.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y + 9.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + (15.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y + 3.0 * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_yyy[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 6.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y * rpc_y + 3.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_y * rpc_y + rpa_z * rpb_y * rpb_y * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y + 3.0 * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpa_z * rpa_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y - 5.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y - 9.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y - 6.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y - 9.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_yyy[i] += fss * b3_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 6.0 * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y * rpc_y - rpa_z * rpa_y * rpa_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-6.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 3.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - rpb_y * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] * (5.0 * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y + 6.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + 9.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + 5.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y * rpc_y + 6.0 * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + rpa_y * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - rpa_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y - 2.0 * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b6_vals[i] * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yyz[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_yyz[i] += fss * b0_vals[i] * (fe_0 * fe_0 * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y - 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y - 2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z);

        fints_yyz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z - 2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z - fe_0 * rpa_y * rpb_y * rpb_y * rpc_y);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z - fe_0 * rpa_y * rpa_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_yyz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - 2.0 * fe_0 * fe_0 * rpa_y * rpb_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_yyz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_y - 2.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpc_y);

        fints_yyz[i] += fss * b1_vals[i] * (-rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * rpc_z - rpa_y * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y + 2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z + 3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z + 3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y + 3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z + 2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z + 3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + 2.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y + fe_0 * rpa_y * rpb_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z + fe_0 * rpa_y * rpa_y * rpb_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z + fe_0 * fe_0 * rpa_y * rpb_y + 3.0 * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 4.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y * rpc_y + 2.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y + rpa_z * rpa_y * rpa_y * rpb_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_y * rpa_y * rpb_y * rpc_z * rpc_y + rpa_z * rpb_z * rpb_y * rpb_y * rpc_y * rpc_y + 2.0 * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + 2.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + rpa_y * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-2.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z - 2.0 * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - fe_0 * rpa_y * rpc_y * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_y * rpc_y - 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - 3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - fe_0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - 3.0 * fe_0 * fe_0 * rpc_y * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-4.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - rpa_z * rpa_y * rpa_y * rpc_z * rpc_y * rpc_y - 2.0 * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y * rpc_y - rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y - 4.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-2.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - rpa_y * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - rpb_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + fe_0 * rpa_y * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y + 3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (fe_0 * rpb_y * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y + rpa_z * rpb_z * rpc_y * rpc_y * rpc_y * rpc_y + 2.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + 2.0 * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 4.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (rpa_y * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y - rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - 2.0 * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - rpb_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b5_vals[i] * (-2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yzz[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y + fe_0 * rpa_y * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y);

        fints_yzz[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + fe_0 * fe_0 * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y + rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z - fe_0 * rpa_z * rpa_y * rpb_z * rpb_z - fe_0 * rpa_z * rpa_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y);

        fints_yzz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y - 2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y - fe_0 * rpa_y * rpb_z * rpb_z * rpc_z);

        fints_yzz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpa_y * rpb_z * rpb_y - fe_0 * rpa_y * rpa_y * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z - fe_0 * fe_0 * rpa_z * rpa_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - 2.0 * fe_0 * fe_0 * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z - fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_y - 2.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z - rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpc_y);

        fints_yzz[i] += fss * b1_vals[i] * (-rpa_y * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z + fe_0 * rpa_z * rpa_y * rpb_y * rpc_y + fe_0 * rpa_z * rpa_y * rpc_z * rpc_z + fe_0 * rpa_z * rpa_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpb_z * rpb_y * rpc_z + 3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_y + 2.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z + 2.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y + fe_0 * rpa_y * rpb_z * rpb_z * rpc_z + 3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpa_y * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y + fe_0 * rpb_z * rpb_y * rpc_z * rpc_z + fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (fe_0 * fe_0 * rpa_y * rpb_z + 3.0 * fe_0 * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y + 3.0 * fe_0 * fe_0 * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 4.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_y * rpc_y + 2.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_z * rpc_y + rpa_z * rpa_y * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * (rpa_z * rpb_z * rpb_z * rpb_y * rpc_y * rpc_y + 2.0 * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y + 2.0 * rpa_y * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z + rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_y * rpc_z * rpc_z - fe_0 * rpa_z * rpa_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y - 2.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z - 2.0 * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - fe_0 * rpa_y * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpc_z * rpc_y - fe_0 * rpb_z * rpb_y * rpc_z * rpc_z - fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y - fe_0 * rpb_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-4.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - rpa_z * rpa_y * rpa_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - rpa_z * rpb_z * rpb_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-4.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y - rpa_y * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + fe_0 * rpa_y * rpc_z * rpc_z * rpc_z + 3.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (fe_0 * rpb_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 4.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y + rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y + 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + rpb_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 2.0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - 2.0 * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yzz[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_zzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y);

        fints_zzz[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpa_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_y - 3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z - rpa_y * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y + 3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_y * rpa_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z + 9.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_y * rpa_y * rpb_z * rpc_z + 3.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * fe_0 * rpa_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpa_y + (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 6.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y + 3.0 * rpa_z * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z + rpa_z * rpb_z * rpb_z * rpb_z * rpc_y * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] * (2.0 * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y + 3.0 * rpa_y * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - 6.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpa_y * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * fe_0 * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y - rpa_z * rpa_y * rpa_y * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y - 6.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - 3.0 * rpa_y * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_zzz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z + 6.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + 2.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b4_vals[i] * (3.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + 6.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + rpa_y * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z + 3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_zzz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z - rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y - 3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_zzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialFF_YZZ_T(      TDoubleArray& buffer_xxx,
                                            TDoubleArray& buffer_xxy,
                                            TDoubleArray& buffer_xxz,
                                            TDoubleArray& buffer_xyy,
                                            TDoubleArray& buffer_xyz,
                                            TDoubleArray& buffer_xzz,
                                            TDoubleArray& buffer_yyy,
                                            TDoubleArray& buffer_yyz,
                                            TDoubleArray& buffer_yzz,
                                            TDoubleArray& buffer_zzz,
                       const double charge,
                       const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

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

    // set up Boys function variables

    const CBoysFunc<6> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<7> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<7>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_xxx,\
                             fints_xxy,\
                             fints_xxz,\
                             fints_xyy,\
                             fints_xyz,\
                             fints_xzz,\
                             fints_yyy,\
                             fints_yyz,\
                             fints_yzz,\
                             fints_zzz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxx[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xxx[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpb_x * rpc_z - 3.0 * rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpc_x - rpa_z * rpa_z * rpb_x * rpb_x * rpb_x * rpc_y);

        fints_xxx[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z + 3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_y);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * (6.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_y + 3.0 * rpa_z * rpa_z * rpa_y * rpb_x * rpc_x * rpc_x + 3.0 * rpa_z * rpa_z * rpb_x * rpb_x * rpc_y * rpc_x + rpa_y * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxx[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 6.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-6.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_x - rpa_z * rpa_z * rpa_y * rpc_x * rpc_x * rpc_x - 3.0 * rpa_z * rpa_z * rpb_x * rpc_y * rpc_x * rpc_x - 3.0 * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x - rpb_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxx[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_x * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpc_z * rpc_x * rpc_x * rpc_x + 6.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * (rpa_z * rpa_z * rpc_y * rpc_x * rpc_x * rpc_x + 3.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_x * rpc_x * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxx[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x * rpc_x;

        fints_xxy[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_xxy[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (1.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x);

        fints_xxy[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_y * rpb_y * rpc_z - fe_0 * rpa_z * rpb_x * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x - fe_0 * rpa_y * rpb_y * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxy[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpc_x - rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpc_y - rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_y * rpc_z + fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + fe_0 * rpa_z * rpb_y * rpc_z * rpc_y + 2.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x + fe_0 * rpa_z * rpb_x * rpb_x * rpc_z);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_y + fe_0 * rpa_z * rpa_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_y * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x + fe_0 * rpa_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (fe_0 * rpb_y * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y + fe_0 * fe_0 * rpa_z * rpc_z);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y + fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0);

        fints_xxy[i] += fss * b2_vals[i] * (4.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_y + rpa_z * rpa_z * rpa_y * rpb_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_z * rpa_y * rpb_x * rpc_y * rpc_x);

        fints_xxy[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_z * rpb_y * rpb_x * rpc_y * rpc_x + rpa_z * rpa_z * rpb_x * rpb_x * rpc_y * rpc_y + rpa_y * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxy[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_y * rpc_z * rpc_y - fe_0 * rpa_z * rpb_y * rpc_z * rpc_y - 2.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x - fe_0 * rpa_z * rpc_z * rpc_y * rpc_y - fe_0 * rpa_z * rpc_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_x * rpc_x - fe_0 * rpa_y * rpb_x * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x - fe_0 * rpb_y * rpb_x * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-fe_0 * rpb_x * rpc_z * rpc_z * rpc_x - fe_0 * rpb_x * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_xxy[i] += fss * b3_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_xxy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (1.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_x * rpc_x - 4.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_x - 4.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y * rpc_y - rpa_z * rpa_z * rpa_y * rpc_y * rpc_x * rpc_x - rpa_z * rpa_z * rpb_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpa_z * rpb_x * rpc_y * rpc_y * rpc_x - 2.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y - rpb_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxy[i] += fss * b4_vals[i] * (fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_xxy[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_x * rpc_x + fe_0 * rpb_x * rpc_z * rpc_z * rpc_x + fe_0 * rpb_x * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (1.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + 2.0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x * rpc_x + 4.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_x + rpa_z * rpa_z * rpc_y * rpc_y * rpc_x * rpc_x + rpa_y * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (2.0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xxy[i] += fss * b5_vals[i] * (-(1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xxy[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x * rpc_x;

        fints_xxz[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_xxz[i] += fss * b0_vals[i] * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x;

        fints_xxz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_y * rpb_z * rpc_z - 2.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x - fe_0 * rpa_z * rpa_y * rpb_x * rpb_x - fe_0 * rpa_z * rpb_x * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_y - fe_0 * rpa_y * rpb_z * rpb_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_y - fe_0 * fe_0 * rpa_z * rpa_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_xxz[i] += fss * b1_vals[i] * (-(1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - 2.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * rpc_x - rpa_z * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z - rpa_z * rpa_z * rpb_z * rpb_x * rpb_x * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_z * rpc_z + 2.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_x + fe_0 * rpa_z * rpa_y * rpc_z * rpc_z + fe_0 * rpa_z * rpa_y * rpc_x * rpc_x + fe_0 * rpa_z * rpb_z * rpc_z * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x + fe_0 * rpa_z * rpb_x * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_z * rpb_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x + 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * (fe_0 * rpb_z * rpb_x * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y + fe_0 * fe_0 * rpa_z * rpc_y);

        fints_xxz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 4.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_y * rpb_x * rpb_x * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_y + rpa_z * rpa_z * rpa_y * rpb_z * rpc_x * rpc_x + 2.0 * rpa_z * rpa_z * rpa_y * rpb_x * rpc_z * rpc_x + 2.0 * rpa_z * rpa_z * rpb_z * rpb_x * rpc_y * rpc_x);

        fints_xxz[i] += fss * b2_vals[i] * (rpa_z * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y + rpa_y * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_y * rpc_z * rpc_z - fe_0 * rpa_z * rpa_y * rpc_x * rpc_x - fe_0 * rpa_z * rpb_z * rpc_z * rpc_y - 2.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x - fe_0 * rpa_z * rpc_z * rpc_z * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_x * rpc_x - 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z - fe_0 * rpb_z * rpb_x * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 2.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x * rpc_x - 4.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_x - 4.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxz[i] += fss * b3_vals[i] * (-rpa_z * rpa_z * rpa_y * rpc_z * rpc_x * rpc_x - rpa_z * rpa_z * rpb_z * rpc_y * rpc_x * rpc_x - 2.0 * rpa_z * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x - 2.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x - rpa_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-rpb_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xxz[i] += fss * b4_vals[i] * (fe_0 * rpa_z * rpc_z * rpc_z * rpc_y + fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_xxz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_x * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xxz[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x * rpc_x + 4.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpa_z * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x + rpa_y * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x + 2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xxz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpa_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - rpb_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x);

        fints_xxz[i] += fss * b5_vals[i] * (-2.0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_xyy[i] += fss * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x + fe_0 * rpa_z * rpa_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x);

        fints_xyy[i] += fss * b0_vals[i] * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x;

        fints_xyy[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_y * rpb_x * rpc_z - 2.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_x - fe_0 * rpa_z * rpa_z * rpb_y * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_z * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_y - fe_0 * rpa_y * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - fe_0 * fe_0 * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_xyy[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - 2.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpc_y - rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpc_x - rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_x * rpc_z + fe_0 * rpa_z * rpa_y * rpc_z * rpc_x + 2.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z + 2.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_x + fe_0 * rpa_z * rpa_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x + fe_0 * rpa_y * rpb_y * rpb_x * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpa_y * rpb_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y + fe_0 * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyy[i] += fss * b2_vals[i] * (fe_0 * rpb_y * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x + fe_0 * fe_0 * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 4.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_x + 2.0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpc_y * rpc_x + rpa_z * rpa_z * rpa_y * rpb_x * rpc_y * rpc_y + 2.0 * rpa_z * rpa_z * rpb_y * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * (rpa_z * rpa_z * rpb_y * rpb_y * rpc_y * rpc_x + rpa_y * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyy[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_y * rpc_z * rpc_x - 2.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-fe_0 * rpa_y * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-fe_0 * rpb_y * rpb_x * rpc_z * rpc_z - fe_0 * rpb_y * rpb_x * rpc_y * rpc_y - fe_0 * rpb_y * rpc_z * rpc_z * rpc_x - fe_0 * rpb_y * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y - (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 4.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y * rpc_y - 4.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y * rpc_y - 2.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-rpa_z * rpa_z * rpa_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpa_z * rpb_y * rpc_y * rpc_y * rpc_x - rpa_z * rpa_z * rpb_x * rpc_y * rpc_y * rpc_y - 2.0 * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y - rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-rpb_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_x + fe_0 * rpb_y * rpc_z * rpc_z * rpc_x + fe_0 * rpb_y * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_x + 4.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y * rpc_y + rpa_z * rpa_z * rpc_y * rpc_y * rpc_y * rpc_x + 2.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b4_vals[i] * (rpa_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y + rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_x - 2.0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - 2.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-rpb_x * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_xyy[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_x;

        fints_xyz[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xyz[i] += fss * b0_vals[i] * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x;

        fints_xyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_y * rpb_y * rpb_x - fe_0 * rpa_z * rpa_y * rpb_y * rpc_x - fe_0 * rpa_z * rpa_y * rpb_x * rpc_y - fe_0 * rpa_z * rpb_z * rpb_x * rpc_z - fe_0 * rpa_z * rpb_y * rpb_x * rpc_y);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_y - fe_0 * fe_0 * rpa_z * rpb_x - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_x * rpc_z - rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * rpc_y - rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z - rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_y * rpc_x + fe_0 * rpa_z * rpa_y * rpb_x * rpc_y + fe_0 * rpa_z * rpa_y * rpc_y * rpc_x + fe_0 * rpa_z * rpb_z * rpb_x * rpc_z + fe_0 * rpa_z * rpb_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpb_y * rpb_x * rpc_y + fe_0 * rpa_z * rpb_y * rpc_y * rpc_x + fe_0 * rpa_z * rpb_x * rpc_z * rpc_z + fe_0 * rpa_z * rpb_x * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x + fe_0 * fe_0 * rpa_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpb_y * rpb_x * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_y + rpa_z * rpa_z * rpa_y * rpb_z * rpc_y * rpc_x + rpa_z * rpa_z * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * (rpa_z * rpa_z * rpa_y * rpb_x * rpc_z * rpc_y + rpa_z * rpa_z * rpb_z * rpb_y * rpc_y * rpc_x + rpa_z * rpa_z * rpb_z * rpb_x * rpc_y * rpc_y + rpa_z * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y + rpa_y * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_y * rpc_y * rpc_x - fe_0 * rpa_z * rpb_z * rpc_z * rpc_x - fe_0 * rpa_z * rpb_y * rpc_y * rpc_x - fe_0 * rpa_z * rpb_x * rpc_z * rpc_z - fe_0 * rpa_z * rpb_x * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpc_z * rpc_z * rpc_x - fe_0 * rpa_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x - (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y - rpa_z * rpa_z * rpa_y * rpc_z * rpc_y * rpc_x - rpa_z * rpa_z * rpb_z * rpc_y * rpc_y * rpc_x - rpa_z * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x - rpa_z * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x - rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - rpa_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z - rpb_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] * (fe_0 * rpa_z * rpc_z * rpc_z * rpc_x + fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 2.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b4_vals[i] * (rpa_z * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x + rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x + rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y + rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpb_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y + rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x);

        fints_xyz[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_xzz[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_x - 2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x - 3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z - 2.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_y - 3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - 2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_x * rpc_z);

        fints_xzz[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z - rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpc_x - rpa_z * rpa_z * rpb_z * rpb_z * rpb_x * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_x + 3.0 * fe_0 * rpa_z * rpa_y * rpb_x * rpc_z + 3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x + 2.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y + 2.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpb_z * rpb_x * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_x + 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z + 3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_y);

        fints_xzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (4.0 * rpa_z * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z + 2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_x + 2.0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_y + 2.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpc_z * rpc_x + rpa_z * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y + rpa_z * rpa_z * rpb_z * rpb_z * rpc_y * rpc_x + rpa_y * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_x - 2.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_x - 3.0 * fe_0 * rpa_y * rpb_x * rpc_z * rpc_z - 3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_x - 3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-4.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_x - 2.0 * rpa_z * rpa_y * rpb_x * rpc_z * rpc_z * rpc_z - 4.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_x - rpa_z * rpa_z * rpa_y * rpc_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x - rpa_z * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z - rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x - rpb_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 2.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_x + 4.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + 2.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y + rpa_z * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (2.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + rpa_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z + 2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y + rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - 2.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - 2.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_yyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_yyy[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z - 3.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpb_y * rpc_z - 3.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpc_y - rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z + 3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + 9.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y + 3.0 * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_y + 3.0 * fe_0 * rpa_z * rpa_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_y + 3.0 * fe_0 * fe_0 * rpa_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y + (9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 6.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_y + 2.0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_y + 3.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * (3.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpc_y * rpc_y + rpa_y * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y - 9.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y - 6.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpa_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_yyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - 3.0 * fe_0 * fe_0 * rpc_y * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y * rpc_y - 6.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-rpa_z * rpa_z * rpa_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_z * rpa_z * rpb_y * rpc_y * rpc_y * rpc_y - 3.0 * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - rpb_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] * (6.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y + 2.0 * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] * (6.0 * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y * rpc_y + rpa_z * rpa_z * rpc_y * rpc_y * rpc_y * rpc_y + 3.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyy[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_y * rpc_y * rpc_y * rpc_y - 2.0 * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y - rpa_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y * rpc_y;

        fints_yyz[i] += fss * b0_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z + fe_0 * rpa_z * rpa_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y);

        fints_yyz[i] += fss * b0_vals[i] * (fe_0 * fe_0 * rpa_z * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y + rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_z * rpa_y * rpb_z * rpc_z - 2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y - fe_0 * rpa_z * rpa_y * rpb_y * rpb_y - 2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z - fe_0 * rpa_z * rpb_y * rpb_y * rpc_y);

        fints_yyz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_z - fe_0 * rpa_z * rpa_z * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_y - fe_0 * rpa_z * rpa_z * rpb_y * rpc_z);

        fints_yyz[i] += fss * b1_vals[i] * (-fe_0 * rpa_y * rpb_z * rpb_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_y * rpb_y - (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y - fe_0 * fe_0 * rpa_z * rpa_y);

        fints_yyz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * fe_0 * rpa_z * rpb_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y - (1.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z - 2.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpb_y * rpc_z - 2.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpc_y - rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b1_vals[i] * (-rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (fe_0 * rpa_z * rpa_y * rpb_z * rpc_z + 2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_y + fe_0 * rpa_z * rpa_y * rpc_z * rpc_z + fe_0 * rpa_z * rpa_y * rpc_y * rpc_y + 2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y + 2.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z + 2.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y + fe_0 * rpa_z * rpb_y * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_y + fe_0 * rpa_z * rpa_z * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y + fe_0 * rpa_y * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpb_y * rpb_y * rpc_z + fe_0 * rpb_z * rpb_y * rpc_z * rpc_z + fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y + fe_0 * fe_0 * rpa_z * rpb_y + 3.0 * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z + (1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y + 3.0 * fe_0 * fe_0 * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 4.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpb_y * rpb_y * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_y + rpa_z * rpa_z * rpa_y * rpb_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_z * rpa_y * rpb_y * rpc_z * rpc_y + 2.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpc_y * rpc_y + rpa_z * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y + rpa_y * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpa_y * rpc_z * rpc_z - fe_0 * rpa_z * rpa_y * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y - 2.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z - 2.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y - fe_0 * rpa_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z - fe_0 * rpb_z * rpb_y * rpc_z * rpc_z - fe_0 * rpb_z * rpb_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y - 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - fe_0 * rpb_y * rpc_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y * rpc_y - 4.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_y - 4.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y * rpc_y - 2.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - rpa_z * rpa_z * rpa_y * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-rpa_z * rpa_z * rpb_z * rpc_y * rpc_y * rpc_y - 2.0 * rpa_z * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y - 2.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - rpa_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z - rpb_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y + fe_0 * rpa_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpc_y * rpc_y * rpc_y + 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + fe_0 * rpb_y * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 2.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y * rpc_y + 4.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + rpa_z * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (rpa_y * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y + 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - rpa_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - rpb_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyz[i] += fss * b5_vals[i] * (-2.0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_yzz[i] += fss * b0_vals[i] * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y + fe_0 * fe_0 * rpa_z * rpb_z);

        fints_yzz[i] += fss * b0_vals[i] * ((1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y - 2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y - 3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z - 2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y - fe_0 * rpa_z * rpb_z * rpb_z * rpc_z);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_y - fe_0 * rpa_z * rpa_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_y);

        fints_yzz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y - 2.0 * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_yzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z - (1.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_yzz[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - (9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_y * rpc_z - 2.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z);

        fints_yzz[i] += fss * b1_vals[i] * (-rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpc_y - rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_y + 3.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpc_z + 3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y + 2.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y + 2.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y + fe_0 * rpa_z * rpb_z * rpb_z * rpc_z + 3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_y + fe_0 * rpa_z * rpa_z * rpb_z * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_z + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_y + 3.0 * fe_0 * rpa_y * rpb_z * rpb_y * rpc_z + 3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_y + 3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y + fe_0 * fe_0 * rpa_z * rpb_z + 3.0 * fe_0 * fe_0 * rpa_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_y);

        fints_yzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_y + 3.0 * fe_0 * fe_0 * rpb_z * rpc_z + (1.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 4.0 * rpa_z * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z + 2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_y + 2.0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpc_z * rpc_y + rpa_z * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z + 2.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y + rpa_z * rpa_z * rpb_z * rpb_z * rpc_y * rpc_y + rpa_y * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_y - 2.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z - 2.0 * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-fe_0 * rpa_z * rpc_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_y - 3.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_y - 3.0 * fe_0 * rpa_y * rpb_y * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - fe_0 * rpb_z * rpc_z * rpc_z * rpc_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_y * rpc_y - 3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_y - 3.0 * fe_0 * fe_0 * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 4.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-2.0 * rpa_z * rpa_y * rpb_y * rpc_z * rpc_z * rpc_z - 4.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y * rpc_y - rpa_z * rpa_z * rpa_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-rpa_z * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y - 2.0 * rpa_y * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z - rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + fe_0 * rpa_z * rpc_z * rpc_z * rpc_z + 3.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y + fe_0 * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_y + 4.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y + rpa_z * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (rpa_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z + 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y + rpb_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z - 2.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y - 2.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yzz[i] += fss * b5_vals[i] * (-rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y;

        fints_zzz[i] += fss * b0_vals[i] * (3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z + (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z);

        fints_zzz[i] += fss * b0_vals[i] * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z;

        fints_zzz[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z - 3.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z - 3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z - (1.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpb_z - (1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * fe_0 * rpa_z * rpa_y);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y - (9.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpb_z - (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - 2.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpb_z * rpc_z);

        fints_zzz[i] += fss * b1_vals[i] * (-3.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z - rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpc_z + 6.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z + 9.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y + 3.0 * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_y * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y + 9.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_y * rpb_z * rpb_z * rpc_z + (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_y + 3.0 * fe_0 * fe_0 * rpa_z * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpb_z + (15.0 / 2.0) * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 6.0 * rpa_z * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z + 2.0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_y + 3.0 * rpa_z * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * (3.0 * rpa_z * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y + rpa_y * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpa_z * rpa_y * rpc_z * rpc_z - 9.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y - 6.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_y * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-5.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z - 9.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpa_y * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_y - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 6.0 * rpa_z * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z - 6.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - rpa_z * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y - 3.0 * rpa_y * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b4_vals[i] * (6.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y + 5.0 * fe_0 * rpa_y * rpc_z * rpc_z * rpc_z + 9.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y + 5.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_zzz[i] += fss * b4_vals[i] * (2.0 * rpa_z * rpa_y * rpc_z * rpc_z * rpc_z * rpc_z + 6.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + rpa_z * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * rpa_y * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z + 3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - 2.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y - rpa_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z - 3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_zzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y;

    }
}

auto
compPrimitiveNuclearPotentialFF_ZZZ_T(      TDoubleArray& buffer_xxx,
                                            TDoubleArray& buffer_xxy,
                                            TDoubleArray& buffer_xxz,
                                            TDoubleArray& buffer_xyy,
                                            TDoubleArray& buffer_xyz,
                                            TDoubleArray& buffer_xzz,
                                            TDoubleArray& buffer_yyy,
                                            TDoubleArray& buffer_yyz,
                                            TDoubleArray& buffer_yzz,
                                            TDoubleArray& buffer_zzz,
                       const double charge,
                       const TPoint3D& point,
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

    // set up coordinates for C center

    const auto c_rx = point[0];

    const auto c_ry = point[1];

    const auto c_rz = point[2];

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

    // set up Boys function variables

    const CBoysFunc<6> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<7> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

    auto b4_vals = bf_values[4].data();

    auto b5_vals = bf_values[5].data();

    auto b6_vals = bf_values[6].data();

    auto targs = bf_args.data();

    // compute Boys function values

    #pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);
    }

    bf_table.compute<7>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_xxx,\
                             fints_xxy,\
                             fints_xxz,\
                             fints_xyy,\
                             fints_xyz,\
                             fints_xzz,\
                             fints_yyy,\
                             fints_yyz,\
                             fints_yzz,\
                             fints_zzz,\
                             ket_fe,\
                             ket_fn,\
                             ket_rx,\
                             ket_ry,\
                             ket_rz : 64)
    for (int64_t i = 0; i < ket_dim; i++)
    {
        const auto ab_x = bra_rx - ket_rx[i];

        const auto ab_y = bra_ry - ket_ry[i];

        const auto ab_z = bra_rz - ket_rz[i];

        const auto fxi_0 = bra_exp + ket_fe[i];

        const auto fe_0 = 1.0 / fxi_0;

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        fints_xxx[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpb_x - (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_x);

        fints_xxx[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - 3.0 * rpa_z * rpa_z * rpb_x * rpb_x * rpb_x * rpc_z);

        fints_xxx[i] += fss * b1_vals[i] * (-3.0 * rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpb_x * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x);

        fints_xxx[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 3.0 * rpa_z * rpb_x * rpb_x * rpb_x * rpc_z * rpc_z + 9.0 * rpa_z * rpa_z * rpb_x * rpb_x * rpc_z * rpc_x + 3.0 * rpa_z * rpa_z * rpa_z * rpb_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x);

        fints_xxx[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (9.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xxx[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 9.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_x - 9.0 * rpa_z * rpa_z * rpb_x * rpc_z * rpc_x * rpc_x - rpa_z * rpa_z * rpa_z * rpc_x * rpc_x * rpc_x - rpb_x * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xxx[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpc_x * rpc_x * rpc_x + (9.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x);

        fints_xxx[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 9.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_x * rpc_x + 3.0 * rpa_z * rpa_z * rpc_z * rpc_x * rpc_x * rpc_x + 3.0 * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xxx[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_x * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - 3.0 * rpa_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x - 3.0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxx[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x * rpc_x;

        fints_xxy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x);

        fints_xxy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y);

        fints_xxy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z);

        fints_xxy[i] += fss * b1_vals[i] * (-3.0 * rpa_z * rpa_z * rpb_y * rpb_x * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * rpc_x - rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x + 3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_y + 3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpb_x * rpb_x * rpc_z);

        fints_xxy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_xxy[i] += fss * b2_vals[i] * (3.0 * rpa_z * rpb_y * rpb_x * rpb_x * rpc_z * rpc_z + 6.0 * rpa_z * rpa_z * rpb_y * rpb_x * rpc_z * rpc_x + 3.0 * rpa_z * rpa_z * rpb_x * rpb_x * rpc_z * rpc_y + rpa_z * rpa_z * rpa_z * rpb_y * rpc_x * rpc_x + 2.0 * rpa_z * rpa_z * rpa_z * rpb_x * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_x * rpc_x - 3.0 * fe_0 * rpa_z * rpb_x * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y - 3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z - 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_x * rpb_x * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 6.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxy[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_y - 3.0 * rpa_z * rpa_z * rpb_y * rpc_z * rpc_x * rpc_x - 6.0 * rpa_z * rpa_z * rpb_x * rpc_z * rpc_y * rpc_x - rpa_z * rpa_z * rpa_z * rpc_y * rpc_x * rpc_x - rpb_y * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xxy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z + 3.0 * fe_0 * rpb_x * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 3.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_x * rpc_x + 6.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b4_vals[i] * (3.0 * rpa_z * rpa_z * rpc_z * rpc_y * rpc_x * rpc_x + 2.0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xxy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - 3.0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x - rpb_y * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xxy[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x * rpc_x;

        fints_xxz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z);

        fints_xxz[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpb_x - (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_z - 3.0 * fe_0 * rpa_z * rpa_z * rpb_x * rpc_x);

        fints_xxz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_xxz[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpb_x);

        fints_xxz[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 3.0 * rpa_z * rpa_z * rpb_z * rpb_x * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * rpc_x - rpa_z * rpa_z * rpa_z * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x + 9.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpb_x * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_z + 3.0 * fe_0 * rpa_z * rpa_z * rpb_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpb_x * rpc_z + 3.0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_x * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 3.0 * rpa_z * rpb_z * rpb_x * rpb_x * rpc_z * rpc_z + 6.0 * rpa_z * rpa_z * rpb_z * rpb_x * rpc_z * rpc_x + 3.0 * rpa_z * rpa_z * rpb_x * rpb_x * rpc_z * rpc_z);

        fints_xxz[i] += fss * b2_vals[i] * (rpa_z * rpa_z * rpa_z * rpb_z * rpc_x * rpc_x + 2.0 * rpa_z * rpa_z * rpa_z * rpb_x * rpc_z * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_x * rpc_x - 9.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_x * rpc_x - 3.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x - 3.0 * fe_0 * rpb_x * rpb_x * rpc_z * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_x);

        fints_xxz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * fe_0 * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_x * rpc_x - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_x - 3.0 * rpa_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xxz[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_x * rpc_x - 6.0 * rpa_z * rpa_z * rpb_x * rpc_z * rpc_z * rpc_x - rpa_z * rpa_z * rpa_z * rpc_z * rpc_x * rpc_x - rpb_z * rpb_x * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xxz[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_x * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z + 6.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_x * rpc_x + 3.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_x * rpc_x);

        fints_xxz[i] += fss * b4_vals[i] * (6.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x + 3.0 * rpa_z * rpa_z * rpc_z * rpc_z * rpc_x * rpc_x + 2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_x + rpb_x * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_xxz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_x * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z - 3.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - rpb_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x - 2.0 * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xxz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x * rpc_x;

        fints_xyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x);

        fints_xyy[i] += fss * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z);

        fints_xyy[i] += fss * b1_vals[i] * (-3.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_x * rpc_z - 2.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * rpc_y - rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_y * rpb_x * rpc_y + 3.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_x + 3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_x * rpc_z);

        fints_xyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xyy[i] += fss * b2_vals[i] * (3.0 * rpa_z * rpb_y * rpb_y * rpb_x * rpc_z * rpc_z + 6.0 * rpa_z * rpa_z * rpb_y * rpb_x * rpc_z * rpc_y + 3.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpc_z * rpc_x + 2.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpc_y * rpc_x + rpa_z * rpa_z * rpa_z * rpb_x * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpa_z * rpb_y * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x);

        fints_xyy[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x - 3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_y - 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-(1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 6.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyy[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_x - 6.0 * rpa_z * rpa_z * rpb_y * rpc_z * rpc_y * rpc_x - 3.0 * rpa_z * rpa_z * rpb_x * rpc_z * rpc_y * rpc_y - rpa_z * rpa_z * rpa_z * rpc_y * rpc_y * rpc_x - rpb_y * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_x + 3.0 * fe_0 * rpb_y * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x + (3.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 6.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_x + 3.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyy[i] += fss * b4_vals[i] * (3.0 * rpa_z * rpa_z * rpc_z * rpc_y * rpc_y * rpc_x + 2.0 * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y + rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xyy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_x - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - 3.0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x - 2.0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_xyy[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_x;

        fints_xyz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x + rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_xyz[i] += fss * b1_vals[i] * (-(3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - 3.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_x * rpc_z - rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpc_x - rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * rpc_y - rpa_z * rpa_z * rpa_z * rpb_y * rpb_x * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x + (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_x * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_x * rpc_z);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y + 3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_x);

        fints_xyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 3.0 * rpa_z * rpb_z * rpb_y * rpb_x * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpc_z * rpc_x + 3.0 * rpa_z * rpa_z * rpb_z * rpb_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b2_vals[i] * (3.0 * rpa_z * rpa_z * rpb_y * rpb_x * rpc_z * rpc_z + rpa_z * rpa_z * rpa_z * rpb_z * rpc_y * rpc_x + rpa_z * rpa_z * rpa_z * rpb_y * rpc_z * rpc_x + rpa_z * rpa_z * rpa_z * rpb_x * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpb_x * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_x * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x - 3.0 * fe_0 * rpb_y * rpb_x * rpc_z * rpc_z - 3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y - (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_x - (3.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_x - 3.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_x);

        fints_xyz[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_y - 3.0 * rpa_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z - 3.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_y * rpc_x - 3.0 * rpa_z * rpa_z * rpb_y * rpc_z * rpc_z * rpc_x - 3.0 * rpa_z * rpa_z * rpb_x * rpc_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b3_vals[i] * (-rpa_z * rpa_z * rpa_z * rpc_z * rpc_y * rpc_x - rpb_z * rpb_y * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xyz[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_x + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_x + 3.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_x + 3.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_y + 3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_x + 3.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_x + 3.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x + 3.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * rpa_z * rpa_z * rpc_z * rpc_z * rpc_y * rpc_x);

        fints_xyz[i] += fss * b4_vals[i] * (rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_x + rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_y + rpb_y * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_xyz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_x - 3.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x - rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_xyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_x;

        fints_xzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x + 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x);

        fints_xzz[i] += fss * b0_vals[i] * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x;

        fints_xzz[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x - 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_x - 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpc_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_x);

        fints_xzz[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - 3.0 * fe_0 * fe_0 * rpb_z * rpb_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x - (15.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - 3.0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_x * rpc_z);

        fints_xzz[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_x * rpc_z - rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_z * rpb_z * rpb_x * rpc_z + 9.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_x + 9.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z + 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_x * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_x + 6.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_x * rpc_z);

        fints_xzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_x + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_x + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_x + 3.0 * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * ((15.0 / 2.0) * fe_0 * fe_0 * rpb_x * rpc_z + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x + 3.0 * rpa_z * rpb_z * rpb_z * rpb_x * rpc_z * rpc_z + 6.0 * rpa_z * rpa_z * rpb_z * rpb_x * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpb_z * rpb_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpc_z * rpc_x + rpa_z * rpa_z * rpa_z * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_x - 9.0 * fe_0 * rpa_z * rpb_x * rpc_z * rpc_z - 9.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x - (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_x - 6.0 * fe_0 * rpb_z * rpb_x * rpc_z * rpc_z);

        fints_xzz[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_x - 5.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_x - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-(15.0 / 4.0) * fe_0 * fe_0 * rpb_x * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_x - 6.0 * rpa_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z - 3.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_x - 6.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpa_z * rpb_x * rpc_z * rpc_z * rpc_z - rpa_z * rpa_z * rpa_z * rpc_z * rpc_z * rpc_x - rpb_z * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z);

        fints_xzz[i] += fss * b4_vals[i] * (9.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_x + 6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_x + 5.0 * fe_0 * rpb_x * rpc_z * rpc_z * rpc_z + 5.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_x);

        fints_xzz[i] += fss * b4_vals[i] * (6.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x + 3.0 * rpa_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpc_z * rpc_z * rpc_z * rpc_x + 2.0 * rpb_z * rpb_x * rpc_z * rpc_z * rpc_z * rpc_z + rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_x);

        fints_xzz[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_x - 3.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - 2.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x - rpb_x * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_xzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_x;

        fints_yyy[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpb_y - (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_y);

        fints_yyy[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - 3.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpb_y * rpc_z);

        fints_yyy[i] += fss * b1_vals[i] * (-3.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_y + (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpb_y * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y);

        fints_yyy[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 3.0 * rpa_z * rpb_y * rpb_y * rpb_y * rpc_z * rpc_z + 9.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpc_z * rpc_y + 3.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y);

        fints_yyy[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z - (9.0 / 2.0) * fe_0 * rpb_y * rpb_y * rpc_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (9.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z);

        fints_yyy[i] += fss * b3_vals[i] * (-(9.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 9.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_y - 9.0 * rpa_z * rpa_z * rpb_y * rpc_z * rpc_y * rpc_y - rpa_z * rpa_z * rpa_z * rpc_y * rpc_y * rpc_y - rpb_y * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yyy[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpc_y * rpc_y * rpc_y + (9.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y);

        fints_yyy[i] += fss * b4_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 9.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y * rpc_y + 3.0 * rpa_z * rpa_z * rpc_z * rpc_y * rpc_y * rpc_y + 3.0 * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yyy[i] += fss * b5_vals[i] * (-(3.0 / 2.0) * fe_0 * rpc_z * rpc_y * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - 3.0 * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y - 3.0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyy[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y * rpc_y;

        fints_yyz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z);

        fints_yyz[i] += fss * b0_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_y * rpb_y - (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_z - 3.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpc_y);

        fints_yyz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z);

        fints_yyz[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpb_y);

        fints_yyz[i] += fss * b1_vals[i] * (-(9.0 / 8.0) * fe_0 * fe_0 * fe_0 - 3.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpb_y * rpc_z - 2.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpc_y - rpa_z * rpa_z * rpa_z * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y + 9.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y + (9.0 / 2.0) * fe_0 * rpa_z * rpb_y * rpb_y * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_z + 3.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * (3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpb_z * rpb_y * rpb_y * rpc_z + 3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z + 3.0 * fe_0 * fe_0 * rpb_y * rpc_y + (3.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * ((3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + (9.0 / 8.0) * fe_0 * fe_0 * fe_0 + 3.0 * rpa_z * rpb_z * rpb_y * rpb_y * rpc_z * rpc_z + 6.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpc_z * rpc_y + 3.0 * rpa_z * rpa_z * rpb_y * rpb_y * rpc_z * rpc_z);

        fints_yyz[i] += fss * b2_vals[i] * (rpa_z * rpa_z * rpa_z * rpb_z * rpc_y * rpc_y + 2.0 * rpa_z * rpa_z * rpa_z * rpb_y * rpc_z * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpc_y * rpc_y - 9.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y - (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_y * rpc_y - 3.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y - 3.0 * fe_0 * rpb_y * rpb_y * rpc_z * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (3.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_y);

        fints_yyz[i] += fss * b3_vals[i] * (-3.0 * fe_0 * fe_0 * rpc_z * rpc_z - (3.0 / 2.0) * fe_0 * fe_0 * rpc_y * rpc_y - (3.0 / 8.0) * fe_0 * fe_0 * fe_0 - 6.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_y - 3.0 * rpa_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yyz[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_y * rpc_y - 6.0 * rpa_z * rpa_z * rpb_y * rpc_z * rpc_z * rpc_y - rpa_z * rpa_z * rpa_z * rpc_z * rpc_y * rpc_y - rpb_z * rpb_y * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yyz[i] += fss * b4_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_y * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z + 6.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y + (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z + (3.0 / 4.0) * fe_0 * fe_0 * rpc_y * rpc_y + 3.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y * rpc_y);

        fints_yyz[i] += fss * b4_vals[i] * (6.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * rpa_z * rpa_z * rpc_z * rpc_z * rpc_y * rpc_y + 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_y + rpb_y * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_yyz[i] += fss * b5_vals[i] * (-3.0 * fe_0 * rpc_z * rpc_z * rpc_y * rpc_y - (1.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z - 3.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - rpb_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y - 2.0 * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yyz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y * rpc_y;

        fints_yzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y + 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y);

        fints_yzz[i] += fss * b0_vals[i] * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y;

        fints_yzz[i] += fss * b1_vals[i] * (-9.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y - 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y - 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpc_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_z - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y - (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_y);

        fints_yzz[i] += fss * b1_vals[i] * (-(9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - 3.0 * fe_0 * fe_0 * rpb_z * rpb_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y - (15.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - 3.0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_yzz[i] += fss * b1_vals[i] * (-2.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_y * rpc_z - rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpa_z * rpb_z * rpb_y * rpc_z + 9.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y + (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_y + 9.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z + 3.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * ((9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_y * rpc_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y + (1.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_y + 6.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_y * rpc_z);

        fints_yzz[i] += fss * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_y + (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_y + (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_y + 3.0 * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * ((15.0 / 2.0) * fe_0 * fe_0 * rpb_y * rpc_z + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y + 3.0 * rpa_z * rpb_z * rpb_z * rpb_y * rpc_z * rpc_z + 6.0 * rpa_z * rpa_z * rpb_z * rpb_y * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpb_z * rpb_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b2_vals[i] * (2.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpc_z * rpc_y + rpa_z * rpa_z * rpa_z * rpb_y * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-9.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_y - 9.0 * fe_0 * rpa_z * rpb_y * rpc_z * rpc_z - 9.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y - (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpc_z * rpc_y - 6.0 * fe_0 * rpb_z * rpb_y * rpc_z * rpc_z);

        fints_yzz[i] += fss * b3_vals[i] * (-6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpc_z * rpc_y - 5.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z - (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_y - (3.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-(15.0 / 4.0) * fe_0 * fe_0 * rpb_y * rpc_z - (15.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_y - 6.0 * rpa_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z - 3.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_y - 6.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b3_vals[i] * (-3.0 * rpa_z * rpa_z * rpb_y * rpc_z * rpc_z * rpc_z - rpa_z * rpa_z * rpa_z * rpc_z * rpc_z * rpc_y - rpb_z * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b4_vals[i] * (9.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_y + 6.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_y + 5.0 * fe_0 * rpb_y * rpc_z * rpc_z * rpc_z + 5.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y + (15.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_y);

        fints_yzz[i] += fss * b4_vals[i] * (6.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y + 3.0 * rpa_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpc_z * rpc_z * rpc_z * rpc_y + 2.0 * rpb_z * rpb_y * rpc_z * rpc_z * rpc_z * rpc_z + rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_y);

        fints_yzz[i] += fss * b5_vals[i] * (-5.0 * fe_0 * rpc_z * rpc_z * rpc_z * rpc_y - 3.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y - 2.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y - rpb_y * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_yzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_y;

        fints_zzz[i] += fss * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z + (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z + (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z);

        fints_zzz[i] += fss * b0_vals[i] * ((9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (15.0 / 8.0) * fe_0 * fe_0 * fe_0 + rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(27.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z - (3.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpb_z - (27.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z - (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_z - (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z - (27.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpb_z - (45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpa_z);

        fints_zzz[i] += fss * b1_vals[i] * (-(45.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (9.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpb_z - (45.0 / 8.0) * fe_0 * fe_0 * fe_0 - 3.0 * rpa_z * rpa_z * rpb_z * rpb_z * rpb_z * rpc_z - 3.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpb_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * (27.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z + (27.0 / 2.0) * fe_0 * rpa_z * rpb_z * rpb_z * rpc_z + (27.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpb_z * rpc_z + 9.0 * fe_0 * rpa_z * rpa_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpa_z * rpa_z * rpa_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * (9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z + (3.0 / 2.0) * fe_0 * rpb_z * rpb_z * rpb_z * rpc_z + (27.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpb_z + (45.0 / 2.0) * fe_0 * fe_0 * rpa_z * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpa_z);

        fints_zzz[i] += fss * b2_vals[i] * ((45.0 / 2.0) * fe_0 * fe_0 * rpb_z * rpc_z + (9.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpb_z + (45.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + (45.0 / 8.0) * fe_0 * fe_0 * fe_0 + 3.0 * rpa_z * rpb_z * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b2_vals[i] * (9.0 * rpa_z * rpa_z * rpb_z * rpb_z * rpc_z * rpc_z + 3.0 * rpa_z * rpa_z * rpa_z * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-27.0 * fe_0 * rpa_z * rpb_z * rpc_z * rpc_z - 15.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z - 9.0 * fe_0 * rpa_z * rpa_z * rpc_z * rpc_z - 15.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z - 9.0 * fe_0 * rpb_z * rpb_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-(45.0 / 4.0) * fe_0 * fe_0 * rpa_z * rpc_z - (45.0 / 4.0) * fe_0 * fe_0 * rpb_z * rpc_z - (45.0 / 2.0) * fe_0 * fe_0 * rpc_z * rpc_z - (15.0 / 8.0) * fe_0 * fe_0 * fe_0 - 9.0 * rpa_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b3_vals[i] * (-9.0 * rpa_z * rpa_z * rpb_z * rpc_z * rpc_z * rpc_z - rpa_z * rpa_z * rpa_z * rpc_z * rpc_z * rpc_z - rpb_z * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b4_vals[i] * (15.0 * fe_0 * rpa_z * rpc_z * rpc_z * rpc_z + 15.0 * fe_0 * rpb_z * rpc_z * rpc_z * rpc_z + (15.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z + (45.0 / 4.0) * fe_0 * fe_0 * rpc_z * rpc_z + 9.0 * rpa_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b4_vals[i] * (3.0 * rpa_z * rpa_z * rpc_z * rpc_z * rpc_z * rpc_z + 3.0 * rpb_z * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b5_vals[i] * (-(15.0 / 2.0) * fe_0 * rpc_z * rpc_z * rpc_z * rpc_z - 3.0 * rpa_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z - 3.0 * rpb_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z);

        fints_zzz[i] += fss * b6_vals[i] * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z * rpc_z;

    }
}

} // npotrec namespace

