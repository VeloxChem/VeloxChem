#include "DipoleRecFF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveDipoleFF_XXX_XXX.hpp"
#include "PrimitiveDipoleFF_XXX_XXY.hpp"
#include "PrimitiveDipoleFF_XXX_XXZ.hpp"
#include "PrimitiveDipoleFF_XXX_XYY.hpp"
#include "PrimitiveDipoleFF_XXX_XYZ.hpp"
#include "PrimitiveDipoleFF_XXX_XZZ.hpp"
#include "PrimitiveDipoleFF_XXX_YYY.hpp"
#include "PrimitiveDipoleFF_XXX_YYZ.hpp"
#include "PrimitiveDipoleFF_XXX_YZZ.hpp"
#include "PrimitiveDipoleFF_XXX_ZZZ.hpp"
#include "PrimitiveDipoleFF_XXY_XXX.hpp"
#include "PrimitiveDipoleFF_XXY_XXY.hpp"
#include "PrimitiveDipoleFF_XXY_XXZ.hpp"
#include "PrimitiveDipoleFF_XXY_XYY.hpp"
#include "PrimitiveDipoleFF_XXY_XYZ.hpp"
#include "PrimitiveDipoleFF_XXY_XZZ.hpp"
#include "PrimitiveDipoleFF_XXY_YYY.hpp"
#include "PrimitiveDipoleFF_XXY_YYZ.hpp"
#include "PrimitiveDipoleFF_XXY_YZZ.hpp"
#include "PrimitiveDipoleFF_XXY_ZZZ.hpp"
#include "PrimitiveDipoleFF_XXZ_XXX.hpp"
#include "PrimitiveDipoleFF_XXZ_XXY.hpp"
#include "PrimitiveDipoleFF_XXZ_XXZ.hpp"
#include "PrimitiveDipoleFF_XXZ_XYY.hpp"
#include "PrimitiveDipoleFF_XXZ_XYZ.hpp"
#include "PrimitiveDipoleFF_XXZ_XZZ.hpp"
#include "PrimitiveDipoleFF_XXZ_YYY.hpp"
#include "PrimitiveDipoleFF_XXZ_YYZ.hpp"
#include "PrimitiveDipoleFF_XXZ_YZZ.hpp"
#include "PrimitiveDipoleFF_XXZ_ZZZ.hpp"
#include "PrimitiveDipoleFF_XYY_XXX.hpp"
#include "PrimitiveDipoleFF_XYY_XXY.hpp"
#include "PrimitiveDipoleFF_XYY_XXZ.hpp"
#include "PrimitiveDipoleFF_XYY_XYY.hpp"
#include "PrimitiveDipoleFF_XYY_XYZ.hpp"
#include "PrimitiveDipoleFF_XYY_XZZ.hpp"
#include "PrimitiveDipoleFF_XYY_YYY.hpp"
#include "PrimitiveDipoleFF_XYY_YYZ.hpp"
#include "PrimitiveDipoleFF_XYY_YZZ.hpp"
#include "PrimitiveDipoleFF_XYY_ZZZ.hpp"
#include "PrimitiveDipoleFF_XYZ_XXX.hpp"
#include "PrimitiveDipoleFF_XYZ_XXY.hpp"
#include "PrimitiveDipoleFF_XYZ_XXZ.hpp"
#include "PrimitiveDipoleFF_XYZ_XYY.hpp"
#include "PrimitiveDipoleFF_XYZ_XYZ.hpp"
#include "PrimitiveDipoleFF_XYZ_XZZ.hpp"
#include "PrimitiveDipoleFF_XYZ_YYY.hpp"
#include "PrimitiveDipoleFF_XYZ_YYZ.hpp"
#include "PrimitiveDipoleFF_XYZ_YZZ.hpp"
#include "PrimitiveDipoleFF_XYZ_ZZZ.hpp"
#include "PrimitiveDipoleFF_XZZ_XXX.hpp"
#include "PrimitiveDipoleFF_XZZ_XXY.hpp"
#include "PrimitiveDipoleFF_XZZ_XXZ.hpp"
#include "PrimitiveDipoleFF_XZZ_XYY.hpp"
#include "PrimitiveDipoleFF_XZZ_XYZ.hpp"
#include "PrimitiveDipoleFF_XZZ_XZZ.hpp"
#include "PrimitiveDipoleFF_XZZ_YYY.hpp"
#include "PrimitiveDipoleFF_XZZ_YYZ.hpp"
#include "PrimitiveDipoleFF_XZZ_YZZ.hpp"
#include "PrimitiveDipoleFF_XZZ_ZZZ.hpp"
#include "PrimitiveDipoleFF_YYY_XXX.hpp"
#include "PrimitiveDipoleFF_YYY_XXY.hpp"
#include "PrimitiveDipoleFF_YYY_XXZ.hpp"
#include "PrimitiveDipoleFF_YYY_XYY.hpp"
#include "PrimitiveDipoleFF_YYY_XYZ.hpp"
#include "PrimitiveDipoleFF_YYY_XZZ.hpp"
#include "PrimitiveDipoleFF_YYY_YYY.hpp"
#include "PrimitiveDipoleFF_YYY_YYZ.hpp"
#include "PrimitiveDipoleFF_YYY_YZZ.hpp"
#include "PrimitiveDipoleFF_YYY_ZZZ.hpp"
#include "PrimitiveDipoleFF_YYZ_XXX.hpp"
#include "PrimitiveDipoleFF_YYZ_XXY.hpp"
#include "PrimitiveDipoleFF_YYZ_XXZ.hpp"
#include "PrimitiveDipoleFF_YYZ_XYY.hpp"
#include "PrimitiveDipoleFF_YYZ_XYZ.hpp"
#include "PrimitiveDipoleFF_YYZ_XZZ.hpp"
#include "PrimitiveDipoleFF_YYZ_YYY.hpp"
#include "PrimitiveDipoleFF_YYZ_YYZ.hpp"
#include "PrimitiveDipoleFF_YYZ_YZZ.hpp"
#include "PrimitiveDipoleFF_YYZ_ZZZ.hpp"
#include "PrimitiveDipoleFF_YZZ_XXX.hpp"
#include "PrimitiveDipoleFF_YZZ_XXY.hpp"
#include "PrimitiveDipoleFF_YZZ_XXZ.hpp"
#include "PrimitiveDipoleFF_YZZ_XYY.hpp"
#include "PrimitiveDipoleFF_YZZ_XYZ.hpp"
#include "PrimitiveDipoleFF_YZZ_XZZ.hpp"
#include "PrimitiveDipoleFF_YZZ_YYY.hpp"
#include "PrimitiveDipoleFF_YZZ_YYZ.hpp"
#include "PrimitiveDipoleFF_YZZ_YZZ.hpp"
#include "PrimitiveDipoleFF_YZZ_ZZZ.hpp"
#include "PrimitiveDipoleFF_ZZZ_XXX.hpp"
#include "PrimitiveDipoleFF_ZZZ_XXY.hpp"
#include "PrimitiveDipoleFF_ZZZ_XXZ.hpp"
#include "PrimitiveDipoleFF_ZZZ_XYY.hpp"
#include "PrimitiveDipoleFF_ZZZ_XYZ.hpp"
#include "PrimitiveDipoleFF_ZZZ_XZZ.hpp"
#include "PrimitiveDipoleFF_ZZZ_YYY.hpp"
#include "PrimitiveDipoleFF_ZZZ_YYZ.hpp"
#include "PrimitiveDipoleFF_ZZZ_YZZ.hpp"
#include "PrimitiveDipoleFF_ZZZ_ZZZ.hpp"
#include "T2CDistributor.hpp"

namespace mpol {  // mpol namespace

auto
compDipoleFF(CSubMatrix*      matrix_x,
             CSubMatrix*      matrix_y,
             CSubMatrix*      matrix_z,
             const TPoint3D&  point,
             const CGtoBlock& gto_block,
             const int64_t    bra_first,
             const int64_t    bra_last) -> void
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

    alignas(64) TDoubleArray buffer_x;

    alignas(64) TDoubleArray buffer_y;

    alignas(64) TDoubleArray buffer_z;

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

            // compute primitive integrals block (XXX_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XXX_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XXX_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XXX_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XXX_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_15, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_15, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_15, gto_indexes, 6, 1, j, ket_first, ket_last);

            // compute primitive integrals block (XXX_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 4.0 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 4.0 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 4.0 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XXX_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XXX_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XXX_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 4.0 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 4.0 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 4.0 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XXX_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 2.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 2.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 2.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            // compute primitive integrals block (XXY_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XXY_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XXY_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XXY_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XXY_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_15, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_15, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_15, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            // compute primitive integrals block (XXY_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XXY_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XXY_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XXY_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XXY_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 2.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 2.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 2.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            // compute primitive integrals block (XXZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XXZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XXZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XXZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XXZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_15, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_15, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_15, gto_indexes, 5, 1, j, ket_first, ket_last);

            // compute primitive integrals block (XXZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XXZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XXZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XXZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XXZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 2.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 2.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 2.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            // compute primitive integrals block (XYY_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XYY_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XYY_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XYY_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f3_5, gto_indexes, 6, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XYY_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_15, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_15, gto_indexes, 6, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_15, gto_indexes, 6, 1, j, ket_first, ket_last);

            // compute primitive integrals block (XYY_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 6, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XYY_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_5, gto_indexes, 6, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XYY_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 3.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 0.5 * f3_15, gto_indexes, 6, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XYY_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 4.0 * f3_3, gto_indexes, 6, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XYY_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 2.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 2.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 2.0, gto_indexes, 6, 3, j, ket_first, ket_last);

            // compute primitive integrals block (XYZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * f3_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * f3_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * f3_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XYZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 3.0 * f3_5, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 3.0 * f3_5, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 3.0 * f3_5, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XYZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 0.5 * f3_15, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 0.5 * f3_15, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 0.5 * f3_15, gto_indexes, 1, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XYZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0 * f3_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0 * f3_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0 * f3_5, gto_indexes, 1, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XYZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * f3_15, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * f3_15, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * f3_15, gto_indexes, 1, 1, j, ket_first, ket_last);

            // compute primitive integrals block (XYZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 4.0 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 4.0 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 4.0 * f3_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XYZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f3_5, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f3_5, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f3_5, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XYZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 0.5 * f3_15, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 0.5 * f3_15, gto_indexes, 1, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 0.5 * f3_15, gto_indexes, 1, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XYZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 4.0 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 4.0 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 4.0 * f3_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XYZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 2.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 2.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 2.0, gto_indexes, 1, 3, j, ket_first, ket_last);

            // compute primitive integrals block (XZZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XZZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XZZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XZZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 4, 6, j, ket_first, ket_last);

            // compute primitive integrals block (XZZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f3_15, gto_indexes, 4, 1, j, ket_first, ket_last);

            // compute primitive integrals block (XZZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            // compute primitive integrals block (XZZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_5, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XZZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 3.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 4, 5, j, ket_first, ket_last);

            // compute primitive integrals block (XZZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XZZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 2.0, gto_indexes, 4, 3, j, ket_first, ket_last);

            // compute primitive integrals block (YYY_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            // compute primitive integrals block (YYY_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            // compute primitive integrals block (YYY_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            // compute primitive integrals block (YYY_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f3_5, gto_indexes, 0, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            // compute primitive integrals block (YYY_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_15, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_15, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_15, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            // compute primitive integrals block (YYY_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 4.0 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 4.0 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 4.0 * f3_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            // compute primitive integrals block (YYY_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_5, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            // compute primitive integrals block (YYY_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 0.5 * f3_15, gto_indexes, 0, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            // compute primitive integrals block (YYY_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 4.0 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 4.0 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 4.0 * f3_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            // compute primitive integrals block (YYY_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 2.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 2.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 2.0, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            // compute primitive integrals block (YYZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            // compute primitive integrals block (YYZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            // compute primitive integrals block (YYZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            // compute primitive integrals block (YYZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f3_5, gto_indexes, 5, 6, j, ket_first, ket_last);

            // compute primitive integrals block (YYZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_15, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_15, gto_indexes, 5, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_15, gto_indexes, 5, 1, j, ket_first, ket_last);

            // compute primitive integrals block (YYZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 4, j, ket_first, ket_last);

            // compute primitive integrals block (YYZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_5, gto_indexes, 5, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            // compute primitive integrals block (YYZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 3.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 0.5 * f3_15, gto_indexes, 5, 5, j, ket_first, ket_last);

            // compute primitive integrals block (YYZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 4.0 * f3_3, gto_indexes, 5, 2, j, ket_first, ket_last);

            // compute primitive integrals block (YYZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 2.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 2.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 2.0, gto_indexes, 5, 3, j, ket_first, ket_last);

            // compute primitive integrals block (YZZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            // compute primitive integrals block (YZZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            // compute primitive integrals block (YZZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            // compute primitive integrals block (YZZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f3_5, gto_indexes, 2, 6, j, ket_first, ket_last);

            // compute primitive integrals block (YZZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f3_15, gto_indexes, 2, 1, j, ket_first, ket_last);

            // compute primitive integrals block (YZZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            // compute primitive integrals block (YZZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_5, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            // compute primitive integrals block (YZZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 3.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 0.5 * f3_15, gto_indexes, 2, 5, j, ket_first, ket_last);

            // compute primitive integrals block (YZZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f3_3, gto_indexes, 2, 2, j, ket_first, ket_last);

            // compute primitive integrals block (YZZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 2.0, gto_indexes, 2, 3, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 3.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0 * f3_5, gto_indexes, 3, 6, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * f3_15, gto_indexes, 3, 1, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 4.0 * f3_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f3_5, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 0.5 * f3_15, gto_indexes, 3, 5, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 4.0 * f3_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            // compute primitive integrals block (ZZZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 2.0, gto_indexes, 3, 3, j, ket_first, ket_last);
        }
    }
}

auto
compDipoleFF(CSubMatrix*      matrix_x,
             CSubMatrix*      matrix_y,
             CSubMatrix*      matrix_z,
             const TPoint3D&  point,
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

    alignas(64) TDoubleArray buffer_x;

    alignas(64) TDoubleArray buffer_y;

    alignas(64) TDoubleArray buffer_z;

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

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXX_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXX_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXY_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXY_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XXZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XXZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 6, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 3.0 * f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 6, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -3.0 * f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 6, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYY_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYY_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 6, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 1, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XYZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XYZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 4, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_XZZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_5, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 3.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 0, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYY_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYY_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_5, bra_gto_indexes, ket_gto_indexes, 5, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15 * 3.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 0.5 * f3_15 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 5, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -0.5 * f3_15 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 5, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YYZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YYZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15 * 2.0, bra_gto_indexes, ket_gto_indexes, 5, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_5, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_x, buffer_x, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -4.0 * f3_3 * 3.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, -4.0 * f3_3 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 2, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(
                matrix_x, buffer_x, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_y, buffer_y, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix_z, buffer_z, 4.0 * f3_3 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_YZZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XXX)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XXX(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XXY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XXY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XXZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XXZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 6, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_XZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_XZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_YYY)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_YYY(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f3_5, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_YYZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_YYZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_x, buffer_x, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 3.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, -2.0 * 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes, 3, 5, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_YZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_YZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZZ_ZZZ)

            simd::zero(buffer_x);

            simd::zero(buffer_y);

            simd::zero(buffer_z);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    mpol::compPrimitiveDipoleFF_ZZZ_ZZZ(buffer_x,
                                                        buffer_y,
                                                        buffer_z,
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

            t2cfunc::distribute(matrix_x, buffer_x, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);
        }
    }
}

}  // namespace mpol
