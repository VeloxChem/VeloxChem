#include "NuclearPotentialGeom010RecFS.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "MathConst.hpp"
#include "BoysFunc.hpp"
#include "T2CDistributor.hpp"

namespace geom_npotrec { // geom_npotrec namespace

auto
compNuclearPotentialGeom010FS(      CSubMatrix* matrix_x,
                                    CSubMatrix* matrix_y,
                                    CSubMatrix* matrix_z,
                              const TPoint3D& dipole,
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

    alignas(64) TDoubleArray buffer_x;

    alignas(64) TDoubleArray buffer_y;

    alignas(64) TDoubleArray buffer_z;

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

            // compute primitive integrals block (XXX_0)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FS_XXX_0(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
                                                                               point,
                                                                               bra_exp,
                                                                               bra_norm,
                                                                               bra_coord,
                                                                               ket_exps,
                                                                               ket_norms,
                                                                               ket_coords_x,
                                                                               ket_coords_y,
                                                                               ket_coords_z,
                                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXY_0)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FS_XXY_0(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
                                                                               point,
                                                                               bra_exp,
                                                                               bra_norm,
                                                                               bra_coord,
                                                                               ket_exps,
                                                                               ket_norms,
                                                                               ket_coords_x,
                                                                               ket_coords_y,
                                                                               ket_coords_z,
                                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x, buffer_x, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XXZ_0)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FS_XXZ_0(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
                                                                               point,
                                                                               bra_exp,
                                                                               bra_norm,
                                                                               bra_coord,
                                                                               ket_exps,
                                                                               ket_norms,
                                                                               ket_coords_x,
                                                                               ket_coords_y,
                                                                               ket_coords_z,
                                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYY_0)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FS_XYY_0(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
                                                                               point,
                                                                               bra_exp,
                                                                               bra_norm,
                                                                               bra_coord,
                                                                               ket_exps,
                                                                               ket_norms,
                                                                               ket_coords_x,
                                                                               ket_coords_y,
                                                                               ket_coords_z,
                                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                6, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XYZ_0)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FS_XYZ_0(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
                                                                               point,
                                                                               bra_exp,
                                                                               bra_norm,
                                                                               bra_coord,
                                                                               ket_exps,
                                                                               ket_norms,
                                                                               ket_coords_x,
                                                                               ket_coords_y,
                                                                               ket_coords_z,
                                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x, buffer_x, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, f3_15, bra_gto_indexes, ket_gto_indexes,
                                1, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (XZZ_0)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FS_XZZ_0(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
                                                                               point,
                                                                               bra_exp,
                                                                               bra_norm,
                                                                               bra_coord,
                                                                               ket_exps,
                                                                               ket_norms,
                                                                               ket_coords_x,
                                                                               ket_coords_y,
                                                                               ket_coords_z,
                                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                4, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYY_0)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FS_YYY_0(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
                                                                               point,
                                                                               bra_exp,
                                                                               bra_norm,
                                                                               bra_coord,
                                                                               ket_exps,
                                                                               ket_norms,
                                                                               ket_coords_x,
                                                                               ket_coords_y,
                                                                               ket_coords_z,
                                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x, buffer_x, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YYZ_0)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FS_YYZ_0(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
                                                                               point,
                                                                               bra_exp,
                                                                               bra_norm,
                                                                               bra_coord,
                                                                               ket_exps,
                                                                               ket_norms,
                                                                               ket_coords_x,
                                                                               ket_coords_y,
                                                                               ket_coords_z,
                                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x, buffer_x, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_x, buffer_x, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -3.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                5, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (YZZ_0)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FS_YZZ_0(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
                                                                               point,
                                                                               bra_exp,
                                                                               bra_norm,
                                                                               bra_coord,
                                                                               ket_exps,
                                                                               ket_norms,
                                                                               ket_coords_x,
                                                                               ket_coords_y,
                                                                               ket_coords_z,
                                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x, buffer_x, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                2, 0, j, ket_first, ket_last, ang_order);

            // compute primitive integrals block (ZZZ_0)

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

                    geom_npotrec::compPrimitiveNuclearPotentialGeom010FS_ZZZ_0(buffer_x,
                                                                               buffer_y,
                                                                               buffer_z,
                                                                               dipole,
                                                                               point,
                                                                               bra_exp,
                                                                               bra_norm,
                                                                               bra_coord,
                                                                               ket_exps,
                                                                               ket_norms,
                                                                               ket_coords_x,
                                                                               ket_coords_y,
                                                                               ket_coords_z,
                                                                               ket_dim);
                }
            }

            t2cfunc::distribute(matrix_x, buffer_x, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_y, buffer_y, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix_z, buffer_z, 2.0, bra_gto_indexes, ket_gto_indexes,
                                3, 0, j, ket_first, ket_last, ang_order);

        }
    }
}

auto
compPrimitiveNuclearPotentialGeom010FS_XXX_0(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<3> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<4> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

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

    bf_table.compute<4>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
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

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 + 3.0 * rpa_x * rpa_x);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 - 6.0 * rpa_x * rpc_x);

        fints_x[i] += dip_x * fss * b3_vals[i] * 3.0 * rpc_x * rpc_x;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x + rpa_x * rpa_x * rpa_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x - (3.0 / 2.0) * fe_0 * rpc_x - 3.0 * rpa_x * rpa_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_x + 3.0 * rpa_x * rpc_x * rpc_x);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-rpc_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x + rpa_x * rpa_x * rpa_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x - (3.0 / 2.0) * fe_0 * rpc_x - 3.0 * rpa_x * rpa_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_x + 3.0 * rpa_x * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-rpc_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_x + rpa_x * rpa_x * rpa_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_x - (3.0 / 2.0) * fe_0 * rpc_x - 3.0 * rpa_x * rpa_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_x + 3.0 * rpa_x * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-rpc_x * rpc_x * rpc_x);

    }
}

auto
compPrimitiveNuclearPotentialGeom010FS_XXY_0(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<3> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<4> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

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

    bf_table.compute<4>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
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

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * 2.0 * rpa_y * rpa_x;

        fints_x[i] += dip_x * fss * b2_vals[i] * (-2.0 * rpa_y * rpc_x - 2.0 * rpa_x * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] * 2.0 * rpc_y * rpc_x;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_x * rpa_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_y * rpa_x * rpc_x - rpa_x * rpa_x * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + rpa_y * rpc_x * rpc_x + 2.0 * rpa_x * rpc_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-rpc_y * rpc_x * rpc_x);

        fints_y[i] += dip_y * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 + rpa_x * rpa_x);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 - 2.0 * rpa_x * rpc_x);

        fints_y[i] += dip_y * fss * b3_vals[i] * rpc_x * rpc_x;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_x * rpa_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_y * rpa_x * rpc_x - rpa_x * rpa_x * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + rpa_y * rpc_x * rpc_x + 2.0 * rpa_x * rpc_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-rpc_y * rpc_x * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_x * rpa_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_y * rpa_x * rpc_x - rpa_x * rpa_x * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + rpa_y * rpc_x * rpc_x + 2.0 * rpa_x * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-rpc_y * rpc_x * rpc_x);

    }
}

auto
compPrimitiveNuclearPotentialGeom010FS_XXZ_0(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<3> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<4> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

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

    bf_table.compute<4>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
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

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * 2.0 * rpa_z * rpa_x;

        fints_x[i] += dip_x * fss * b2_vals[i] * (-2.0 * rpa_z * rpc_x - 2.0 * rpa_x * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * 2.0 * rpc_z * rpc_x;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_x * rpa_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_z * rpa_x * rpc_x - rpa_x * rpa_x * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + rpa_z * rpc_x * rpc_x + 2.0 * rpa_x * rpc_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-rpc_z * rpc_x * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_x * rpa_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_z * rpa_x * rpc_x - rpa_x * rpa_x * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + rpa_z * rpc_x * rpc_x + 2.0 * rpa_x * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-rpc_z * rpc_x * rpc_x);

        fints_z[i] += dip_z * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 + rpa_x * rpa_x);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 - 2.0 * rpa_x * rpc_x);

        fints_z[i] += dip_z * fss * b3_vals[i] * rpc_x * rpc_x;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_x * rpa_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_z * rpa_x * rpc_x - rpa_x * rpa_x * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + rpa_z * rpc_x * rpc_x + 2.0 * rpa_x * rpc_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-rpc_z * rpc_x * rpc_x);

    }
}

auto
compPrimitiveNuclearPotentialGeom010FS_XYY_0(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<3> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<4> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

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

    bf_table.compute<4>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
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

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 + rpa_y * rpa_y);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 - 2.0 * rpa_y * rpc_y);

        fints_x[i] += dip_x * fss * b3_vals[i] * rpc_y * rpc_y;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_y * rpa_y * rpa_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_y * rpa_x * rpc_y - rpa_y * rpa_y * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_y * rpc_y * rpc_x + rpa_x * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-rpc_y * rpc_y * rpc_x);

        fints_y[i] += dip_y * fss * b1_vals[i] * 2.0 * rpa_y * rpa_x;

        fints_y[i] += dip_y * fss * b2_vals[i] * (-2.0 * rpa_y * rpc_x - 2.0 * rpa_x * rpc_y);

        fints_y[i] += dip_y * fss * b3_vals[i] * 2.0 * rpc_y * rpc_x;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_y * rpa_y * rpa_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_y * rpa_x * rpc_y - rpa_y * rpa_y * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_y * rpc_y * rpc_x + rpa_x * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-rpc_y * rpc_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_y * rpa_y * rpa_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_y * rpa_x * rpc_y - rpa_y * rpa_y * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_y * rpc_y * rpc_x + rpa_x * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-rpc_y * rpc_y * rpc_x);

    }
}

auto
compPrimitiveNuclearPotentialGeom010FS_XYZ_0(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<3> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<4> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

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

    bf_table.compute<4>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
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

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * rpa_z * rpa_y;

        fints_x[i] += dip_x * fss * b2_vals[i] * (-rpa_z * rpc_y - rpa_y * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * rpc_z * rpc_y;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * rpa_z * rpa_y * rpa_x;

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-rpa_z * rpa_y * rpc_x - rpa_z * rpa_x * rpc_y - rpa_y * rpa_x * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * (rpa_z * rpc_y * rpc_x + rpa_y * rpc_z * rpc_x + rpa_x * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-rpc_z * rpc_y * rpc_x);

        fints_y[i] += dip_y * fss * b1_vals[i] * rpa_z * rpa_x;

        fints_y[i] += dip_y * fss * b2_vals[i] * (-rpa_z * rpc_x - rpa_x * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * rpc_z * rpc_x;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * rpa_z * rpa_y * rpa_x;

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-rpa_z * rpa_y * rpc_x - rpa_z * rpa_x * rpc_y - rpa_y * rpa_x * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * (rpa_z * rpc_y * rpc_x + rpa_y * rpc_z * rpc_x + rpa_x * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-rpc_z * rpc_y * rpc_x);

        fints_z[i] += dip_z * fss * b1_vals[i] * rpa_y * rpa_x;

        fints_z[i] += dip_z * fss * b2_vals[i] * (-rpa_y * rpc_x - rpa_x * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] * rpc_y * rpc_x;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * rpa_z * rpa_y * rpa_x;

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-rpa_z * rpa_y * rpc_x - rpa_z * rpa_x * rpc_y - rpa_y * rpa_x * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * (rpa_z * rpc_y * rpc_x + rpa_y * rpc_z * rpc_x + rpa_x * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-rpc_z * rpc_y * rpc_x);

    }
}

auto
compPrimitiveNuclearPotentialGeom010FS_XZZ_0(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<3> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<4> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

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

    bf_table.compute<4>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
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

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        fints_x[i] += dip_x * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 + rpa_z * rpa_z);

        fints_x[i] += dip_x * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 - 2.0 * rpa_z * rpc_z);

        fints_x[i] += dip_x * fss * b3_vals[i] * rpc_z * rpc_z;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_z * rpa_z * rpa_x);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_z * rpa_x * rpc_z - rpa_z * rpa_z * rpc_x);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_z * rpc_z * rpc_x + rpa_x * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-rpc_z * rpc_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_z * rpa_z * rpa_x);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_z * rpa_x * rpc_z - rpa_z * rpa_z * rpc_x);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_z * rpc_z * rpc_x + rpa_x * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-rpc_z * rpc_z * rpc_x);

        fints_z[i] += dip_z * fss * b1_vals[i] * 2.0 * rpa_z * rpa_x;

        fints_z[i] += dip_z * fss * b2_vals[i] * (-2.0 * rpa_z * rpc_x - 2.0 * rpa_x * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * 2.0 * rpc_z * rpc_x;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_z * rpa_z * rpa_x);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_x - (1.0 / 2.0) * fe_0 * rpc_x - 2.0 * rpa_z * rpa_x * rpc_z - rpa_z * rpa_z * rpc_x);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_x + 2.0 * rpa_z * rpc_z * rpc_x + rpa_x * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-rpc_z * rpc_z * rpc_x);

    }
}

auto
compPrimitiveNuclearPotentialGeom010FS_YYY_0(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<3> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<4> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

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

    bf_table.compute<4>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
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

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_y * rpa_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y - (3.0 / 2.0) * fe_0 * rpc_y - 3.0 * rpa_y * rpa_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y + 3.0 * rpa_y * rpc_y * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-rpc_y * rpc_y * rpc_y);

        fints_y[i] += dip_y * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 + 3.0 * rpa_y * rpa_y);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 - 6.0 * rpa_y * rpc_y);

        fints_y[i] += dip_y * fss * b3_vals[i] * 3.0 * rpc_y * rpc_y;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_y * rpa_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y - (3.0 / 2.0) * fe_0 * rpc_y - 3.0 * rpa_y * rpa_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y + 3.0 * rpa_y * rpc_y * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-rpc_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_y * rpa_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_y - (3.0 / 2.0) * fe_0 * rpc_y - 3.0 * rpa_y * rpa_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_y + 3.0 * rpa_y * rpc_y * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-rpc_y * rpc_y * rpc_y);

    }
}

auto
compPrimitiveNuclearPotentialGeom010FS_YYZ_0(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<3> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<4> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

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

    bf_table.compute<4>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
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

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_y * rpa_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_z * rpa_y * rpc_y - rpa_y * rpa_y * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + rpa_z * rpc_y * rpc_y + 2.0 * rpa_y * rpc_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-rpc_z * rpc_y * rpc_y);

        fints_y[i] += dip_y * fss * b1_vals[i] * 2.0 * rpa_z * rpa_y;

        fints_y[i] += dip_y * fss * b2_vals[i] * (-2.0 * rpa_z * rpc_y - 2.0 * rpa_y * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * 2.0 * rpc_z * rpc_y;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_y * rpa_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_z * rpa_y * rpc_y - rpa_y * rpa_y * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + rpa_z * rpc_y * rpc_y + 2.0 * rpa_y * rpc_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-rpc_z * rpc_y * rpc_y);

        fints_z[i] += dip_z * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 + rpa_y * rpa_y);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 - 2.0 * rpa_y * rpc_y);

        fints_z[i] += dip_z * fss * b3_vals[i] * rpc_y * rpc_y;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_y * rpa_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_z - (1.0 / 2.0) * fe_0 * rpc_z - 2.0 * rpa_z * rpa_y * rpc_y - rpa_y * rpa_y * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_z + rpa_z * rpc_y * rpc_y + 2.0 * rpa_y * rpc_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-rpc_z * rpc_y * rpc_y);

    }
}

auto
compPrimitiveNuclearPotentialGeom010FS_YZZ_0(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<3> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<4> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

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

    bf_table.compute<4>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
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

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_z * rpa_z * rpa_y);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_z * rpa_y * rpc_z - rpa_z * rpa_z * rpc_y);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + 2.0 * rpa_z * rpc_z * rpc_y + rpa_y * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-rpc_z * rpc_z * rpc_y);

        fints_y[i] += dip_y * fss * b1_vals[i] * ((1.0 / 2.0) * fe_0 + rpa_z * rpa_z);

        fints_y[i] += dip_y * fss * b2_vals[i] * (-(1.0 / 2.0) * fe_0 - 2.0 * rpa_z * rpc_z);

        fints_y[i] += dip_y * fss * b3_vals[i] * rpc_z * rpc_z;

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_z * rpa_z * rpa_y);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_z * rpa_y * rpc_z - rpa_z * rpa_z * rpc_y);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + 2.0 * rpa_z * rpc_z * rpc_y + rpa_y * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-rpc_z * rpc_z * rpc_y);

        fints_z[i] += dip_z * fss * b1_vals[i] * 2.0 * rpa_z * rpa_y;

        fints_z[i] += dip_z * fss * b2_vals[i] * (-2.0 * rpa_z * rpc_y - 2.0 * rpa_y * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * 2.0 * rpc_z * rpc_y;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_z * rpa_z * rpa_y);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(1.0 / 2.0) * fe_0 * rpa_y - (1.0 / 2.0) * fe_0 * rpc_y - 2.0 * rpa_z * rpa_y * rpc_z - rpa_z * rpa_z * rpc_y);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((1.0 / 2.0) * fe_0 * rpc_y + 2.0 * rpa_z * rpc_z * rpc_y + rpa_y * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-rpc_z * rpc_z * rpc_y);

    }
}

auto
compPrimitiveNuclearPotentialGeom010FS_ZZZ_0(      TDoubleArray& buffer_x,
                                                   TDoubleArray& buffer_y,
                                                   TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
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

    // set up dipole components

    const auto dip_x = dipole[0];

    const auto dip_y = dipole[1];

    const auto dip_z = dipole[2];

    // set up pointer to integrals buffer(s)

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

    // set up Boys function variables

    const CBoysFunc<3> bf_table;

    alignas(64) TDoubleArray bf_args;

    TDoubleArray2D<4> bf_values;

    auto b0_vals = bf_values[0].data();

    auto b1_vals = bf_values[1].data();

    auto b2_vals = bf_values[2].data();

    auto b3_vals = bf_values[3].data();

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

    bf_table.compute<4>(bf_values, bf_args, ket_dim);

    #pragma omp simd aligned(fints_x,\
                             fints_y,\
                             fints_z,\
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

        const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);

        const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);

        const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);

        const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;

        const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;

        const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        fints_x[i] += dip_x * faa_x * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_z * rpa_z);

        fints_x[i] += dip_x * faa_x * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z - (3.0 / 2.0) * fe_0 * rpc_z - 3.0 * rpa_z * rpa_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z + 3.0 * rpa_z * rpc_z * rpc_z);

        fints_x[i] += dip_x * faa_x * b3_vals[i] * (-rpc_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_z * rpa_z);

        fints_y[i] += dip_y * faa_y * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z - (3.0 / 2.0) * fe_0 * rpc_z - 3.0 * rpa_z * rpa_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z + 3.0 * rpa_z * rpc_z * rpc_z);

        fints_y[i] += dip_y * faa_y * b3_vals[i] * (-rpc_z * rpc_z * rpc_z);

        fints_z[i] += dip_z * fss * b1_vals[i] * ((3.0 / 2.0) * fe_0 + 3.0 * rpa_z * rpa_z);

        fints_z[i] += dip_z * fss * b2_vals[i] * (-(3.0 / 2.0) * fe_0 - 6.0 * rpa_z * rpc_z);

        fints_z[i] += dip_z * fss * b3_vals[i] * 3.0 * rpc_z * rpc_z;

        fints_z[i] += dip_z * faa_z * b0_vals[i] * ((3.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_z * rpa_z);

        fints_z[i] += dip_z * faa_z * b1_vals[i] * (-(3.0 / 2.0) * fe_0 * rpa_z - (3.0 / 2.0) * fe_0 * rpc_z - 3.0 * rpa_z * rpa_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b2_vals[i] * ((3.0 / 2.0) * fe_0 * rpc_z + 3.0 * rpa_z * rpc_z * rpc_z);

        fints_z[i] += dip_z * faa_z * b3_vals[i] * (-rpc_z * rpc_z * rpc_z);

    }
}

} // geom_npotrec namespace

