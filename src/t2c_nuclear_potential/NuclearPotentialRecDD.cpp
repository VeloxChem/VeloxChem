#include "NuclearPotentialRecDD.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "PrimitiveNuclearPotentialDD_XX_T.hpp"
#include "PrimitiveNuclearPotentialDD_XY_T.hpp"
#include "PrimitiveNuclearPotentialDD_XZ_T.hpp"
#include "PrimitiveNuclearPotentialDD_YY_T.hpp"
#include "PrimitiveNuclearPotentialDD_YZ_T.hpp"
#include "PrimitiveNuclearPotentialDD_ZZ_T.hpp"
#include "T2CDistributor.hpp"

namespace npotrec {  // npotrec namespace

auto
compNuclearPotentialDD(CSubMatrix*      matrix,
                       const double     charge,
                       const TPoint3D&  point,
                       const CGtoBlock& gto_block,
                       const int64_t    bra_first,
                       const int64_t    bra_last) -> void
{
    // spherical transformation factors

    const double f2_3 = 2.0 * std::sqrt(3.0);

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

    alignas(64) TDoubleArray buffer_xx;

    alignas(64) TDoubleArray buffer_xy;

    alignas(64) TDoubleArray buffer_xz;

    alignas(64) TDoubleArray buffer_yy;

    alignas(64) TDoubleArray buffer_yz;

    alignas(64) TDoubleArray buffer_zz;

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

            // compute primitive integrals block (XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_XX_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xx, -0.5 * f2_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xx, -0.5 * f2_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xx, 0.5 * f2_3 * 0.5 * f2_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xy, -f2_3, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xy, 0.5 * f2_3 * f2_3, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xz, -f2_3, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xz, 0.5 * f2_3 * f2_3, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, -0.5 * f2_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, 0.5 * f2_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, -0.5 * f2_3 * 0.5 * f2_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yz, -f2_3, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yz, 0.5 * f2_3 * f2_3, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zz, -2.0, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zz, 0.5 * f2_3 * 2.0, gto_indexes, 4, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_XY_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, -f2_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xx, f2_3 * 0.5 * f2_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xy, f2_3 * f2_3, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xz, f2_3 * f2_3, gto_indexes, 0, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3 * 0.5 * f2_3, gto_indexes, 0, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yz, f2_3 * f2_3, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zz, f2_3 * 2.0, gto_indexes, 0, 2, j, ket_first, ket_last);

            // compute primitive integrals block (XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_XZ_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, -f2_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xx, f2_3 * 0.5 * f2_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xy, f2_3 * f2_3, gto_indexes, 3, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xz, f2_3 * f2_3, gto_indexes, 3, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3, gto_indexes, 3, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3 * 0.5 * f2_3, gto_indexes, 3, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yz, f2_3 * f2_3, gto_indexes, 3, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zz, f2_3 * 2.0, gto_indexes, 3, 2, j, ket_first, ket_last);

            // compute primitive integrals block (YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_YY_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xx, 0.5 * f2_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xx, -0.5 * f2_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xx, -0.5 * f2_3 * 0.5 * f2_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xy, -f2_3, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xy, -0.5 * f2_3 * f2_3, gto_indexes, 4, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xz, -f2_3, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xz, -0.5 * f2_3 * f2_3, gto_indexes, 4, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, 0.5 * f2_3, gto_indexes, 4, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, 0.5 * f2_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, 0.5 * f2_3 * 0.5 * f2_3, gto_indexes, 4, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yz, -f2_3, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yz, -0.5 * f2_3 * f2_3, gto_indexes, 4, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zz, -2.0, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zz, -0.5 * f2_3 * 2.0, gto_indexes, 4, 2, j, ket_first, ket_last);

            // compute primitive integrals block (YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_YZ_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, -f2_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xx, f2_3 * 0.5 * f2_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xy, f2_3 * f2_3, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xz, f2_3 * f2_3, gto_indexes, 1, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3 * 0.5 * f2_3, gto_indexes, 1, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yz, f2_3 * f2_3, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zz, f2_3 * 2.0, gto_indexes, 1, 2, j, ket_first, ket_last);

            // compute primitive integrals block (ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto bra_index = l * ncgtos + j;

                    const auto bra_exp = gto_exps[bra_index];

                    const auto bra_norm = gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_ZZ_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, -2.0, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xx, 2.0 * 0.5 * f2_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xy, 2.0 * f2_3, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_xz, 2.0 * f2_3, gto_indexes, 2, 3, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, -2.0, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yy, -2.0 * 0.5 * f2_3, gto_indexes, 2, 4, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_yz, 2.0 * f2_3, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_zz, 2.0 * 2.0, gto_indexes, 2, 2, j, ket_first, ket_last);
        }
    }
}

auto
compNuclearPotentialDD(CSubMatrix*      matrix,
                       const double     charge,
                       const TPoint3D&  point,
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

    alignas(64) TDoubleArray buffer_xx;

    alignas(64) TDoubleArray buffer_xy;

    alignas(64) TDoubleArray buffer_xz;

    alignas(64) TDoubleArray buffer_yy;

    alignas(64) TDoubleArray buffer_yz;

    alignas(64) TDoubleArray buffer_zz;

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

            // compute primitive integrals block (XX)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_XX_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xx, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xy, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_yy, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yz, 0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zz, 0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_XY_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xx, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (XZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_XZ_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xx, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 3, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 3, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YY)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_YY_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xx, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xx, -0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(
                matrix, buffer_xx, -0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xy, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xy, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, 0.5 * f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yz, -f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yz, -0.5 * f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 4, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zz, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zz, -0.5 * f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 4, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (YZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_YZ_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xx, f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xy, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, -f2_3 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yz, f2_3 * f2_3, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zz, f2_3 * 2.0, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (ZZ)

            simd::zero(buffer_xx);

            simd::zero(buffer_xy);

            simd::zero(buffer_xz);

            simd::zero(buffer_yy);

            simd::zero(buffer_yz);

            simd::zero(buffer_zz);

            for (int64_t k = 0; k < ket_npgtos; k++)
            {
                simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);

                simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);

                for (int64_t l = 0; l < bra_npgtos; l++)
                {
                    const auto bra_index = l * bra_ncgtos + j;

                    const auto bra_exp = bra_gto_exps[bra_index];

                    const auto bra_norm = bra_gto_norms[bra_index];

                    npotrec::compPrimitiveNuclearPotentialDD_ZZ_T(buffer_xx,
                                                                  buffer_xy,
                                                                  buffer_xz,
                                                                  buffer_yy,
                                                                  buffer_yz,
                                                                  buffer_zz,
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

            t2cfunc::distribute(matrix, buffer_xx, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xx, 2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xy, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_xz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 3, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, -2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yy, -2.0 * 0.5 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 4, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_yz, 2.0 * f2_3, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_zz, 2.0 * 2.0, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);
        }
    }
}

}  // namespace npotrec
