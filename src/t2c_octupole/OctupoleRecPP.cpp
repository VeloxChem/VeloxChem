#include "OctupoleRecPP.hpp"

#include "BatchFunc.hpp"
#include "PrimitiveOctupolePP_X_X.hpp"
#include "PrimitiveOctupolePP_X_Y.hpp"
#include "PrimitiveOctupolePP_X_Z.hpp"
#include "PrimitiveOctupolePP_Y_X.hpp"
#include "PrimitiveOctupolePP_Y_Y.hpp"
#include "PrimitiveOctupolePP_Y_Z.hpp"
#include "PrimitiveOctupolePP_Z_X.hpp"
#include "PrimitiveOctupolePP_Z_Y.hpp"
#include "PrimitiveOctupolePP_Z_Z.hpp"
#include "T2CDistributor.hpp"

namespace mpol {  // mpol namespace

auto
compOctupolePP(CSubMatrix*      matrix_xxx,
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
               const CGtoBlock& gto_block,
               const int64_t    bra_first,
               const int64_t    bra_last) -> void
{
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

            // compute primitive integrals block (X_X)

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

                    mpol::compPrimitiveOctupolePP_X_X(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, gto_indexes, 2, 2, j, ket_first, ket_last);

            // compute primitive integrals block (X_Y)

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

                    mpol::compPrimitiveOctupolePP_X_Y(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, gto_indexes, 2, 0, j, ket_first, ket_last);

            // compute primitive integrals block (X_Z)

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

                    mpol::compPrimitiveOctupolePP_X_Z(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, gto_indexes, 2, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, gto_indexes, 2, 1, j, ket_first, ket_last);

            // compute primitive integrals block (Y_X)

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

                    mpol::compPrimitiveOctupolePP_Y_X(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, gto_indexes, 0, 2, j, ket_first, ket_last);

            // compute primitive integrals block (Y_Y)

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

                    mpol::compPrimitiveOctupolePP_Y_Y(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, gto_indexes, 0, 0, j, ket_first, ket_last);

            // compute primitive integrals block (Y_Z)

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

                    mpol::compPrimitiveOctupolePP_Y_Z(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, gto_indexes, 0, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, gto_indexes, 0, 1, j, ket_first, ket_last);

            // compute primitive integrals block (Z_X)

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

                    mpol::compPrimitiveOctupolePP_Z_X(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, gto_indexes, 1, 2, j, ket_first, ket_last);

            // compute primitive integrals block (Z_Y)

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

                    mpol::compPrimitiveOctupolePP_Z_Y(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, gto_indexes, 1, 0, j, ket_first, ket_last);

            // compute primitive integrals block (Z_Z)

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

                    mpol::compPrimitiveOctupolePP_Z_Z(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, gto_indexes, 1, 1, j, ket_first, ket_last);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, gto_indexes, 1, 1, j, ket_first, ket_last);
        }
    }
}

auto
compOctupolePP(CSubMatrix*      matrix_xxx,
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
               const int64_t    bra_first,
               const int64_t    bra_last,
               const mat_t      mat_type) -> void
{
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

            // compute primitive integrals block (X_X)

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

                    mpol::compPrimitiveOctupolePP_X_X(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (X_Y)

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

                    mpol::compPrimitiveOctupolePP_X_Y(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (X_Z)

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

                    mpol::compPrimitiveOctupolePP_X_Z(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Y_X)

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

                    mpol::compPrimitiveOctupolePP_Y_X(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Y_Y)

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

                    mpol::compPrimitiveOctupolePP_Y_Y(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Y_Z)

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

                    mpol::compPrimitiveOctupolePP_Y_Z(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Z_X)

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

                    mpol::compPrimitiveOctupolePP_Z_X(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Z_Y)

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

                    mpol::compPrimitiveOctupolePP_Z_Y(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Z_Z)

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

                    mpol::compPrimitiveOctupolePP_Z_Z(buffer_xxx,
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

            t2cfunc::distribute(matrix_xxx, buffer_xxx, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxy, buffer_xxy, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xxz, buffer_xxz, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyy, buffer_xyy, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xyz, buffer_xyz, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_xzz, buffer_xzz, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyy, buffer_yyy, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yyz, buffer_yyz, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_yzz, buffer_yzz, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix_zzz, buffer_zzz, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);
        }
    }
}

}  // namespace mpol
