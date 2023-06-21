#include "OverlapRecSF.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "MathConst.hpp"
#include "T2CDistributor.hpp"

namespace ovlrec { // ovlrec namespace

auto
compOverlapSF(      CSubMatrix* matrix,
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

            // compute primitive integrals block

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

                    ovlrec::compPrimitiveOverlapSF(buffer_xxx,
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

            t2cfunc::distribute(matrix, buffer_xxx, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxx, f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, 3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xxz, 0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyy, -3.0 * f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 6, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xyz, f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 1, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_xzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 4, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_5, bra_gto_indexes, ket_gto_indexes,
                                0, 0, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyy, -f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -3.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yyz, -0.5 * f3_15, bra_gto_indexes, ket_gto_indexes,
                                0, 5, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_yzz, 4.0 * f3_3, bra_gto_indexes, ket_gto_indexes,
                                0, 2, j, ket_first, ket_last, ang_order);

            t2cfunc::distribute(matrix, buffer_zzz, 2.0, bra_gto_indexes, ket_gto_indexes,
                                0, 3, j, ket_first, ket_last, ang_order);

        }
    }
}

auto
compPrimitiveOverlapSF(      TDoubleArray& buffer_xxx,
                             TDoubleArray& buffer_xxy,
                             TDoubleArray& buffer_xxz,
                             TDoubleArray& buffer_xyy,
                             TDoubleArray& buffer_xyz,
                             TDoubleArray& buffer_xzz,
                             TDoubleArray& buffer_yyy,
                             TDoubleArray& buffer_yyz,
                             TDoubleArray& buffer_yzz,
                             TDoubleArray& buffer_zzz,
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

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        fints_xxx[i] += fss * ((3.0 / 2.0) * fe_0 * rpb_x + rpb_x * rpb_x * rpb_x);

        fints_xxy[i] += fss * ((1.0 / 2.0) * fe_0 * rpb_y + rpb_y * rpb_x * rpb_x);

        fints_xxz[i] += fss * ((1.0 / 2.0) * fe_0 * rpb_z + rpb_z * rpb_x * rpb_x);

        fints_xyy[i] += fss * ((1.0 / 2.0) * fe_0 * rpb_x + rpb_y * rpb_y * rpb_x);

        fints_xyz[i] += fss * rpb_z * rpb_y * rpb_x;

        fints_xzz[i] += fss * ((1.0 / 2.0) * fe_0 * rpb_x + rpb_z * rpb_z * rpb_x);

        fints_yyy[i] += fss * ((3.0 / 2.0) * fe_0 * rpb_y + rpb_y * rpb_y * rpb_y);

        fints_yyz[i] += fss * ((1.0 / 2.0) * fe_0 * rpb_z + rpb_z * rpb_y * rpb_y);

        fints_yzz[i] += fss * ((1.0 / 2.0) * fe_0 * rpb_y + rpb_z * rpb_z * rpb_y);

        fints_zzz[i] += fss * ((3.0 / 2.0) * fe_0 * rpb_z + rpb_z * rpb_z * rpb_z);
    }
}

} // ovlrec namespace

