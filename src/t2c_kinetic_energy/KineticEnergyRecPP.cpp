#include "KineticEnergyRecPP.hpp"

#include <cmath>

#include "BatchFunc.hpp"
#include "MathConst.hpp"
#include "T2CDistributor.hpp"

namespace kinrec {  // kinrec namespace

auto
compKineticEnergyPP(CSubMatrix* matrix, const CGtoBlock& gto_block, const int64_t bra_first, const int64_t bra_last) -> void
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

            // compute primitive integrals block (X)

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

                    kinrec::compPrimitiveKineticEnergyPP_X_T(buffer_x,
                                                             buffer_y,
                                                             buffer_z,
                                                             bra_exp,
                                                             bra_norm,
                                                             bra_coord,
                                                             ket_exps,
                                                             ket_norms,
                                                             ket_coords_x,
                                                             ket_coords_y,
                                                             ket_coords_z,
                                                             ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, gto_indexes, 2, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_y, gto_indexes, 2, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_z, gto_indexes, 2, 1, j, ket_first, ket_last);

            // compute primitive integrals block (Y)

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

                    kinrec::compPrimitiveKineticEnergyPP_Y_T(buffer_x,
                                                             buffer_y,
                                                             buffer_z,
                                                             bra_exp,
                                                             bra_norm,
                                                             bra_coord,
                                                             ket_exps,
                                                             ket_norms,
                                                             ket_coords_x,
                                                             ket_coords_y,
                                                             ket_coords_z,
                                                             ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, gto_indexes, 0, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_y, gto_indexes, 0, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_z, gto_indexes, 0, 1, j, ket_first, ket_last);

            // compute primitive integrals block (Z)

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

                    kinrec::compPrimitiveKineticEnergyPP_Z_T(buffer_x,
                                                             buffer_y,
                                                             buffer_z,
                                                             bra_exp,
                                                             bra_norm,
                                                             bra_coord,
                                                             ket_exps,
                                                             ket_norms,
                                                             ket_coords_x,
                                                             ket_coords_y,
                                                             ket_coords_z,
                                                             ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, gto_indexes, 1, 2, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_y, gto_indexes, 1, 0, j, ket_first, ket_last);

            t2cfunc::distribute(matrix, buffer_z, gto_indexes, 1, 1, j, ket_first, ket_last);
        }
    }
}

auto
compKineticEnergyPP(CSubMatrix*      matrix,
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

            // compute primitive integrals block (X)

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

                    kinrec::compPrimitiveKineticEnergyPP_X_T(buffer_x,
                                                             buffer_y,
                                                             buffer_z,
                                                             bra_exp,
                                                             bra_norm,
                                                             bra_coord,
                                                             ket_exps,
                                                             ket_norms,
                                                             ket_coords_x,
                                                             ket_coords_y,
                                                             ket_coords_z,
                                                             ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, bra_gto_indexes, ket_gto_indexes, 2, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_y, bra_gto_indexes, ket_gto_indexes, 2, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_z, bra_gto_indexes, ket_gto_indexes, 2, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Y)

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

                    kinrec::compPrimitiveKineticEnergyPP_Y_T(buffer_x,
                                                             buffer_y,
                                                             buffer_z,
                                                             bra_exp,
                                                             bra_norm,
                                                             bra_coord,
                                                             ket_exps,
                                                             ket_norms,
                                                             ket_coords_x,
                                                             ket_coords_y,
                                                             ket_coords_z,
                                                             ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, bra_gto_indexes, ket_gto_indexes, 0, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_y, bra_gto_indexes, ket_gto_indexes, 0, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_z, bra_gto_indexes, ket_gto_indexes, 0, 1, j, ket_first, ket_last, mat_type);

            // compute primitive integrals block (Z)

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

                    kinrec::compPrimitiveKineticEnergyPP_Z_T(buffer_x,
                                                             buffer_y,
                                                             buffer_z,
                                                             bra_exp,
                                                             bra_norm,
                                                             bra_coord,
                                                             ket_exps,
                                                             ket_norms,
                                                             ket_coords_x,
                                                             ket_coords_y,
                                                             ket_coords_z,
                                                             ket_dim);
                }
            }

            t2cfunc::distribute(matrix, buffer_x, bra_gto_indexes, ket_gto_indexes, 1, 2, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_y, bra_gto_indexes, ket_gto_indexes, 1, 0, j, ket_first, ket_last, mat_type);

            t2cfunc::distribute(matrix, buffer_z, bra_gto_indexes, ket_gto_indexes, 1, 1, j, ket_first, ket_last, mat_type);
        }
    }
}

auto
compPrimitiveKineticEnergyPP_X_T(TDoubleArray&       buffer_x,
                                 TDoubleArray&       buffer_y,
                                 TDoubleArray&       buffer_z,
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

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
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

        fints_x[i] += fss * (fe_0 * fz_0 + 4.0 * rpa_x * rpb_x * fz_0);

        fints_x[i] += ftt * ((1.0 / 2.0) * fe_0 + rpa_x * rpb_x);

        fints_y[i] += fss * 4.0 * rpa_x * rpb_y * fz_0;

        fints_y[i] += ftt * rpa_x * rpb_y;

        fints_z[i] += fss * 4.0 * rpa_x * rpb_z * fz_0;

        fints_z[i] += ftt * rpa_x * rpb_z;
    }
}

auto
compPrimitiveKineticEnergyPP_Y_T(TDoubleArray&       buffer_x,
                                 TDoubleArray&       buffer_y,
                                 TDoubleArray&       buffer_z,
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

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
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

        fints_x[i] += fss * 4.0 * rpa_y * rpb_x * fz_0;

        fints_x[i] += ftt * rpa_y * rpb_x;

        fints_y[i] += fss * (fe_0 * fz_0 + 4.0 * rpa_y * rpb_y * fz_0);

        fints_y[i] += ftt * ((1.0 / 2.0) * fe_0 + rpa_y * rpb_y);

        fints_z[i] += fss * 4.0 * rpa_y * rpb_z * fz_0;

        fints_z[i] += ftt * rpa_y * rpb_z;
    }
}

auto
compPrimitiveKineticEnergyPP_Z_T(TDoubleArray&       buffer_x,
                                 TDoubleArray&       buffer_y,
                                 TDoubleArray&       buffer_z,
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

    auto fints_x = buffer_x.data();

    auto fints_y = buffer_y.data();

    auto fints_z = buffer_z.data();

#pragma omp simd aligned(fints_x, fints_y, fints_z, ket_fe, ket_fn, ket_rx, ket_ry, ket_rz : 64)
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

        fints_x[i] += fss * 4.0 * rpa_z * rpb_x * fz_0;

        fints_x[i] += ftt * rpa_z * rpb_x;

        fints_y[i] += fss * 4.0 * rpa_z * rpb_y * fz_0;

        fints_y[i] += ftt * rpa_z * rpb_y;

        fints_z[i] += fss * (fe_0 * fz_0 + 4.0 * rpa_z * rpb_z * fz_0);

        fints_z[i] += ftt * ((1.0 / 2.0) * fe_0 + rpa_z * rpb_z);
    }
}

}  // namespace kinrec
