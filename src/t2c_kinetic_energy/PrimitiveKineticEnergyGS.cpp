#include "PrimitiveKineticEnergyGS.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace kinrec {  // kinrec namespace

auto
compPrimitiveKineticEnergyGS(TDoubleArray&       buffer_xxxx,
                             TDoubleArray&       buffer_xxxy,
                             TDoubleArray&       buffer_xxxz,
                             TDoubleArray&       buffer_xxyy,
                             TDoubleArray&       buffer_xxyz,
                             TDoubleArray&       buffer_xxzz,
                             TDoubleArray&       buffer_xyyy,
                             TDoubleArray&       buffer_xyyz,
                             TDoubleArray&       buffer_xyzz,
                             TDoubleArray&       buffer_xzzz,
                             TDoubleArray&       buffer_yyyy,
                             TDoubleArray&       buffer_yyyz,
                             TDoubleArray&       buffer_yyzz,
                             TDoubleArray&       buffer_yzzz,
                             TDoubleArray&       buffer_zzzz,
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

    auto fints_xxxx = buffer_xxxx.data();

    auto fints_xxxy = buffer_xxxy.data();

    auto fints_xxxz = buffer_xxxz.data();

    auto fints_xxyy = buffer_xxyy.data();

    auto fints_xxyz = buffer_xxyz.data();

    auto fints_xxzz = buffer_xxzz.data();

    auto fints_xyyy = buffer_xyyy.data();

    auto fints_xyyz = buffer_xyyz.data();

    auto fints_xyzz = buffer_xyzz.data();

    auto fints_xzzz = buffer_xzzz.data();

    auto fints_yyyy = buffer_yyyy.data();

    auto fints_yyyz = buffer_yyyz.data();

    auto fints_yyzz = buffer_yyzz.data();

    auto fints_yzzz = buffer_yzzz.data();

    auto fints_zzzz = buffer_zzzz.data();

#pragma omp simd aligned(fints_xxxx,     \
                             fints_xxxy, \
                             fints_xxxz, \
                             fints_xxyy, \
                             fints_xxyz, \
                             fints_xxzz, \
                             fints_xyyy, \
                             fints_xyyz, \
                             fints_xyzz, \
                             fints_xzzz, \
                             fints_yyyy, \
                             fints_yyyz, \
                             fints_yyzz, \
                             fints_yzzz, \
                             fints_zzzz, \
                             ket_fe,     \
                             ket_fn,     \
                             ket_rx,     \
                             ket_ry,     \
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

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxxx[i] += fss * (-3.0 * fbe_0 * fe_0 * fz_0 - 6.0 * fbe_0 * rpa_x * rpa_x * fz_0 + 18.0 * fe_0 * rpa_x * rpa_x * fz_0 +
                                3.0 * fe_0 * fe_0 * fz_0 + 8.0 * rpa_x * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxxx[i] += ftt * (3.0 * fe_0 * rpa_x * rpa_x + (3.0 / 4.0) * fe_0 * fe_0 + rpa_x * rpa_x * rpa_x * rpa_x);

        fints_xxxy[i] += fss * (-3.0 * fbe_0 * rpa_y * rpa_x * fz_0 + 9.0 * fe_0 * rpa_y * rpa_x * fz_0 + 8.0 * rpa_y * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxxy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x + rpa_y * rpa_x * rpa_x * rpa_x);

        fints_xxxz[i] += fss * (-3.0 * fbe_0 * rpa_z * rpa_x * fz_0 + 9.0 * fe_0 * rpa_z * rpa_x * fz_0 + 8.0 * rpa_z * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxxz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x + rpa_z * rpa_x * rpa_x * rpa_x);

        fints_xxyy[i] += fss * (-fbe_0 * fe_0 * fz_0 - fbe_0 * rpa_y * rpa_y * fz_0 - fbe_0 * rpa_x * rpa_x * fz_0 +
                                3.0 * fe_0 * rpa_y * rpa_y * fz_0 + 3.0 * fe_0 * rpa_x * rpa_x * fz_0);

        fints_xxyy[i] += fss * (fe_0 * fe_0 * fz_0 + 8.0 * rpa_y * rpa_y * rpa_x * rpa_x * fz_0);

        fints_xxyy[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_y + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 +
                                rpa_y * rpa_y * rpa_x * rpa_x);

        fints_xxyz[i] += fss * (-fbe_0 * rpa_z * rpa_y * fz_0 + 3.0 * fe_0 * rpa_z * rpa_y * fz_0 + 8.0 * rpa_z * rpa_y * rpa_x * rpa_x * fz_0);

        fints_xxyz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_y + rpa_z * rpa_y * rpa_x * rpa_x);

        fints_xxzz[i] += fss * (-fbe_0 * fe_0 * fz_0 - fbe_0 * rpa_z * rpa_z * fz_0 - fbe_0 * rpa_x * rpa_x * fz_0 +
                                3.0 * fe_0 * rpa_z * rpa_z * fz_0 + 3.0 * fe_0 * rpa_x * rpa_x * fz_0);

        fints_xxzz[i] += fss * (fe_0 * fe_0 * fz_0 + 8.0 * rpa_z * rpa_z * rpa_x * rpa_x * fz_0);

        fints_xxzz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z + (1.0 / 2.0) * fe_0 * rpa_x * rpa_x + (1.0 / 4.0) * fe_0 * fe_0 +
                                rpa_z * rpa_z * rpa_x * rpa_x);

        fints_xyyy[i] += fss * (-3.0 * fbe_0 * rpa_y * rpa_x * fz_0 + 9.0 * fe_0 * rpa_y * rpa_x * fz_0 + 8.0 * rpa_y * rpa_y * rpa_y * rpa_x * fz_0);

        fints_xyyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y * rpa_x + rpa_y * rpa_y * rpa_y * rpa_x);

        fints_xyyz[i] += fss * (-fbe_0 * rpa_z * rpa_x * fz_0 + 3.0 * fe_0 * rpa_z * rpa_x * fz_0 + 8.0 * rpa_z * rpa_y * rpa_y * rpa_x * fz_0);

        fints_xyyz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_x + rpa_z * rpa_y * rpa_y * rpa_x);

        fints_xyzz[i] += fss * (-fbe_0 * rpa_y * rpa_x * fz_0 + 3.0 * fe_0 * rpa_y * rpa_x * fz_0 + 8.0 * rpa_z * rpa_z * rpa_y * rpa_x * fz_0);

        fints_xyzz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_y * rpa_x + rpa_z * rpa_z * rpa_y * rpa_x);

        fints_xzzz[i] += fss * (-3.0 * fbe_0 * rpa_z * rpa_x * fz_0 + 9.0 * fe_0 * rpa_z * rpa_x * fz_0 + 8.0 * rpa_z * rpa_z * rpa_z * rpa_x * fz_0);

        fints_xzzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_x + rpa_z * rpa_z * rpa_z * rpa_x);

        fints_yyyy[i] += fss * (-3.0 * fbe_0 * fe_0 * fz_0 - 6.0 * fbe_0 * rpa_y * rpa_y * fz_0 + 18.0 * fe_0 * rpa_y * rpa_y * fz_0 +
                                3.0 * fe_0 * fe_0 * fz_0 + 8.0 * rpa_y * rpa_y * rpa_y * rpa_y * fz_0);

        fints_yyyy[i] += ftt * (3.0 * fe_0 * rpa_y * rpa_y + (3.0 / 4.0) * fe_0 * fe_0 + rpa_y * rpa_y * rpa_y * rpa_y);

        fints_yyyz[i] += fss * (-3.0 * fbe_0 * rpa_z * rpa_y * fz_0 + 9.0 * fe_0 * rpa_z * rpa_y * fz_0 + 8.0 * rpa_z * rpa_y * rpa_y * rpa_y * fz_0);

        fints_yyyz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y + rpa_z * rpa_y * rpa_y * rpa_y);

        fints_yyzz[i] += fss * (-fbe_0 * fe_0 * fz_0 - fbe_0 * rpa_z * rpa_z * fz_0 - fbe_0 * rpa_y * rpa_y * fz_0 +
                                3.0 * fe_0 * rpa_z * rpa_z * fz_0 + 3.0 * fe_0 * rpa_y * rpa_y * fz_0);

        fints_yyzz[i] += fss * (fe_0 * fe_0 * fz_0 + 8.0 * rpa_z * rpa_z * rpa_y * rpa_y * fz_0);

        fints_yyzz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z * rpa_z + (1.0 / 2.0) * fe_0 * rpa_y * rpa_y + (1.0 / 4.0) * fe_0 * fe_0 +
                                rpa_z * rpa_z * rpa_y * rpa_y);

        fints_yzzz[i] += fss * (-3.0 * fbe_0 * rpa_z * rpa_y * fz_0 + 9.0 * fe_0 * rpa_z * rpa_y * fz_0 + 8.0 * rpa_z * rpa_z * rpa_z * rpa_y * fz_0);

        fints_yzzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z * rpa_y + rpa_z * rpa_z * rpa_z * rpa_y);

        fints_zzzz[i] += fss * (-3.0 * fbe_0 * fe_0 * fz_0 - 6.0 * fbe_0 * rpa_z * rpa_z * fz_0 + 18.0 * fe_0 * rpa_z * rpa_z * fz_0 +
                                3.0 * fe_0 * fe_0 * fz_0 + 8.0 * rpa_z * rpa_z * rpa_z * rpa_z * fz_0);

        fints_zzzz[i] += ftt * (3.0 * fe_0 * rpa_z * rpa_z + (3.0 / 4.0) * fe_0 * fe_0 + rpa_z * rpa_z * rpa_z * rpa_z);
    }
}

}  // namespace kinrec
