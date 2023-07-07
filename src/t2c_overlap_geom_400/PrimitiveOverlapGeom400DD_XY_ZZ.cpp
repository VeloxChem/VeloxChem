#include "PrimitiveOverlapGeom400DD_XY_ZZ.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compPrimitiveOverlapGeom400DD_XY_ZZ(TDoubleArray&       buffer_xxxx,
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

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto tbe_0 = bra_exp;

        fints_xxxx[i] +=
            fss * (30.0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 - 20.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   40.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 - 40.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxxx[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxxx[i] +=
            fss * (16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   60.0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 - 40.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   80.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxxx[i] += fss * (-80.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                32.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxxx[i] += fss * (16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                32.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxxy[i] += fss * (-3.0 * fe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * tbe_0 * tbe_0 +
                                6.0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 + 12.0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0);

        fints_xxxy[i] +=
            fss * (-6.0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 - 3.0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 -
                   12.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 - 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxxy[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
             3.0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tbe_0 - 24.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
             4.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxxy[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxxy[i] +=
            fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 - 6.0 * rpb_z * rpb_z * tbe_0);

        fints_xxxy[i] += fss * (6.0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 + 12.0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 +
                                12.0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 + 24.0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 -
                                12.0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxxy[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 - 24.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
             24.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
             4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 - 8.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxxy[i] += fss * (-12.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 -
                                48.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxxy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxxy[i] += fss * (16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                24.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxxz[i] +=
            fss * (6.0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 + 12.0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 -
                   12.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 - 24.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxxz[i] += fss * (-48.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxxz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                24.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxxz[i] +=
            fss * (16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 - 24.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   48.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxxz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxyy[i] +=
            fss * (18.0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   12.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   12.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxyy[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                12.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxyy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxyy[i] +=
            fss * (8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   36.0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 - 12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxyy[i] += fss * (-24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxyy[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxyy[i] += fss * (16.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * rpa_y * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxyz[i] +=
            fss * (6.0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 + 12.0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 -
                   6.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 - 12.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xxyz[i] +=
            fss * (-4.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   8.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 - 12.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   24.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxyz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxyz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxyz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_xxyz[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxyz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxyz[i] += fss * 16.0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0;

        fints_xxzz[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 + 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
             8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
             8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0);

        fints_xxzz[i] +=
            fss * (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 - 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxzz[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxzz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxzz[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_xxzz[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xxzz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xxzz[i] += fss * 16.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0;

        fints_xyyy[i] += fss * (-3.0 * fe_0 * tbe_0 + 6.0 * fe_0 * fe_0 * tbe_0 * tbe_0 + 3.0 * fe_0 * fe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 - 3.0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyyy[i] +=
            fss * (6.0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 -
                   6.0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0);

        fints_xyyy[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 - 12.0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 +
             3.0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 -
             24.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xyyy[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyyy[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 - 6.0 * rpb_z * rpb_z * tbe_0);

        fints_xyyy[i] += fss * (12.0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 + 6.0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 +
                                24.0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_xyyy[i] +=
            fss * (-12.0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   8.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   12.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   24.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyyy[i] += fss * (-24.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * rpa_y * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                48.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyyy[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyyy[i] += fss * (24.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * rpa_y * rpa_y * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyyz[i] +=
            fss * (6.0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 + 12.0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0);

        fints_xyyz[i] +=
            fss * (-8.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 - 6.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 -
                   12.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 -
                   8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyyz[i] += fss * (-12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyyz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyyz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_xyyz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                12.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyyz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyyz[i] += fss * 16.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0;

        fints_xyzz[i] += fss * (2.0 * fe_0 * fe_0 * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 -
                                2.0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 -
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 - fe_0 * tbe_0);

        fints_xyzz[i] += fss * (fe_0 * fe_0 * tbe_0 * tbe_0 + fe_0 * fe_0 * tbe_0 * tbe_0 + fe_0 * fe_0 * tbe_0 * tbe_0 +
                                2.0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 + 4.0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0);

        fints_xyzz[i] += fss * (4.0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 + 2.0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 -
                                fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 + 2.0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 -
                                fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] +=
            fss * (-fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] +=
            fss *
            (-2.0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 + fe_0 * fe_0 * fe_0 * fe_0 * tbe_0 * tbe_0 * tbe_0 * tbe_0 -
             4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 - 8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 -
             8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] +=
            fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 - 2.0 * rpb_z * rpb_z * tbe_0);

        fints_xyzz[i] += fss * (2.0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 + 2.0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 +
                                2.0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 + 4.0 * rpa_z * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 +
                                4.0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_xyzz[i] +=
            fss * (-2.0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 + 4.0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 -
                   2.0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   4.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] +=
            fss *
            (-4.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
             4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
             4.0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] += fss * (2.0 * fe_0 * fe_0 * fe_0 * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xyzz[i] += fss * 16.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0;

        fints_xzzz[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 + 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
             4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xzzz[i] +=
            fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0);

        fints_xzzz[i] +=
            fss * (12.0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 - 8.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xzzz[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 - 12.0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 -
             4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 - 8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 -
             8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xzzz[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                12.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xzzz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xzzz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xzzz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_xzzz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                12.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_xzzz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_xzzz[i] += fss * 16.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0;

        fints_yyyy[i] +=
            fss * (30.0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 - 20.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   40.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 - 40.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                   6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yyyy[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yyyy[i] +=
            fss * (16.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   60.0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 - 40.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   80.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yyyy[i] += fss * (-80.0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                32.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yyyy[i] += fss * (16.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                24.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                32.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * rpa_y * rpa_y * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yyyz[i] +=
            fss * (6.0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 + 12.0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 -
                   12.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 - 24.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   24.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0);

        fints_yyyz[i] += fss * (-48.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                6.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yyyz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                24.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yyyz[i] +=
            fss * (16.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 - 24.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   48.0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yyyz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                24.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yyzz[i] +=
            fss *
            (-12.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 + 4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
             8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
             8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0);

        fints_yyzz[i] +=
            fss * (-6.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 - 12.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yyzz[i] += fss * (-24.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yyzz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yyzz[i] += fss * (16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_yyzz[i] += fss * (-12.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                4.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yyzz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yyzz[i] += fss * 16.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0;

        fints_yzzz[i] +=
            fss *
            (-4.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 - 4.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
             4.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 + 4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
             4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yzzz[i] +=
            fss * (4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 + 6.0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0);

        fints_yzzz[i] +=
            fss * (12.0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 - 2.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 - 8.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   4.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yzzz[i] +=
            fss *
            (-6.0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 - 12.0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
             4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 - 8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
             8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yzzz[i] += fss * (-8.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                12.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                                2.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yzzz[i] += fss * (8.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yzzz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yzzz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0);

        fints_yzzz[i] += fss * (-4.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                12.0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                8.0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_yzzz[i] += fss * (4.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_yzzz[i] += fss * 16.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0;

        fints_zzzz[i] += fss * (-24.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 +
                                12.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_zzzz[i] += fss * (8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_zzzz[i] += fss * (6.0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 - 12.0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 -
                                48.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                                48.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0);

        fints_zzzz[i] += fss * (6.0 * fe_0 * fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                4.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_zzzz[i] += fss * (12.0 * fe_0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                24.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                8.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_zzzz[i] += fss * (8.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_zzzz[i] +=
            fss * (12.0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 - 24.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 -
                   48.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 +
                   12.0 * fe_0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                   8.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);

        fints_zzzz[i] += fss * (16.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                24.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0 +
                                16.0 * rpa_z * rpa_z * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_z * tbe_0 * tbe_0 * tbe_0 * tbe_0);
    }
}

}  // namespace ovlrec
