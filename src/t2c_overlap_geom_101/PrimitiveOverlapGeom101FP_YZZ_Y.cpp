#include "PrimitiveOverlapGeom101FP_YZZ_Y.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ovlrec {  // ovlrec namespace

auto
compPrimitiveOverlapGeom101FP_YZZ_Y(TDoubleArray&       buffer_x_x,
                                    TDoubleArray&       buffer_x_y,
                                    TDoubleArray&       buffer_x_z,
                                    TDoubleArray&       buffer_y_x,
                                    TDoubleArray&       buffer_y_y,
                                    TDoubleArray&       buffer_y_z,
                                    TDoubleArray&       buffer_z_x,
                                    TDoubleArray&       buffer_z_y,
                                    TDoubleArray&       buffer_z_z,
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

    auto fints_x_x = buffer_x_x.data();

    auto fints_x_y = buffer_x_y.data();

    auto fints_x_z = buffer_x_z.data();

    auto fints_y_x = buffer_y_x.data();

    auto fints_y_y = buffer_y_y.data();

    auto fints_y_z = buffer_y_z.data();

    auto fints_z_x = buffer_z_x.data();

    auto fints_z_y = buffer_z_y.data();

    auto fints_z_z = buffer_z_z.data();

#pragma omp simd aligned(fints_x_x,     \
                             fints_x_y, \
                             fints_x_z, \
                             fints_y_x, \
                             fints_y_y, \
                             fints_y_z, \
                             fints_z_x, \
                             fints_z_y, \
                             fints_z_z, \
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

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto rpb_x = bra_exp * ab_x * fe_0;

        const auto rpb_y = bra_exp * ab_y * fe_0;

        const auto rpb_z = bra_exp * ab_z * fe_0;

        const auto tke_0 = ket_fe[i];

        const auto tbe_0 = bra_exp;

        fints_x_x[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tke_0 +
                               fe_0 * fe_0 * rpa_x * rpb_x * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tke_0 +
                               2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_x * tbe_0 * tke_0);

        fints_x_x[i] +=
            fss * (2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0 +
                   4.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_x_y[i] +=
            fss * (-fe_0 * rpa_y * rpa_x * tbe_0 - 2.0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 + fe_0 * fe_0 * rpa_y * rpa_x * tbe_0 * tke_0 +
                   2.0 * fe_0 * fe_0 * rpa_x * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_x * tbe_0 * tke_0);

        fints_x_y[i] +=
            fss * (4.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0 +
                   4.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_x_z[i] += fss * (fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_z * rpa_x * tbe_0 * tke_0 +
                               fe_0 * fe_0 * rpa_x * rpb_z * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0 +
                               2.0 * fe_0 * rpa_z * rpa_y * rpa_x * rpb_y * tbe_0 * tke_0);

        fints_x_z[i] +=
            fss * (2.0 * fe_0 * rpa_z * rpa_z * rpa_x * rpb_z * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0 +
                   4.0 * rpa_z * rpa_z * rpa_y * rpa_x * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_y_x[i] += fss * (fe_0 * fe_0 * rpa_y * rpb_x * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_y * rpb_x * tbe_0 * tke_0 +
                               2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * tbe_0 * tke_0 +
                               2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_x * tbe_0 * tke_0 - fe_0 * rpb_y * rpb_x * tke_0);

        fints_y_x[i] +=
            fss * (-2.0 * rpa_z * rpa_z * rpb_y * rpb_x * tke_0 + fe_0 * fe_0 * rpb_y * rpb_x * tbe_0 * tke_0 +
                   2.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_x * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 +
                   4.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_y_y[i] += fss * ((1.0 / 2.0) * fe_0 + rpa_z * rpa_z - (1.0 / 2.0) * fe_0 * fe_0 * tbe_0 - fe_0 * rpa_y * rpa_y * tbe_0 -
                               fe_0 * rpa_z * rpa_z * tbe_0);

        fints_y_y[i] += fss * (fe_0 * fe_0 * fe_0 * tbe_0 * tke_0 - 2.0 * rpa_z * rpa_z * rpa_y * rpa_y * tbe_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tke_0 - (1.0 / 2.0) * fe_0 * fe_0 * tke_0 - fe_0 * rpa_z * rpa_z * tke_0);

        fints_y_y[i] += fss * ((1.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tke_0 +
                               fe_0 * fe_0 * rpa_y * rpa_y * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tke_0 +
                               2.0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tke_0);

        fints_y_y[i] +=
            fss *
            (2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpa_y * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tke_0 +
             4.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tke_0 - fe_0 * rpb_y * rpb_y * tke_0 - 2.0 * rpa_z * rpa_z * rpb_y * rpb_y * tke_0);

        fints_y_y[i] +=
            fss * (fe_0 * fe_0 * rpb_y * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpa_z * rpb_y * rpb_y * tbe_0 * tke_0 +
                   2.0 * fe_0 * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0 + 4.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_y_z[i] +=
            fss * (fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tke_0 +
                   fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tke_0 - fe_0 * rpa_z * rpb_y * tke_0);

        fints_y_z[i] +=
            fss * (-fe_0 * rpa_z * rpb_y * tke_0 + fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 +
                   fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_y * rpb_z * tbe_0 * tke_0);

        fints_y_z[i] +=
            fss * (2.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpa_y * rpa_y * rpb_y * tbe_0 * tke_0 +
                   2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_z * tbe_0 * tke_0 -
                   fe_0 * rpb_z * rpb_y * tke_0);

        fints_y_z[i] +=
            fss * (-2.0 * rpa_z * rpa_z * rpb_z * rpb_y * tke_0 + fe_0 * fe_0 * rpb_z * rpb_y * tbe_0 * tke_0 +
                   2.0 * fe_0 * rpa_y * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpa_z * rpb_z * rpb_y * tbe_0 * tke_0 +
                   4.0 * rpa_z * rpa_z * rpa_y * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0);

        fints_z_x[i] += fss * (-2.0 * fe_0 * rpa_z * rpb_x * tke_0 + fe_0 * fe_0 * rpa_z * rpb_x * tbe_0 * tke_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpb_x * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_x * tbe_0 * tke_0 -
                               4.0 * rpa_z * rpa_y * rpb_y * rpb_x * tke_0);

        fints_z_x[i] +=
            fss * (2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0 +
                   4.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_x * tbe_0 * tke_0);

        fints_z_y[i] += fss * (2.0 * rpa_z * rpa_y - fe_0 * rpa_z * rpa_y * tbe_0 - 2.0 * fe_0 * rpa_z * rpa_y * tbe_0 -
                               2.0 * rpa_z * rpa_z * rpa_z * rpa_y * tbe_0 - 2.0 * fe_0 * rpa_z * rpa_y * tke_0);

        fints_z_y[i] += fss * (-4.0 * fe_0 * rpa_z * rpb_y * tke_0 + fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tke_0 +
                               2.0 * fe_0 * fe_0 * rpa_z * rpa_y * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0 +
                               4.0 * fe_0 * fe_0 * rpa_z * rpb_y * tbe_0 * tke_0);

        fints_z_y[i] +=
            fss * (2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpa_y * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_y * tbe_0 * tke_0 -
                   4.0 * rpa_z * rpa_y * rpb_y * rpb_y * tke_0 + 2.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0 +
                   4.0 * fe_0 * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0);

        fints_z_y[i] += fss * 4.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_y * rpb_y * tbe_0 * tke_0;

        fints_z_z[i] += fss * (-fe_0 * fe_0 * tke_0 + fe_0 * fe_0 * fe_0 * tbe_0 * tke_0 + (1.0 / 2.0) * fe_0 * fe_0 * fe_0 * tbe_0 * tke_0 +
                               fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tke_0);

        fints_z_z[i] += fss * (fe_0 * fe_0 * rpa_z * rpa_z * tbe_0 * tke_0 - 2.0 * fe_0 * rpa_y * rpb_y * tke_0 - 2.0 * fe_0 * rpa_z * rpb_z * tke_0 +
                               2.0 * fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tke_0 + fe_0 * fe_0 * rpa_y * rpb_y * tbe_0 * tke_0);

        fints_z_z[i] +=
            fss * (fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tke_0 + 2.0 * fe_0 * fe_0 * rpa_z * rpb_z * tbe_0 * tke_0 +
                   2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tke_0 + 2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tke_0 +
                   2.0 * fe_0 * rpa_z * rpa_z * rpa_y * rpb_y * tbe_0 * tke_0);

        fints_z_z[i] +=
            fss * (2.0 * fe_0 * rpa_z * rpa_z * rpa_z * rpb_z * tbe_0 * tke_0 - 4.0 * rpa_z * rpa_y * rpb_z * rpb_y * tke_0 +
                   2.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0 + 4.0 * fe_0 * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0 +
                   4.0 * rpa_z * rpa_z * rpa_z * rpa_y * rpb_z * rpb_y * tbe_0 * tke_0);
    }
}

}  // namespace ovlrec
