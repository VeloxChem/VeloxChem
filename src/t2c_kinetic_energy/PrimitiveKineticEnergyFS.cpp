#include "PrimitiveKineticEnergyFS.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace kinrec {  // kinrec namespace

auto
compPrimitiveKineticEnergyFS(TDoubleArray&       buffer_xxx,
                             TDoubleArray&       buffer_xxy,
                             TDoubleArray&       buffer_xxz,
                             TDoubleArray&       buffer_xyy,
                             TDoubleArray&       buffer_xyz,
                             TDoubleArray&       buffer_xzz,
                             TDoubleArray&       buffer_yyy,
                             TDoubleArray&       buffer_yyz,
                             TDoubleArray&       buffer_yzz,
                             TDoubleArray&       buffer_zzz,
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

#pragma omp simd aligned(fints_xxx,     \
                             fints_xxy, \
                             fints_xxz, \
                             fints_xyy, \
                             fints_xyz, \
                             fints_xzz, \
                             fints_yyy, \
                             fints_yyz, \
                             fints_yzz, \
                             fints_zzz, \
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

        const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

        const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);

        const auto fz_0 = bra_exp * ket_fe[i] * fe_0;

        const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);

        const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;

        const auto rpa_x = -ket_fe[i] * ab_x * fe_0;

        const auto rpa_y = -ket_fe[i] * ab_y * fe_0;

        const auto rpa_z = -ket_fe[i] * ab_z * fe_0;

        const auto fbe_0 = 1.0 / bra_exp;

        fints_xxx[i] += fss * (-3.0 * fbe_0 * rpa_x * fz_0 + 6.0 * fe_0 * rpa_x * fz_0 + 6.0 * rpa_x * rpa_x * rpa_x * fz_0);

        fints_xxx[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_x + rpa_x * rpa_x * rpa_x);

        fints_xxy[i] += fss * (-fbe_0 * rpa_y * fz_0 + 2.0 * fe_0 * rpa_y * fz_0 + 6.0 * rpa_y * rpa_x * rpa_x * fz_0);

        fints_xxy[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_x * rpa_x);

        fints_xxz[i] += fss * (-fbe_0 * rpa_z * fz_0 + 2.0 * fe_0 * rpa_z * fz_0 + 6.0 * rpa_z * rpa_x * rpa_x * fz_0);

        fints_xxz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_x * rpa_x);

        fints_xyy[i] += fss * (-fbe_0 * rpa_x * fz_0 + 2.0 * fe_0 * rpa_x * fz_0 + 6.0 * rpa_y * rpa_y * rpa_x * fz_0);

        fints_xyy[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_y * rpa_y * rpa_x);

        fints_xyz[i] += fss * 6.0 * rpa_z * rpa_y * rpa_x * fz_0;

        fints_xyz[i] += ftt * rpa_z * rpa_y * rpa_x;

        fints_xzz[i] += fss * (-fbe_0 * rpa_x * fz_0 + 2.0 * fe_0 * rpa_x * fz_0 + 6.0 * rpa_z * rpa_z * rpa_x * fz_0);

        fints_xzz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_x + rpa_z * rpa_z * rpa_x);

        fints_yyy[i] += fss * (-3.0 * fbe_0 * rpa_y * fz_0 + 6.0 * fe_0 * rpa_y * fz_0 + 6.0 * rpa_y * rpa_y * rpa_y * fz_0);

        fints_yyy[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_y + rpa_y * rpa_y * rpa_y);

        fints_yyz[i] += fss * (-fbe_0 * rpa_z * fz_0 + 2.0 * fe_0 * rpa_z * fz_0 + 6.0 * rpa_z * rpa_y * rpa_y * fz_0);

        fints_yyz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_y * rpa_y);

        fints_yzz[i] += fss * (-fbe_0 * rpa_y * fz_0 + 2.0 * fe_0 * rpa_y * fz_0 + 6.0 * rpa_z * rpa_z * rpa_y * fz_0);

        fints_yzz[i] += ftt * ((1.0 / 2.0) * fe_0 * rpa_y + rpa_z * rpa_z * rpa_y);

        fints_zzz[i] += fss * (-3.0 * fbe_0 * rpa_z * fz_0 + 6.0 * fe_0 * rpa_z * fz_0 + 6.0 * rpa_z * rpa_z * rpa_z * fz_0);

        fints_zzz[i] += ftt * ((3.0 / 2.0) * fe_0 * rpa_z + rpa_z * rpa_z * rpa_z);
    }
}

}  // namespace kinrec
