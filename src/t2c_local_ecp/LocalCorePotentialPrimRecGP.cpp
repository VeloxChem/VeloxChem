#include "LocalCorePotentialPrimRecGP.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_gp(CSimdArray<double>& pbuffer, 
                                  const size_t idx_gp,
                                  const size_t idx_dp,
                                  const size_t idx_fs,
                                  const size_t idx_fp,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : DP

    auto tg_xx_x = pbuffer.data(idx_dp);

    auto tg_xx_y = pbuffer.data(idx_dp + 1);

    auto tg_xx_z = pbuffer.data(idx_dp + 2);

    auto tg_xy_y = pbuffer.data(idx_dp + 4);

    auto tg_xz_z = pbuffer.data(idx_dp + 8);

    auto tg_yy_x = pbuffer.data(idx_dp + 9);

    auto tg_yy_y = pbuffer.data(idx_dp + 10);

    auto tg_yy_z = pbuffer.data(idx_dp + 11);

    auto tg_yz_z = pbuffer.data(idx_dp + 14);

    auto tg_zz_x = pbuffer.data(idx_dp + 15);

    auto tg_zz_y = pbuffer.data(idx_dp + 16);

    auto tg_zz_z = pbuffer.data(idx_dp + 17);

    // Set up components of auxiliary buffer : FS

    auto tg_xxx_0 = pbuffer.data(idx_fs);

    auto tg_yyy_0 = pbuffer.data(idx_fs + 6);

    auto tg_zzz_0 = pbuffer.data(idx_fs + 9);

    // Set up components of auxiliary buffer : FP

    auto tg_xxx_x = pbuffer.data(idx_fp);

    auto tg_xxx_y = pbuffer.data(idx_fp + 1);

    auto tg_xxx_z = pbuffer.data(idx_fp + 2);

    auto tg_xxy_x = pbuffer.data(idx_fp + 3);

    auto tg_xxy_y = pbuffer.data(idx_fp + 4);

    auto tg_xxz_x = pbuffer.data(idx_fp + 6);

    auto tg_xxz_z = pbuffer.data(idx_fp + 8);

    auto tg_xyy_x = pbuffer.data(idx_fp + 9);

    auto tg_xyy_y = pbuffer.data(idx_fp + 10);

    auto tg_xyy_z = pbuffer.data(idx_fp + 11);

    auto tg_xzz_x = pbuffer.data(idx_fp + 15);

    auto tg_xzz_y = pbuffer.data(idx_fp + 16);

    auto tg_xzz_z = pbuffer.data(idx_fp + 17);

    auto tg_yyy_x = pbuffer.data(idx_fp + 18);

    auto tg_yyy_y = pbuffer.data(idx_fp + 19);

    auto tg_yyy_z = pbuffer.data(idx_fp + 20);

    auto tg_yyz_y = pbuffer.data(idx_fp + 22);

    auto tg_yyz_z = pbuffer.data(idx_fp + 23);

    auto tg_yzz_x = pbuffer.data(idx_fp + 24);

    auto tg_yzz_y = pbuffer.data(idx_fp + 25);

    auto tg_yzz_z = pbuffer.data(idx_fp + 26);

    auto tg_zzz_x = pbuffer.data(idx_fp + 27);

    auto tg_zzz_y = pbuffer.data(idx_fp + 28);

    auto tg_zzz_z = pbuffer.data(idx_fp + 29);

    // Set up components of targeted buffer : GP

    auto tg_xxxx_x = pbuffer.data(idx_gp);

    auto tg_xxxx_y = pbuffer.data(idx_gp + 1);

    auto tg_xxxx_z = pbuffer.data(idx_gp + 2);

    auto tg_xxxy_x = pbuffer.data(idx_gp + 3);

    auto tg_xxxy_y = pbuffer.data(idx_gp + 4);

    auto tg_xxxy_z = pbuffer.data(idx_gp + 5);

    auto tg_xxxz_x = pbuffer.data(idx_gp + 6);

    auto tg_xxxz_y = pbuffer.data(idx_gp + 7);

    auto tg_xxxz_z = pbuffer.data(idx_gp + 8);

    auto tg_xxyy_x = pbuffer.data(idx_gp + 9);

    auto tg_xxyy_y = pbuffer.data(idx_gp + 10);

    auto tg_xxyy_z = pbuffer.data(idx_gp + 11);

    auto tg_xxyz_x = pbuffer.data(idx_gp + 12);

    auto tg_xxyz_y = pbuffer.data(idx_gp + 13);

    auto tg_xxyz_z = pbuffer.data(idx_gp + 14);

    auto tg_xxzz_x = pbuffer.data(idx_gp + 15);

    auto tg_xxzz_y = pbuffer.data(idx_gp + 16);

    auto tg_xxzz_z = pbuffer.data(idx_gp + 17);

    auto tg_xyyy_x = pbuffer.data(idx_gp + 18);

    auto tg_xyyy_y = pbuffer.data(idx_gp + 19);

    auto tg_xyyy_z = pbuffer.data(idx_gp + 20);

    auto tg_xyyz_x = pbuffer.data(idx_gp + 21);

    auto tg_xyyz_y = pbuffer.data(idx_gp + 22);

    auto tg_xyyz_z = pbuffer.data(idx_gp + 23);

    auto tg_xyzz_x = pbuffer.data(idx_gp + 24);

    auto tg_xyzz_y = pbuffer.data(idx_gp + 25);

    auto tg_xyzz_z = pbuffer.data(idx_gp + 26);

    auto tg_xzzz_x = pbuffer.data(idx_gp + 27);

    auto tg_xzzz_y = pbuffer.data(idx_gp + 28);

    auto tg_xzzz_z = pbuffer.data(idx_gp + 29);

    auto tg_yyyy_x = pbuffer.data(idx_gp + 30);

    auto tg_yyyy_y = pbuffer.data(idx_gp + 31);

    auto tg_yyyy_z = pbuffer.data(idx_gp + 32);

    auto tg_yyyz_x = pbuffer.data(idx_gp + 33);

    auto tg_yyyz_y = pbuffer.data(idx_gp + 34);

    auto tg_yyyz_z = pbuffer.data(idx_gp + 35);

    auto tg_yyzz_x = pbuffer.data(idx_gp + 36);

    auto tg_yyzz_y = pbuffer.data(idx_gp + 37);

    auto tg_yyzz_z = pbuffer.data(idx_gp + 38);

    auto tg_yzzz_x = pbuffer.data(idx_gp + 39);

    auto tg_yzzz_y = pbuffer.data(idx_gp + 40);

    auto tg_yzzz_z = pbuffer.data(idx_gp + 41);

    auto tg_zzzz_x = pbuffer.data(idx_gp + 42);

    auto tg_zzzz_y = pbuffer.data(idx_gp + 43);

    auto tg_zzzz_z = pbuffer.data(idx_gp + 44);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xx_x, tg_xx_y, tg_xx_z, tg_xxx_0, tg_xxx_x, tg_xxx_y, tg_xxx_z, tg_xxxx_x, tg_xxxx_y, tg_xxxx_z, tg_xxxy_x, tg_xxxy_y, tg_xxxy_z, tg_xxxz_x, tg_xxxz_y, tg_xxxz_z, tg_xxy_x, tg_xxy_y, tg_xxyy_x, tg_xxyy_y, tg_xxyy_z, tg_xxyz_x, tg_xxyz_y, tg_xxyz_z, tg_xxz_x, tg_xxz_z, tg_xxzz_x, tg_xxzz_y, tg_xxzz_z, tg_xy_y, tg_xyy_x, tg_xyy_y, tg_xyy_z, tg_xyyy_x, tg_xyyy_y, tg_xyyy_z, tg_xyyz_x, tg_xyyz_y, tg_xyyz_z, tg_xyzz_x, tg_xyzz_y, tg_xyzz_z, tg_xz_z, tg_xzz_x, tg_xzz_y, tg_xzz_z, tg_xzzz_x, tg_xzzz_y, tg_xzzz_z, tg_yy_x, tg_yy_y, tg_yy_z, tg_yyy_0, tg_yyy_x, tg_yyy_y, tg_yyy_z, tg_yyyy_x, tg_yyyy_y, tg_yyyy_z, tg_yyyz_x, tg_yyyz_y, tg_yyyz_z, tg_yyz_y, tg_yyz_z, tg_yyzz_x, tg_yyzz_y, tg_yyzz_z, tg_yz_z, tg_yzz_x, tg_yzz_y, tg_yzz_z, tg_yzzz_x, tg_yzzz_y, tg_yzzz_z, tg_zz_x, tg_zz_y, tg_zz_z, tg_zzz_0, tg_zzz_x, tg_zzz_y, tg_zzz_z, tg_zzzz_x, tg_zzzz_y, tg_zzzz_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxx_x[i] = 3.0 * tg_xx_x[i] * fxi[i] + tg_xxx_0[i] * fxi[i] + tg_xxx_x[i] * ra_x[i];

        tg_xxxx_y[i] = 3.0 * tg_xx_y[i] * fxi[i] + tg_xxx_y[i] * ra_x[i];

        tg_xxxx_z[i] = 3.0 * tg_xx_z[i] * fxi[i] + tg_xxx_z[i] * ra_x[i];

        tg_xxxy_x[i] = tg_xxx_x[i] * ra_y[i];

        tg_xxxy_y[i] = 2.0 * tg_xy_y[i] * fxi[i] + tg_xxy_y[i] * ra_x[i];

        tg_xxxy_z[i] = tg_xxx_z[i] * ra_y[i];

        tg_xxxz_x[i] = tg_xxx_x[i] * ra_z[i];

        tg_xxxz_y[i] = tg_xxx_y[i] * ra_z[i];

        tg_xxxz_z[i] = 2.0 * tg_xz_z[i] * fxi[i] + tg_xxz_z[i] * ra_x[i];

        tg_xxyy_x[i] = tg_xx_x[i] * fxi[i] + tg_xxy_x[i] * ra_y[i];

        tg_xxyy_y[i] = tg_yy_y[i] * fxi[i] + tg_xyy_y[i] * ra_x[i];

        tg_xxyy_z[i] = tg_yy_z[i] * fxi[i] + tg_xyy_z[i] * ra_x[i];

        tg_xxyz_x[i] = tg_xxz_x[i] * ra_y[i];

        tg_xxyz_y[i] = tg_xxy_y[i] * ra_z[i];

        tg_xxyz_z[i] = tg_xxz_z[i] * ra_y[i];

        tg_xxzz_x[i] = tg_xx_x[i] * fxi[i] + tg_xxz_x[i] * ra_z[i];

        tg_xxzz_y[i] = tg_zz_y[i] * fxi[i] + tg_xzz_y[i] * ra_x[i];

        tg_xxzz_z[i] = tg_zz_z[i] * fxi[i] + tg_xzz_z[i] * ra_x[i];

        tg_xyyy_x[i] = tg_yyy_0[i] * fxi[i] + tg_yyy_x[i] * ra_x[i];

        tg_xyyy_y[i] = tg_yyy_y[i] * ra_x[i];

        tg_xyyy_z[i] = tg_yyy_z[i] * ra_x[i];

        tg_xyyz_x[i] = tg_xyy_x[i] * ra_z[i];

        tg_xyyz_y[i] = tg_yyz_y[i] * ra_x[i];

        tg_xyyz_z[i] = tg_yyz_z[i] * ra_x[i];

        tg_xyzz_x[i] = tg_xzz_x[i] * ra_y[i];

        tg_xyzz_y[i] = tg_yzz_y[i] * ra_x[i];

        tg_xyzz_z[i] = tg_yzz_z[i] * ra_x[i];

        tg_xzzz_x[i] = tg_zzz_0[i] * fxi[i] + tg_zzz_x[i] * ra_x[i];

        tg_xzzz_y[i] = tg_zzz_y[i] * ra_x[i];

        tg_xzzz_z[i] = tg_zzz_z[i] * ra_x[i];

        tg_yyyy_x[i] = 3.0 * tg_yy_x[i] * fxi[i] + tg_yyy_x[i] * ra_y[i];

        tg_yyyy_y[i] = 3.0 * tg_yy_y[i] * fxi[i] + tg_yyy_0[i] * fxi[i] + tg_yyy_y[i] * ra_y[i];

        tg_yyyy_z[i] = 3.0 * tg_yy_z[i] * fxi[i] + tg_yyy_z[i] * ra_y[i];

        tg_yyyz_x[i] = tg_yyy_x[i] * ra_z[i];

        tg_yyyz_y[i] = tg_yyy_y[i] * ra_z[i];

        tg_yyyz_z[i] = 2.0 * tg_yz_z[i] * fxi[i] + tg_yyz_z[i] * ra_y[i];

        tg_yyzz_x[i] = tg_zz_x[i] * fxi[i] + tg_yzz_x[i] * ra_y[i];

        tg_yyzz_y[i] = tg_yy_y[i] * fxi[i] + tg_yyz_y[i] * ra_z[i];

        tg_yyzz_z[i] = tg_zz_z[i] * fxi[i] + tg_yzz_z[i] * ra_y[i];

        tg_yzzz_x[i] = tg_zzz_x[i] * ra_y[i];

        tg_yzzz_y[i] = tg_zzz_0[i] * fxi[i] + tg_zzz_y[i] * ra_y[i];

        tg_yzzz_z[i] = tg_zzz_z[i] * ra_y[i];

        tg_zzzz_x[i] = 3.0 * tg_zz_x[i] * fxi[i] + tg_zzz_x[i] * ra_z[i];

        tg_zzzz_y[i] = 3.0 * tg_zz_y[i] * fxi[i] + tg_zzz_y[i] * ra_z[i];

        tg_zzzz_z[i] = 3.0 * tg_zz_z[i] * fxi[i] + tg_zzz_0[i] * fxi[i] + tg_zzz_z[i] * ra_z[i];
    }
}

} // t2lecp namespace

