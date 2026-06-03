#include "GeometricalDerivatives110ForDP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_dp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_dp,
                         const int idx_op_sp,
                         const int idx_op_ps,
                         const int idx_op_pd,
                         const int idx_op_dp,
                         const int idx_op_fs,
                         const int idx_op_fd,
                         const int idx_op_gp,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SP

    auto tr_0_x = pbuffer.data(idx_op_sp);

    auto tr_0_y = pbuffer.data(idx_op_sp + 1);

    auto tr_0_z = pbuffer.data(idx_op_sp + 2);

    // Set up components of auxiliary buffer : PS

    auto tr_x_0 = pbuffer.data(idx_op_ps);

    auto tr_y_0 = pbuffer.data(idx_op_ps + 1);

    auto tr_z_0 = pbuffer.data(idx_op_ps + 2);

    // Set up components of auxiliary buffer : PD

    auto tr_x_xx = pbuffer.data(idx_op_pd);

    auto tr_x_xy = pbuffer.data(idx_op_pd + 1);

    auto tr_x_xz = pbuffer.data(idx_op_pd + 2);

    auto tr_x_yy = pbuffer.data(idx_op_pd + 3);

    auto tr_x_yz = pbuffer.data(idx_op_pd + 4);

    auto tr_x_zz = pbuffer.data(idx_op_pd + 5);

    auto tr_y_xx = pbuffer.data(idx_op_pd + 6);

    auto tr_y_xy = pbuffer.data(idx_op_pd + 7);

    auto tr_y_xz = pbuffer.data(idx_op_pd + 8);

    auto tr_y_yy = pbuffer.data(idx_op_pd + 9);

    auto tr_y_yz = pbuffer.data(idx_op_pd + 10);

    auto tr_y_zz = pbuffer.data(idx_op_pd + 11);

    auto tr_z_xx = pbuffer.data(idx_op_pd + 12);

    auto tr_z_xy = pbuffer.data(idx_op_pd + 13);

    auto tr_z_xz = pbuffer.data(idx_op_pd + 14);

    auto tr_z_yy = pbuffer.data(idx_op_pd + 15);

    auto tr_z_yz = pbuffer.data(idx_op_pd + 16);

    auto tr_z_zz = pbuffer.data(idx_op_pd + 17);

    // Set up components of auxiliary buffer : DP

    auto tr_xx_x = pbuffer.data(idx_op_dp);

    auto tr_xx_y = pbuffer.data(idx_op_dp + 1);

    auto tr_xx_z = pbuffer.data(idx_op_dp + 2);

    auto tr_xy_x = pbuffer.data(idx_op_dp + 3);

    auto tr_xy_y = pbuffer.data(idx_op_dp + 4);

    auto tr_xy_z = pbuffer.data(idx_op_dp + 5);

    auto tr_xz_x = pbuffer.data(idx_op_dp + 6);

    auto tr_xz_y = pbuffer.data(idx_op_dp + 7);

    auto tr_xz_z = pbuffer.data(idx_op_dp + 8);

    auto tr_yy_x = pbuffer.data(idx_op_dp + 9);

    auto tr_yy_y = pbuffer.data(idx_op_dp + 10);

    auto tr_yy_z = pbuffer.data(idx_op_dp + 11);

    auto tr_yz_x = pbuffer.data(idx_op_dp + 12);

    auto tr_yz_y = pbuffer.data(idx_op_dp + 13);

    auto tr_yz_z = pbuffer.data(idx_op_dp + 14);

    auto tr_zz_x = pbuffer.data(idx_op_dp + 15);

    auto tr_zz_y = pbuffer.data(idx_op_dp + 16);

    auto tr_zz_z = pbuffer.data(idx_op_dp + 17);

    // Set up components of auxiliary buffer : FS

    auto tr_xxx_0 = pbuffer.data(idx_op_fs);

    auto tr_xxy_0 = pbuffer.data(idx_op_fs + 1);

    auto tr_xxz_0 = pbuffer.data(idx_op_fs + 2);

    auto tr_xyy_0 = pbuffer.data(idx_op_fs + 3);

    auto tr_xyz_0 = pbuffer.data(idx_op_fs + 4);

    auto tr_xzz_0 = pbuffer.data(idx_op_fs + 5);

    auto tr_yyy_0 = pbuffer.data(idx_op_fs + 6);

    auto tr_yyz_0 = pbuffer.data(idx_op_fs + 7);

    auto tr_yzz_0 = pbuffer.data(idx_op_fs + 8);

    auto tr_zzz_0 = pbuffer.data(idx_op_fs + 9);

    // Set up components of auxiliary buffer : FD

    auto tr_xxx_xx = pbuffer.data(idx_op_fd);

    auto tr_xxx_xy = pbuffer.data(idx_op_fd + 1);

    auto tr_xxx_xz = pbuffer.data(idx_op_fd + 2);

    auto tr_xxx_yy = pbuffer.data(idx_op_fd + 3);

    auto tr_xxx_yz = pbuffer.data(idx_op_fd + 4);

    auto tr_xxx_zz = pbuffer.data(idx_op_fd + 5);

    auto tr_xxy_xx = pbuffer.data(idx_op_fd + 6);

    auto tr_xxy_xy = pbuffer.data(idx_op_fd + 7);

    auto tr_xxy_xz = pbuffer.data(idx_op_fd + 8);

    auto tr_xxy_yy = pbuffer.data(idx_op_fd + 9);

    auto tr_xxy_yz = pbuffer.data(idx_op_fd + 10);

    auto tr_xxy_zz = pbuffer.data(idx_op_fd + 11);

    auto tr_xxz_xx = pbuffer.data(idx_op_fd + 12);

    auto tr_xxz_xy = pbuffer.data(idx_op_fd + 13);

    auto tr_xxz_xz = pbuffer.data(idx_op_fd + 14);

    auto tr_xxz_yy = pbuffer.data(idx_op_fd + 15);

    auto tr_xxz_yz = pbuffer.data(idx_op_fd + 16);

    auto tr_xxz_zz = pbuffer.data(idx_op_fd + 17);

    auto tr_xyy_xx = pbuffer.data(idx_op_fd + 18);

    auto tr_xyy_xy = pbuffer.data(idx_op_fd + 19);

    auto tr_xyy_xz = pbuffer.data(idx_op_fd + 20);

    auto tr_xyy_yy = pbuffer.data(idx_op_fd + 21);

    auto tr_xyy_yz = pbuffer.data(idx_op_fd + 22);

    auto tr_xyy_zz = pbuffer.data(idx_op_fd + 23);

    auto tr_xyz_xx = pbuffer.data(idx_op_fd + 24);

    auto tr_xyz_xy = pbuffer.data(idx_op_fd + 25);

    auto tr_xyz_xz = pbuffer.data(idx_op_fd + 26);

    auto tr_xyz_yy = pbuffer.data(idx_op_fd + 27);

    auto tr_xyz_yz = pbuffer.data(idx_op_fd + 28);

    auto tr_xyz_zz = pbuffer.data(idx_op_fd + 29);

    auto tr_xzz_xx = pbuffer.data(idx_op_fd + 30);

    auto tr_xzz_xy = pbuffer.data(idx_op_fd + 31);

    auto tr_xzz_xz = pbuffer.data(idx_op_fd + 32);

    auto tr_xzz_yy = pbuffer.data(idx_op_fd + 33);

    auto tr_xzz_yz = pbuffer.data(idx_op_fd + 34);

    auto tr_xzz_zz = pbuffer.data(idx_op_fd + 35);

    auto tr_yyy_xx = pbuffer.data(idx_op_fd + 36);

    auto tr_yyy_xy = pbuffer.data(idx_op_fd + 37);

    auto tr_yyy_xz = pbuffer.data(idx_op_fd + 38);

    auto tr_yyy_yy = pbuffer.data(idx_op_fd + 39);

    auto tr_yyy_yz = pbuffer.data(idx_op_fd + 40);

    auto tr_yyy_zz = pbuffer.data(idx_op_fd + 41);

    auto tr_yyz_xx = pbuffer.data(idx_op_fd + 42);

    auto tr_yyz_xy = pbuffer.data(idx_op_fd + 43);

    auto tr_yyz_xz = pbuffer.data(idx_op_fd + 44);

    auto tr_yyz_yy = pbuffer.data(idx_op_fd + 45);

    auto tr_yyz_yz = pbuffer.data(idx_op_fd + 46);

    auto tr_yyz_zz = pbuffer.data(idx_op_fd + 47);

    auto tr_yzz_xx = pbuffer.data(idx_op_fd + 48);

    auto tr_yzz_xy = pbuffer.data(idx_op_fd + 49);

    auto tr_yzz_xz = pbuffer.data(idx_op_fd + 50);

    auto tr_yzz_yy = pbuffer.data(idx_op_fd + 51);

    auto tr_yzz_yz = pbuffer.data(idx_op_fd + 52);

    auto tr_yzz_zz = pbuffer.data(idx_op_fd + 53);

    auto tr_zzz_xx = pbuffer.data(idx_op_fd + 54);

    auto tr_zzz_xy = pbuffer.data(idx_op_fd + 55);

    auto tr_zzz_xz = pbuffer.data(idx_op_fd + 56);

    auto tr_zzz_yy = pbuffer.data(idx_op_fd + 57);

    auto tr_zzz_yz = pbuffer.data(idx_op_fd + 58);

    auto tr_zzz_zz = pbuffer.data(idx_op_fd + 59);

    // Set up components of auxiliary buffer : GP

    auto tr_xxxx_x = pbuffer.data(idx_op_gp);

    auto tr_xxxx_y = pbuffer.data(idx_op_gp + 1);

    auto tr_xxxx_z = pbuffer.data(idx_op_gp + 2);

    auto tr_xxxy_x = pbuffer.data(idx_op_gp + 3);

    auto tr_xxxy_y = pbuffer.data(idx_op_gp + 4);

    auto tr_xxxy_z = pbuffer.data(idx_op_gp + 5);

    auto tr_xxxz_x = pbuffer.data(idx_op_gp + 6);

    auto tr_xxxz_y = pbuffer.data(idx_op_gp + 7);

    auto tr_xxxz_z = pbuffer.data(idx_op_gp + 8);

    auto tr_xxyy_x = pbuffer.data(idx_op_gp + 9);

    auto tr_xxyy_y = pbuffer.data(idx_op_gp + 10);

    auto tr_xxyy_z = pbuffer.data(idx_op_gp + 11);

    auto tr_xxyz_x = pbuffer.data(idx_op_gp + 12);

    auto tr_xxyz_y = pbuffer.data(idx_op_gp + 13);

    auto tr_xxyz_z = pbuffer.data(idx_op_gp + 14);

    auto tr_xxzz_x = pbuffer.data(idx_op_gp + 15);

    auto tr_xxzz_y = pbuffer.data(idx_op_gp + 16);

    auto tr_xxzz_z = pbuffer.data(idx_op_gp + 17);

    auto tr_xyyy_x = pbuffer.data(idx_op_gp + 18);

    auto tr_xyyy_y = pbuffer.data(idx_op_gp + 19);

    auto tr_xyyy_z = pbuffer.data(idx_op_gp + 20);

    auto tr_xyyz_x = pbuffer.data(idx_op_gp + 21);

    auto tr_xyyz_y = pbuffer.data(idx_op_gp + 22);

    auto tr_xyyz_z = pbuffer.data(idx_op_gp + 23);

    auto tr_xyzz_x = pbuffer.data(idx_op_gp + 24);

    auto tr_xyzz_y = pbuffer.data(idx_op_gp + 25);

    auto tr_xyzz_z = pbuffer.data(idx_op_gp + 26);

    auto tr_xzzz_x = pbuffer.data(idx_op_gp + 27);

    auto tr_xzzz_y = pbuffer.data(idx_op_gp + 28);

    auto tr_xzzz_z = pbuffer.data(idx_op_gp + 29);

    auto tr_yyyy_x = pbuffer.data(idx_op_gp + 30);

    auto tr_yyyy_y = pbuffer.data(idx_op_gp + 31);

    auto tr_yyyy_z = pbuffer.data(idx_op_gp + 32);

    auto tr_yyyz_x = pbuffer.data(idx_op_gp + 33);

    auto tr_yyyz_y = pbuffer.data(idx_op_gp + 34);

    auto tr_yyyz_z = pbuffer.data(idx_op_gp + 35);

    auto tr_yyzz_x = pbuffer.data(idx_op_gp + 36);

    auto tr_yyzz_y = pbuffer.data(idx_op_gp + 37);

    auto tr_yyzz_z = pbuffer.data(idx_op_gp + 38);

    auto tr_yzzz_x = pbuffer.data(idx_op_gp + 39);

    auto tr_yzzz_y = pbuffer.data(idx_op_gp + 40);

    auto tr_yzzz_z = pbuffer.data(idx_op_gp + 41);

    auto tr_zzzz_x = pbuffer.data(idx_op_gp + 42);

    auto tr_zzzz_y = pbuffer.data(idx_op_gp + 43);

    auto tr_zzzz_z = pbuffer.data(idx_op_gp + 44);

    // Set up 0-3 components of targeted buffer : DP

    auto tr_x_0_x_xx_x = pbuffer.data(idx_op_geom_110_dp);

    auto tr_x_0_x_xx_y = pbuffer.data(idx_op_geom_110_dp + 1);

    auto tr_x_0_x_xx_z = pbuffer.data(idx_op_geom_110_dp + 2);

    #pragma omp simd aligned(tr_0_x, tr_0_y, tr_0_z, tr_x_0, tr_x_0_x_xx_x, tr_x_0_x_xx_y, tr_x_0_x_xx_z, tr_x_xx, tr_x_xy, tr_x_xz, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxx_0, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxxx_x, tr_xxxx_y, tr_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xx_x[i] = 2.0 * tr_0_x[i] + 2.0 * tr_x_0[i] - 4.0 * tr_x_xx[i] * tke_0 - 10.0 * tr_xx_x[i] * tbe_0 - 2.0 * tr_xxx_0[i] * tbe_0 + 4.0 * tr_xxx_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_y[i] = 2.0 * tr_0_y[i] - 4.0 * tr_x_xy[i] * tke_0 - 10.0 * tr_xx_y[i] * tbe_0 + 4.0 * tr_xxx_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_z[i] = 2.0 * tr_0_z[i] - 4.0 * tr_x_xz[i] * tke_0 - 10.0 * tr_xx_z[i] * tbe_0 + 4.0 * tr_xxx_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_z[i] * tbe_0 * tbe_0;
    }

    // Set up 3-6 components of targeted buffer : DP

    auto tr_x_0_x_xy_x = pbuffer.data(idx_op_geom_110_dp + 3);

    auto tr_x_0_x_xy_y = pbuffer.data(idx_op_geom_110_dp + 4);

    auto tr_x_0_x_xy_z = pbuffer.data(idx_op_geom_110_dp + 5);

    #pragma omp simd aligned(tr_x_0_x_xy_x, tr_x_0_x_xy_y, tr_x_0_x_xy_z, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxy_0, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xy_x, tr_xy_y, tr_xy_z, tr_y_0, tr_y_xx, tr_y_xy, tr_y_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xy_x[i] = tr_y_0[i] - 2.0 * tr_y_xx[i] * tke_0 - 6.0 * tr_xy_x[i] * tbe_0 - 2.0 * tr_xxy_0[i] * tbe_0 + 4.0 * tr_xxy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_y[i] = -2.0 * tr_y_xy[i] * tke_0 - 6.0 * tr_xy_y[i] * tbe_0 + 4.0 * tr_xxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_z[i] = -2.0 * tr_y_xz[i] * tke_0 - 6.0 * tr_xy_z[i] * tbe_0 + 4.0 * tr_xxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 6-9 components of targeted buffer : DP

    auto tr_x_0_x_xz_x = pbuffer.data(idx_op_geom_110_dp + 6);

    auto tr_x_0_x_xz_y = pbuffer.data(idx_op_geom_110_dp + 7);

    auto tr_x_0_x_xz_z = pbuffer.data(idx_op_geom_110_dp + 8);

    #pragma omp simd aligned(tr_x_0_x_xz_x, tr_x_0_x_xz_y, tr_x_0_x_xz_z, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xxz_0, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xz_x, tr_xz_y, tr_xz_z, tr_z_0, tr_z_xx, tr_z_xy, tr_z_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xz_x[i] = tr_z_0[i] - 2.0 * tr_z_xx[i] * tke_0 - 6.0 * tr_xz_x[i] * tbe_0 - 2.0 * tr_xxz_0[i] * tbe_0 + 4.0 * tr_xxz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_y[i] = -2.0 * tr_z_xy[i] * tke_0 - 6.0 * tr_xz_y[i] * tbe_0 + 4.0 * tr_xxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_z[i] = -2.0 * tr_z_xz[i] * tke_0 - 6.0 * tr_xz_z[i] * tbe_0 + 4.0 * tr_xxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 9-12 components of targeted buffer : DP

    auto tr_x_0_x_yy_x = pbuffer.data(idx_op_geom_110_dp + 9);

    auto tr_x_0_x_yy_y = pbuffer.data(idx_op_geom_110_dp + 10);

    auto tr_x_0_x_yy_z = pbuffer.data(idx_op_geom_110_dp + 11);

    #pragma omp simd aligned(tr_x_0_x_yy_x, tr_x_0_x_yy_y, tr_x_0_x_yy_z, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xyy_0, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_yy_x, tr_yy_y, tr_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yy_x[i] = -2.0 * tr_yy_x[i] * tbe_0 - 2.0 * tr_xyy_0[i] * tbe_0 + 4.0 * tr_xyy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_y[i] = -2.0 * tr_yy_y[i] * tbe_0 + 4.0 * tr_xyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_z[i] = -2.0 * tr_yy_z[i] * tbe_0 + 4.0 * tr_xyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 12-15 components of targeted buffer : DP

    auto tr_x_0_x_yz_x = pbuffer.data(idx_op_geom_110_dp + 12);

    auto tr_x_0_x_yz_y = pbuffer.data(idx_op_geom_110_dp + 13);

    auto tr_x_0_x_yz_z = pbuffer.data(idx_op_geom_110_dp + 14);

    #pragma omp simd aligned(tr_x_0_x_yz_x, tr_x_0_x_yz_y, tr_x_0_x_yz_z, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xyz_0, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_yz_x, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yz_x[i] = -2.0 * tr_yz_x[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_y[i] = -2.0 * tr_yz_y[i] * tbe_0 + 4.0 * tr_xyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_z[i] = -2.0 * tr_yz_z[i] * tbe_0 + 4.0 * tr_xyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 15-18 components of targeted buffer : DP

    auto tr_x_0_x_zz_x = pbuffer.data(idx_op_geom_110_dp + 15);

    auto tr_x_0_x_zz_y = pbuffer.data(idx_op_geom_110_dp + 16);

    auto tr_x_0_x_zz_z = pbuffer.data(idx_op_geom_110_dp + 17);

    #pragma omp simd aligned(tr_x_0_x_zz_x, tr_x_0_x_zz_y, tr_x_0_x_zz_z, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_xzz_0, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_zz_x[i] = -2.0 * tr_zz_x[i] * tbe_0 - 2.0 * tr_xzz_0[i] * tbe_0 + 4.0 * tr_xzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_y[i] = -2.0 * tr_zz_y[i] * tbe_0 + 4.0 * tr_xzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_z[i] = -2.0 * tr_zz_z[i] * tbe_0 + 4.0 * tr_xzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 18-21 components of targeted buffer : DP

    auto tr_x_0_y_xx_x = pbuffer.data(idx_op_geom_110_dp + 18);

    auto tr_x_0_y_xx_y = pbuffer.data(idx_op_geom_110_dp + 19);

    auto tr_x_0_y_xx_z = pbuffer.data(idx_op_geom_110_dp + 20);

    #pragma omp simd aligned(tr_x_0, tr_x_0_y_xx_x, tr_x_0_y_xx_y, tr_x_0_y_xx_z, tr_x_xy, tr_x_yy, tr_x_yz, tr_xxx_0, tr_xxx_xy, tr_xxx_yy, tr_xxx_yz, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xy_x, tr_xy_y, tr_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xx_x[i] = -4.0 * tr_x_xy[i] * tke_0 - 4.0 * tr_xy_x[i] * tbe_0 + 4.0 * tr_xxx_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_y[i] = 2.0 * tr_x_0[i] - 4.0 * tr_x_yy[i] * tke_0 - 4.0 * tr_xy_y[i] * tbe_0 - 2.0 * tr_xxx_0[i] * tbe_0 + 4.0 * tr_xxx_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_z[i] = -4.0 * tr_x_yz[i] * tke_0 - 4.0 * tr_xy_z[i] * tbe_0 + 4.0 * tr_xxx_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 21-24 components of targeted buffer : DP

    auto tr_x_0_y_xy_x = pbuffer.data(idx_op_geom_110_dp + 21);

    auto tr_x_0_y_xy_y = pbuffer.data(idx_op_geom_110_dp + 22);

    auto tr_x_0_y_xy_z = pbuffer.data(idx_op_geom_110_dp + 23);

    #pragma omp simd aligned(tr_0_x, tr_0_y, tr_0_z, tr_x_0_y_xy_x, tr_x_0_y_xy_y, tr_x_0_y_xy_z, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxy_0, tr_xxy_xy, tr_xxy_yy, tr_xxy_yz, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_y_0, tr_y_xy, tr_y_yy, tr_y_yz, tr_yy_x, tr_yy_y, tr_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xy_x[i] = tr_0_x[i] - 2.0 * tr_y_xy[i] * tke_0 - 2.0 * tr_yy_x[i] * tbe_0 - 2.0 * tr_xx_x[i] * tbe_0 + 4.0 * tr_xxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_y[i] = tr_0_y[i] + tr_y_0[i] - 2.0 * tr_y_yy[i] * tke_0 - 2.0 * tr_yy_y[i] * tbe_0 - 2.0 * tr_xx_y[i] * tbe_0 - 2.0 * tr_xxy_0[i] * tbe_0 + 4.0 * tr_xxy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_z[i] = tr_0_z[i] - 2.0 * tr_y_yz[i] * tke_0 - 2.0 * tr_yy_z[i] * tbe_0 - 2.0 * tr_xx_z[i] * tbe_0 + 4.0 * tr_xxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 24-27 components of targeted buffer : DP

    auto tr_x_0_y_xz_x = pbuffer.data(idx_op_geom_110_dp + 24);

    auto tr_x_0_y_xz_y = pbuffer.data(idx_op_geom_110_dp + 25);

    auto tr_x_0_y_xz_z = pbuffer.data(idx_op_geom_110_dp + 26);

    #pragma omp simd aligned(tr_x_0_y_xz_x, tr_x_0_y_xz_y, tr_x_0_y_xz_z, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxz_0, tr_xxz_xy, tr_xxz_yy, tr_xxz_yz, tr_yz_x, tr_yz_y, tr_yz_z, tr_z_0, tr_z_xy, tr_z_yy, tr_z_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xz_x[i] = -2.0 * tr_z_xy[i] * tke_0 - 2.0 * tr_yz_x[i] * tbe_0 + 4.0 * tr_xxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_y[i] = tr_z_0[i] - 2.0 * tr_z_yy[i] * tke_0 - 2.0 * tr_yz_y[i] * tbe_0 - 2.0 * tr_xxz_0[i] * tbe_0 + 4.0 * tr_xxz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_z[i] = -2.0 * tr_z_yz[i] * tke_0 - 2.0 * tr_yz_z[i] * tbe_0 + 4.0 * tr_xxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 27-30 components of targeted buffer : DP

    auto tr_x_0_y_yy_x = pbuffer.data(idx_op_geom_110_dp + 27);

    auto tr_x_0_y_yy_y = pbuffer.data(idx_op_geom_110_dp + 28);

    auto tr_x_0_y_yy_z = pbuffer.data(idx_op_geom_110_dp + 29);

    #pragma omp simd aligned(tr_x_0_y_yy_x, tr_x_0_y_yy_y, tr_x_0_y_yy_z, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyy_0, tr_xyy_xy, tr_xyy_yy, tr_xyy_yz, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yy_x[i] = -4.0 * tr_xy_x[i] * tbe_0 + 4.0 * tr_xyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_y[i] = -4.0 * tr_xy_y[i] * tbe_0 - 2.0 * tr_xyy_0[i] * tbe_0 + 4.0 * tr_xyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_z[i] = -4.0 * tr_xy_z[i] * tbe_0 + 4.0 * tr_xyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 30-33 components of targeted buffer : DP

    auto tr_x_0_y_yz_x = pbuffer.data(idx_op_geom_110_dp + 30);

    auto tr_x_0_y_yz_y = pbuffer.data(idx_op_geom_110_dp + 31);

    auto tr_x_0_y_yz_z = pbuffer.data(idx_op_geom_110_dp + 32);

    #pragma omp simd aligned(tr_x_0_y_yz_x, tr_x_0_y_yz_y, tr_x_0_y_yz_z, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyz_0, tr_xyz_xy, tr_xyz_yy, tr_xyz_yz, tr_xz_x, tr_xz_y, tr_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yz_x[i] = -2.0 * tr_xz_x[i] * tbe_0 + 4.0 * tr_xyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_y[i] = -2.0 * tr_xz_y[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_z[i] = -2.0 * tr_xz_z[i] * tbe_0 + 4.0 * tr_xyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 33-36 components of targeted buffer : DP

    auto tr_x_0_y_zz_x = pbuffer.data(idx_op_geom_110_dp + 33);

    auto tr_x_0_y_zz_y = pbuffer.data(idx_op_geom_110_dp + 34);

    auto tr_x_0_y_zz_z = pbuffer.data(idx_op_geom_110_dp + 35);

    #pragma omp simd aligned(tr_x_0_y_zz_x, tr_x_0_y_zz_y, tr_x_0_y_zz_z, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xzz_0, tr_xzz_xy, tr_xzz_yy, tr_xzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_zz_x[i] = 4.0 * tr_xzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_y[i] = -2.0 * tr_xzz_0[i] * tbe_0 + 4.0 * tr_xzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_z[i] = 4.0 * tr_xzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 36-39 components of targeted buffer : DP

    auto tr_x_0_z_xx_x = pbuffer.data(idx_op_geom_110_dp + 36);

    auto tr_x_0_z_xx_y = pbuffer.data(idx_op_geom_110_dp + 37);

    auto tr_x_0_z_xx_z = pbuffer.data(idx_op_geom_110_dp + 38);

    #pragma omp simd aligned(tr_x_0, tr_x_0_z_xx_x, tr_x_0_z_xx_y, tr_x_0_z_xx_z, tr_x_xz, tr_x_yz, tr_x_zz, tr_xxx_0, tr_xxx_xz, tr_xxx_yz, tr_xxx_zz, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xz_x, tr_xz_y, tr_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xx_x[i] = -4.0 * tr_x_xz[i] * tke_0 - 4.0 * tr_xz_x[i] * tbe_0 + 4.0 * tr_xxx_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_y[i] = -4.0 * tr_x_yz[i] * tke_0 - 4.0 * tr_xz_y[i] * tbe_0 + 4.0 * tr_xxx_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_z[i] = 2.0 * tr_x_0[i] - 4.0 * tr_x_zz[i] * tke_0 - 4.0 * tr_xz_z[i] * tbe_0 - 2.0 * tr_xxx_0[i] * tbe_0 + 4.0 * tr_xxx_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 39-42 components of targeted buffer : DP

    auto tr_x_0_z_xy_x = pbuffer.data(idx_op_geom_110_dp + 39);

    auto tr_x_0_z_xy_y = pbuffer.data(idx_op_geom_110_dp + 40);

    auto tr_x_0_z_xy_z = pbuffer.data(idx_op_geom_110_dp + 41);

    #pragma omp simd aligned(tr_x_0_z_xy_x, tr_x_0_z_xy_y, tr_x_0_z_xy_z, tr_xxy_0, tr_xxy_xz, tr_xxy_yz, tr_xxy_zz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_y_0, tr_y_xz, tr_y_yz, tr_y_zz, tr_yz_x, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xy_x[i] = -2.0 * tr_y_xz[i] * tke_0 - 2.0 * tr_yz_x[i] * tbe_0 + 4.0 * tr_xxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_y[i] = -2.0 * tr_y_yz[i] * tke_0 - 2.0 * tr_yz_y[i] * tbe_0 + 4.0 * tr_xxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_z[i] = tr_y_0[i] - 2.0 * tr_y_zz[i] * tke_0 - 2.0 * tr_yz_z[i] * tbe_0 - 2.0 * tr_xxy_0[i] * tbe_0 + 4.0 * tr_xxy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 42-45 components of targeted buffer : DP

    auto tr_x_0_z_xz_x = pbuffer.data(idx_op_geom_110_dp + 42);

    auto tr_x_0_z_xz_y = pbuffer.data(idx_op_geom_110_dp + 43);

    auto tr_x_0_z_xz_z = pbuffer.data(idx_op_geom_110_dp + 44);

    #pragma omp simd aligned(tr_0_x, tr_0_y, tr_0_z, tr_x_0_z_xz_x, tr_x_0_z_xz_y, tr_x_0_z_xz_z, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxz_0, tr_xxz_xz, tr_xxz_yz, tr_xxz_zz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_z_0, tr_z_xz, tr_z_yz, tr_z_zz, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xz_x[i] = tr_0_x[i] - 2.0 * tr_z_xz[i] * tke_0 - 2.0 * tr_zz_x[i] * tbe_0 - 2.0 * tr_xx_x[i] * tbe_0 + 4.0 * tr_xxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_y[i] = tr_0_y[i] - 2.0 * tr_z_yz[i] * tke_0 - 2.0 * tr_zz_y[i] * tbe_0 - 2.0 * tr_xx_y[i] * tbe_0 + 4.0 * tr_xxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_z[i] = tr_0_z[i] + tr_z_0[i] - 2.0 * tr_z_zz[i] * tke_0 - 2.0 * tr_zz_z[i] * tbe_0 - 2.0 * tr_xx_z[i] * tbe_0 - 2.0 * tr_xxz_0[i] * tbe_0 + 4.0 * tr_xxz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 45-48 components of targeted buffer : DP

    auto tr_x_0_z_yy_x = pbuffer.data(idx_op_geom_110_dp + 45);

    auto tr_x_0_z_yy_y = pbuffer.data(idx_op_geom_110_dp + 46);

    auto tr_x_0_z_yy_z = pbuffer.data(idx_op_geom_110_dp + 47);

    #pragma omp simd aligned(tr_x_0_z_yy_x, tr_x_0_z_yy_y, tr_x_0_z_yy_z, tr_xyy_0, tr_xyy_xz, tr_xyy_yz, tr_xyy_zz, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yy_x[i] = 4.0 * tr_xyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_y[i] = 4.0 * tr_xyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_z[i] = -2.0 * tr_xyy_0[i] * tbe_0 + 4.0 * tr_xyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 48-51 components of targeted buffer : DP

    auto tr_x_0_z_yz_x = pbuffer.data(idx_op_geom_110_dp + 48);

    auto tr_x_0_z_yz_y = pbuffer.data(idx_op_geom_110_dp + 49);

    auto tr_x_0_z_yz_z = pbuffer.data(idx_op_geom_110_dp + 50);

    #pragma omp simd aligned(tr_x_0_z_yz_x, tr_x_0_z_yz_y, tr_x_0_z_yz_z, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyz_0, tr_xyz_xz, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yz_x[i] = -2.0 * tr_xy_x[i] * tbe_0 + 4.0 * tr_xyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_y[i] = -2.0 * tr_xy_y[i] * tbe_0 + 4.0 * tr_xyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_z[i] = -2.0 * tr_xy_z[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 51-54 components of targeted buffer : DP

    auto tr_x_0_z_zz_x = pbuffer.data(idx_op_geom_110_dp + 51);

    auto tr_x_0_z_zz_y = pbuffer.data(idx_op_geom_110_dp + 52);

    auto tr_x_0_z_zz_z = pbuffer.data(idx_op_geom_110_dp + 53);

    #pragma omp simd aligned(tr_x_0_z_zz_x, tr_x_0_z_zz_y, tr_x_0_z_zz_z, tr_xz_x, tr_xz_y, tr_xz_z, tr_xzz_0, tr_xzz_xz, tr_xzz_yz, tr_xzz_zz, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_zz_x[i] = -4.0 * tr_xz_x[i] * tbe_0 + 4.0 * tr_xzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_y[i] = -4.0 * tr_xz_y[i] * tbe_0 + 4.0 * tr_xzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_z[i] = -4.0 * tr_xz_z[i] * tbe_0 - 2.0 * tr_xzz_0[i] * tbe_0 + 4.0 * tr_xzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 54-57 components of targeted buffer : DP

    auto tr_y_0_x_xx_x = pbuffer.data(idx_op_geom_110_dp + 54);

    auto tr_y_0_x_xx_y = pbuffer.data(idx_op_geom_110_dp + 55);

    auto tr_y_0_x_xx_z = pbuffer.data(idx_op_geom_110_dp + 56);

    #pragma omp simd aligned(tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxy_0, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xy_x, tr_xy_y, tr_xy_z, tr_y_0_x_xx_x, tr_y_0_x_xx_y, tr_y_0_x_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xx_x[i] = -4.0 * tr_xy_x[i] * tbe_0 - 2.0 * tr_xxy_0[i] * tbe_0 + 4.0 * tr_xxy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_y[i] = -4.0 * tr_xy_y[i] * tbe_0 + 4.0 * tr_xxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_z[i] = -4.0 * tr_xy_z[i] * tbe_0 + 4.0 * tr_xxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 57-60 components of targeted buffer : DP

    auto tr_y_0_x_xy_x = pbuffer.data(idx_op_geom_110_dp + 57);

    auto tr_y_0_x_xy_y = pbuffer.data(idx_op_geom_110_dp + 58);

    auto tr_y_0_x_xy_z = pbuffer.data(idx_op_geom_110_dp + 59);

    #pragma omp simd aligned(tr_0_x, tr_0_y, tr_0_z, tr_x_0, tr_x_xx, tr_x_xy, tr_x_xz, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xyy_0, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_y_0_x_xy_x, tr_y_0_x_xy_y, tr_y_0_x_xy_z, tr_yy_x, tr_yy_y, tr_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xy_x[i] = tr_0_x[i] - 2.0 * tr_yy_x[i] * tbe_0 + tr_x_0[i] - 2.0 * tr_x_xx[i] * tke_0 - 2.0 * tr_xyy_0[i] * tbe_0 + 4.0 * tr_xyy_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xx_x[i] * tbe_0 + 4.0 * tr_xxyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_y[i] = tr_0_y[i] - 2.0 * tr_yy_y[i] * tbe_0 - 2.0 * tr_x_xy[i] * tke_0 + 4.0 * tr_xyy_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_y[i] * tbe_0 + 4.0 * tr_xxyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_z[i] = tr_0_z[i] - 2.0 * tr_yy_z[i] * tbe_0 - 2.0 * tr_x_xz[i] * tke_0 + 4.0 * tr_xyy_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_z[i] * tbe_0 + 4.0 * tr_xxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 60-63 components of targeted buffer : DP

    auto tr_y_0_x_xz_x = pbuffer.data(idx_op_geom_110_dp + 60);

    auto tr_y_0_x_xz_y = pbuffer.data(idx_op_geom_110_dp + 61);

    auto tr_y_0_x_xz_z = pbuffer.data(idx_op_geom_110_dp + 62);

    #pragma omp simd aligned(tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xyz_0, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_y_0_x_xz_x, tr_y_0_x_xz_y, tr_y_0_x_xz_z, tr_yz_x, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xz_x[i] = -2.0 * tr_yz_x[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_y[i] = -2.0 * tr_yz_y[i] * tbe_0 + 4.0 * tr_xyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_z[i] = -2.0 * tr_yz_z[i] * tbe_0 + 4.0 * tr_xyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 63-66 components of targeted buffer : DP

    auto tr_y_0_x_yy_x = pbuffer.data(idx_op_geom_110_dp + 63);

    auto tr_y_0_x_yy_y = pbuffer.data(idx_op_geom_110_dp + 64);

    auto tr_y_0_x_yy_z = pbuffer.data(idx_op_geom_110_dp + 65);

    #pragma omp simd aligned(tr_xy_x, tr_xy_y, tr_xy_z, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_y_0, tr_y_0_x_yy_x, tr_y_0_x_yy_y, tr_y_0_x_yy_z, tr_y_xx, tr_y_xy, tr_y_xz, tr_yyy_0, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yy_x[i] = 2.0 * tr_y_0[i] - 4.0 * tr_y_xx[i] * tke_0 - 2.0 * tr_yyy_0[i] * tbe_0 + 4.0 * tr_yyy_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xy_x[i] * tbe_0 + 4.0 * tr_xyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_y[i] = -4.0 * tr_y_xy[i] * tke_0 + 4.0 * tr_yyy_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xy_y[i] * tbe_0 + 4.0 * tr_xyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_z[i] = -4.0 * tr_y_xz[i] * tke_0 + 4.0 * tr_yyy_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_z[i] * tbe_0 + 4.0 * tr_xyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 66-69 components of targeted buffer : DP

    auto tr_y_0_x_yz_x = pbuffer.data(idx_op_geom_110_dp + 66);

    auto tr_y_0_x_yz_y = pbuffer.data(idx_op_geom_110_dp + 67);

    auto tr_y_0_x_yz_z = pbuffer.data(idx_op_geom_110_dp + 68);

    #pragma omp simd aligned(tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xz_x, tr_xz_y, tr_xz_z, tr_y_0_x_yz_x, tr_y_0_x_yz_y, tr_y_0_x_yz_z, tr_yyz_0, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_z_0, tr_z_xx, tr_z_xy, tr_z_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yz_x[i] = tr_z_0[i] - 2.0 * tr_z_xx[i] * tke_0 - 2.0 * tr_yyz_0[i] * tbe_0 + 4.0 * tr_yyz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xz_x[i] * tbe_0 + 4.0 * tr_xyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_y[i] = -2.0 * tr_z_xy[i] * tke_0 + 4.0 * tr_yyz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_y[i] * tbe_0 + 4.0 * tr_xyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_z[i] = -2.0 * tr_z_xz[i] * tke_0 + 4.0 * tr_yyz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_z[i] * tbe_0 + 4.0 * tr_xyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 69-72 components of targeted buffer : DP

    auto tr_y_0_x_zz_x = pbuffer.data(idx_op_geom_110_dp + 69);

    auto tr_y_0_x_zz_y = pbuffer.data(idx_op_geom_110_dp + 70);

    auto tr_y_0_x_zz_z = pbuffer.data(idx_op_geom_110_dp + 71);

    #pragma omp simd aligned(tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_y_0_x_zz_x, tr_y_0_x_zz_y, tr_y_0_x_zz_z, tr_yzz_0, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_zz_x[i] = -2.0 * tr_yzz_0[i] * tbe_0 + 4.0 * tr_yzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_y[i] = 4.0 * tr_yzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_z[i] = 4.0 * tr_yzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 72-75 components of targeted buffer : DP

    auto tr_y_0_y_xx_x = pbuffer.data(idx_op_geom_110_dp + 72);

    auto tr_y_0_y_xx_y = pbuffer.data(idx_op_geom_110_dp + 73);

    auto tr_y_0_y_xx_z = pbuffer.data(idx_op_geom_110_dp + 74);

    #pragma omp simd aligned(tr_xx_x, tr_xx_y, tr_xx_z, tr_xxy_0, tr_xxy_xy, tr_xxy_yy, tr_xxy_yz, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_y_0_y_xx_x, tr_y_0_y_xx_y, tr_y_0_y_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xx_x[i] = -2.0 * tr_xx_x[i] * tbe_0 + 4.0 * tr_xxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_y[i] = -2.0 * tr_xx_y[i] * tbe_0 - 2.0 * tr_xxy_0[i] * tbe_0 + 4.0 * tr_xxy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_z[i] = -2.0 * tr_xx_z[i] * tbe_0 + 4.0 * tr_xxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 75-78 components of targeted buffer : DP

    auto tr_y_0_y_xy_x = pbuffer.data(idx_op_geom_110_dp + 75);

    auto tr_y_0_y_xy_y = pbuffer.data(idx_op_geom_110_dp + 76);

    auto tr_y_0_y_xy_z = pbuffer.data(idx_op_geom_110_dp + 77);

    #pragma omp simd aligned(tr_x_0, tr_x_xy, tr_x_yy, tr_x_yz, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyy_0, tr_xyy_xy, tr_xyy_yy, tr_xyy_yz, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_y_0_y_xy_x, tr_y_0_y_xy_y, tr_y_0_y_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xy_x[i] = -2.0 * tr_x_xy[i] * tke_0 - 6.0 * tr_xy_x[i] * tbe_0 + 4.0 * tr_xyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_y[i] = tr_x_0[i] - 2.0 * tr_x_yy[i] * tke_0 - 6.0 * tr_xy_y[i] * tbe_0 - 2.0 * tr_xyy_0[i] * tbe_0 + 4.0 * tr_xyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_z[i] = -2.0 * tr_x_yz[i] * tke_0 - 6.0 * tr_xy_z[i] * tbe_0 + 4.0 * tr_xyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 78-81 components of targeted buffer : DP

    auto tr_y_0_y_xz_x = pbuffer.data(idx_op_geom_110_dp + 78);

    auto tr_y_0_y_xz_y = pbuffer.data(idx_op_geom_110_dp + 79);

    auto tr_y_0_y_xz_z = pbuffer.data(idx_op_geom_110_dp + 80);

    #pragma omp simd aligned(tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyz_0, tr_xyz_xy, tr_xyz_yy, tr_xyz_yz, tr_xz_x, tr_xz_y, tr_xz_z, tr_y_0_y_xz_x, tr_y_0_y_xz_y, tr_y_0_y_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xz_x[i] = -2.0 * tr_xz_x[i] * tbe_0 + 4.0 * tr_xyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_y[i] = -2.0 * tr_xz_y[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_z[i] = -2.0 * tr_xz_z[i] * tbe_0 + 4.0 * tr_xyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 81-84 components of targeted buffer : DP

    auto tr_y_0_y_yy_x = pbuffer.data(idx_op_geom_110_dp + 81);

    auto tr_y_0_y_yy_y = pbuffer.data(idx_op_geom_110_dp + 82);

    auto tr_y_0_y_yy_z = pbuffer.data(idx_op_geom_110_dp + 83);

    #pragma omp simd aligned(tr_0_x, tr_0_y, tr_0_z, tr_y_0, tr_y_0_y_yy_x, tr_y_0_y_yy_y, tr_y_0_y_yy_z, tr_y_xy, tr_y_yy, tr_y_yz, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyy_0, tr_yyy_xy, tr_yyy_yy, tr_yyy_yz, tr_yyyy_x, tr_yyyy_y, tr_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yy_x[i] = 2.0 * tr_0_x[i] - 4.0 * tr_y_xy[i] * tke_0 - 10.0 * tr_yy_x[i] * tbe_0 + 4.0 * tr_yyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_y[i] = 2.0 * tr_0_y[i] + 2.0 * tr_y_0[i] - 4.0 * tr_y_yy[i] * tke_0 - 10.0 * tr_yy_y[i] * tbe_0 - 2.0 * tr_yyy_0[i] * tbe_0 + 4.0 * tr_yyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_z[i] = 2.0 * tr_0_z[i] - 4.0 * tr_y_yz[i] * tke_0 - 10.0 * tr_yy_z[i] * tbe_0 + 4.0 * tr_yyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 84-87 components of targeted buffer : DP

    auto tr_y_0_y_yz_x = pbuffer.data(idx_op_geom_110_dp + 84);

    auto tr_y_0_y_yz_y = pbuffer.data(idx_op_geom_110_dp + 85);

    auto tr_y_0_y_yz_z = pbuffer.data(idx_op_geom_110_dp + 86);

    #pragma omp simd aligned(tr_y_0_y_yz_x, tr_y_0_y_yz_y, tr_y_0_y_yz_z, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yyz_0, tr_yyz_xy, tr_yyz_yy, tr_yyz_yz, tr_yz_x, tr_yz_y, tr_yz_z, tr_z_0, tr_z_xy, tr_z_yy, tr_z_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yz_x[i] = -2.0 * tr_z_xy[i] * tke_0 - 6.0 * tr_yz_x[i] * tbe_0 + 4.0 * tr_yyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_y[i] = tr_z_0[i] - 2.0 * tr_z_yy[i] * tke_0 - 6.0 * tr_yz_y[i] * tbe_0 - 2.0 * tr_yyz_0[i] * tbe_0 + 4.0 * tr_yyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_z[i] = -2.0 * tr_z_yz[i] * tke_0 - 6.0 * tr_yz_z[i] * tbe_0 + 4.0 * tr_yyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 87-90 components of targeted buffer : DP

    auto tr_y_0_y_zz_x = pbuffer.data(idx_op_geom_110_dp + 87);

    auto tr_y_0_y_zz_y = pbuffer.data(idx_op_geom_110_dp + 88);

    auto tr_y_0_y_zz_z = pbuffer.data(idx_op_geom_110_dp + 89);

    #pragma omp simd aligned(tr_y_0_y_zz_x, tr_y_0_y_zz_y, tr_y_0_y_zz_z, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_yzz_0, tr_yzz_xy, tr_yzz_yy, tr_yzz_yz, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_zz_x[i] = -2.0 * tr_zz_x[i] * tbe_0 + 4.0 * tr_yzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_y[i] = -2.0 * tr_zz_y[i] * tbe_0 - 2.0 * tr_yzz_0[i] * tbe_0 + 4.0 * tr_yzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_z[i] = -2.0 * tr_zz_z[i] * tbe_0 + 4.0 * tr_yzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 90-93 components of targeted buffer : DP

    auto tr_y_0_z_xx_x = pbuffer.data(idx_op_geom_110_dp + 90);

    auto tr_y_0_z_xx_y = pbuffer.data(idx_op_geom_110_dp + 91);

    auto tr_y_0_z_xx_z = pbuffer.data(idx_op_geom_110_dp + 92);

    #pragma omp simd aligned(tr_xxy_0, tr_xxy_xz, tr_xxy_yz, tr_xxy_zz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_y_0_z_xx_x, tr_y_0_z_xx_y, tr_y_0_z_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xx_x[i] = 4.0 * tr_xxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_y[i] = 4.0 * tr_xxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_z[i] = -2.0 * tr_xxy_0[i] * tbe_0 + 4.0 * tr_xxy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 93-96 components of targeted buffer : DP

    auto tr_y_0_z_xy_x = pbuffer.data(idx_op_geom_110_dp + 93);

    auto tr_y_0_z_xy_y = pbuffer.data(idx_op_geom_110_dp + 94);

    auto tr_y_0_z_xy_z = pbuffer.data(idx_op_geom_110_dp + 95);

    #pragma omp simd aligned(tr_x_0, tr_x_xz, tr_x_yz, tr_x_zz, tr_xyy_0, tr_xyy_xz, tr_xyy_yz, tr_xyy_zz, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xz_x, tr_xz_y, tr_xz_z, tr_y_0_z_xy_x, tr_y_0_z_xy_y, tr_y_0_z_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xy_x[i] = -2.0 * tr_x_xz[i] * tke_0 - 2.0 * tr_xz_x[i] * tbe_0 + 4.0 * tr_xyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_y[i] = -2.0 * tr_x_yz[i] * tke_0 - 2.0 * tr_xz_y[i] * tbe_0 + 4.0 * tr_xyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_z[i] = tr_x_0[i] - 2.0 * tr_x_zz[i] * tke_0 - 2.0 * tr_xz_z[i] * tbe_0 - 2.0 * tr_xyy_0[i] * tbe_0 + 4.0 * tr_xyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 96-99 components of targeted buffer : DP

    auto tr_y_0_z_xz_x = pbuffer.data(idx_op_geom_110_dp + 96);

    auto tr_y_0_z_xz_y = pbuffer.data(idx_op_geom_110_dp + 97);

    auto tr_y_0_z_xz_z = pbuffer.data(idx_op_geom_110_dp + 98);

    #pragma omp simd aligned(tr_xy_x, tr_xy_y, tr_xy_z, tr_xyz_0, tr_xyz_xz, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_y_0_z_xz_x, tr_y_0_z_xz_y, tr_y_0_z_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xz_x[i] = -2.0 * tr_xy_x[i] * tbe_0 + 4.0 * tr_xyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_y[i] = -2.0 * tr_xy_y[i] * tbe_0 + 4.0 * tr_xyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_z[i] = -2.0 * tr_xy_z[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 99-102 components of targeted buffer : DP

    auto tr_y_0_z_yy_x = pbuffer.data(idx_op_geom_110_dp + 99);

    auto tr_y_0_z_yy_y = pbuffer.data(idx_op_geom_110_dp + 100);

    auto tr_y_0_z_yy_z = pbuffer.data(idx_op_geom_110_dp + 101);

    #pragma omp simd aligned(tr_y_0, tr_y_0_z_yy_x, tr_y_0_z_yy_y, tr_y_0_z_yy_z, tr_y_xz, tr_y_yz, tr_y_zz, tr_yyy_0, tr_yyy_xz, tr_yyy_yz, tr_yyy_zz, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yz_x, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yy_x[i] = -4.0 * tr_y_xz[i] * tke_0 - 4.0 * tr_yz_x[i] * tbe_0 + 4.0 * tr_yyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_y[i] = -4.0 * tr_y_yz[i] * tke_0 - 4.0 * tr_yz_y[i] * tbe_0 + 4.0 * tr_yyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_z[i] = 2.0 * tr_y_0[i] - 4.0 * tr_y_zz[i] * tke_0 - 4.0 * tr_yz_z[i] * tbe_0 - 2.0 * tr_yyy_0[i] * tbe_0 + 4.0 * tr_yyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 102-105 components of targeted buffer : DP

    auto tr_y_0_z_yz_x = pbuffer.data(idx_op_geom_110_dp + 102);

    auto tr_y_0_z_yz_y = pbuffer.data(idx_op_geom_110_dp + 103);

    auto tr_y_0_z_yz_z = pbuffer.data(idx_op_geom_110_dp + 104);

    #pragma omp simd aligned(tr_0_x, tr_0_y, tr_0_z, tr_y_0_z_yz_x, tr_y_0_z_yz_y, tr_y_0_z_yz_z, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyz_0, tr_yyz_xz, tr_yyz_yz, tr_yyz_zz, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_z_0, tr_z_xz, tr_z_yz, tr_z_zz, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yz_x[i] = tr_0_x[i] - 2.0 * tr_z_xz[i] * tke_0 - 2.0 * tr_zz_x[i] * tbe_0 - 2.0 * tr_yy_x[i] * tbe_0 + 4.0 * tr_yyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_y[i] = tr_0_y[i] - 2.0 * tr_z_yz[i] * tke_0 - 2.0 * tr_zz_y[i] * tbe_0 - 2.0 * tr_yy_y[i] * tbe_0 + 4.0 * tr_yyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_z[i] = tr_0_z[i] + tr_z_0[i] - 2.0 * tr_z_zz[i] * tke_0 - 2.0 * tr_zz_z[i] * tbe_0 - 2.0 * tr_yy_z[i] * tbe_0 - 2.0 * tr_yyz_0[i] * tbe_0 + 4.0 * tr_yyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 105-108 components of targeted buffer : DP

    auto tr_y_0_z_zz_x = pbuffer.data(idx_op_geom_110_dp + 105);

    auto tr_y_0_z_zz_y = pbuffer.data(idx_op_geom_110_dp + 106);

    auto tr_y_0_z_zz_z = pbuffer.data(idx_op_geom_110_dp + 107);

    #pragma omp simd aligned(tr_y_0_z_zz_x, tr_y_0_z_zz_y, tr_y_0_z_zz_z, tr_yz_x, tr_yz_y, tr_yz_z, tr_yzz_0, tr_yzz_xz, tr_yzz_yz, tr_yzz_zz, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_zz_x[i] = -4.0 * tr_yz_x[i] * tbe_0 + 4.0 * tr_yzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_y[i] = -4.0 * tr_yz_y[i] * tbe_0 + 4.0 * tr_yzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_z[i] = -4.0 * tr_yz_z[i] * tbe_0 - 2.0 * tr_yzz_0[i] * tbe_0 + 4.0 * tr_yzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 108-111 components of targeted buffer : DP

    auto tr_z_0_x_xx_x = pbuffer.data(idx_op_geom_110_dp + 108);

    auto tr_z_0_x_xx_y = pbuffer.data(idx_op_geom_110_dp + 109);

    auto tr_z_0_x_xx_z = pbuffer.data(idx_op_geom_110_dp + 110);

    #pragma omp simd aligned(tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xxz_0, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xz_x, tr_xz_y, tr_xz_z, tr_z_0_x_xx_x, tr_z_0_x_xx_y, tr_z_0_x_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xx_x[i] = -4.0 * tr_xz_x[i] * tbe_0 - 2.0 * tr_xxz_0[i] * tbe_0 + 4.0 * tr_xxz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_y[i] = -4.0 * tr_xz_y[i] * tbe_0 + 4.0 * tr_xxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_z[i] = -4.0 * tr_xz_z[i] * tbe_0 + 4.0 * tr_xxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 111-114 components of targeted buffer : DP

    auto tr_z_0_x_xy_x = pbuffer.data(idx_op_geom_110_dp + 111);

    auto tr_z_0_x_xy_y = pbuffer.data(idx_op_geom_110_dp + 112);

    auto tr_z_0_x_xy_z = pbuffer.data(idx_op_geom_110_dp + 113);

    #pragma omp simd aligned(tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xyz_0, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_yz_x, tr_yz_y, tr_yz_z, tr_z_0_x_xy_x, tr_z_0_x_xy_y, tr_z_0_x_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xy_x[i] = -2.0 * tr_yz_x[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_y[i] = -2.0 * tr_yz_y[i] * tbe_0 + 4.0 * tr_xyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_z[i] = -2.0 * tr_yz_z[i] * tbe_0 + 4.0 * tr_xyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 114-117 components of targeted buffer : DP

    auto tr_z_0_x_xz_x = pbuffer.data(idx_op_geom_110_dp + 114);

    auto tr_z_0_x_xz_y = pbuffer.data(idx_op_geom_110_dp + 115);

    auto tr_z_0_x_xz_z = pbuffer.data(idx_op_geom_110_dp + 116);

    #pragma omp simd aligned(tr_0_x, tr_0_y, tr_0_z, tr_x_0, tr_x_xx, tr_x_xy, tr_x_xz, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_xzz_0, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_z_0_x_xz_x, tr_z_0_x_xz_y, tr_z_0_x_xz_z, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xz_x[i] = tr_0_x[i] - 2.0 * tr_zz_x[i] * tbe_0 + tr_x_0[i] - 2.0 * tr_x_xx[i] * tke_0 - 2.0 * tr_xzz_0[i] * tbe_0 + 4.0 * tr_xzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xx_x[i] * tbe_0 + 4.0 * tr_xxzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_y[i] = tr_0_y[i] - 2.0 * tr_zz_y[i] * tbe_0 - 2.0 * tr_x_xy[i] * tke_0 + 4.0 * tr_xzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_y[i] * tbe_0 + 4.0 * tr_xxzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_z[i] = tr_0_z[i] - 2.0 * tr_zz_z[i] * tbe_0 - 2.0 * tr_x_xz[i] * tke_0 + 4.0 * tr_xzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_z[i] * tbe_0 + 4.0 * tr_xxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 117-120 components of targeted buffer : DP

    auto tr_z_0_x_yy_x = pbuffer.data(idx_op_geom_110_dp + 117);

    auto tr_z_0_x_yy_y = pbuffer.data(idx_op_geom_110_dp + 118);

    auto tr_z_0_x_yy_z = pbuffer.data(idx_op_geom_110_dp + 119);

    #pragma omp simd aligned(tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_yyz_0, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_z_0_x_yy_x, tr_z_0_x_yy_y, tr_z_0_x_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yy_x[i] = -2.0 * tr_yyz_0[i] * tbe_0 + 4.0 * tr_yyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_y[i] = 4.0 * tr_yyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_z[i] = 4.0 * tr_yyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 120-123 components of targeted buffer : DP

    auto tr_z_0_x_yz_x = pbuffer.data(idx_op_geom_110_dp + 120);

    auto tr_z_0_x_yz_y = pbuffer.data(idx_op_geom_110_dp + 121);

    auto tr_z_0_x_yz_z = pbuffer.data(idx_op_geom_110_dp + 122);

    #pragma omp simd aligned(tr_xy_x, tr_xy_y, tr_xy_z, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_y_0, tr_y_xx, tr_y_xy, tr_y_xz, tr_yzz_0, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_z_0_x_yz_x, tr_z_0_x_yz_y, tr_z_0_x_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yz_x[i] = tr_y_0[i] - 2.0 * tr_y_xx[i] * tke_0 - 2.0 * tr_yzz_0[i] * tbe_0 + 4.0 * tr_yzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xy_x[i] * tbe_0 + 4.0 * tr_xyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_y[i] = -2.0 * tr_y_xy[i] * tke_0 + 4.0 * tr_yzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_y[i] * tbe_0 + 4.0 * tr_xyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_z[i] = -2.0 * tr_y_xz[i] * tke_0 + 4.0 * tr_yzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_z[i] * tbe_0 + 4.0 * tr_xyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 123-126 components of targeted buffer : DP

    auto tr_z_0_x_zz_x = pbuffer.data(idx_op_geom_110_dp + 123);

    auto tr_z_0_x_zz_y = pbuffer.data(idx_op_geom_110_dp + 124);

    auto tr_z_0_x_zz_z = pbuffer.data(idx_op_geom_110_dp + 125);

    #pragma omp simd aligned(tr_xz_x, tr_xz_y, tr_xz_z, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_z_0, tr_z_0_x_zz_x, tr_z_0_x_zz_y, tr_z_0_x_zz_z, tr_z_xx, tr_z_xy, tr_z_xz, tr_zzz_0, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_zz_x[i] = 2.0 * tr_z_0[i] - 4.0 * tr_z_xx[i] * tke_0 - 2.0 * tr_zzz_0[i] * tbe_0 + 4.0 * tr_zzz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xz_x[i] * tbe_0 + 4.0 * tr_xzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_y[i] = -4.0 * tr_z_xy[i] * tke_0 + 4.0 * tr_zzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_y[i] * tbe_0 + 4.0 * tr_xzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_z[i] = -4.0 * tr_z_xz[i] * tke_0 + 4.0 * tr_zzz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_z[i] * tbe_0 + 4.0 * tr_xzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 126-129 components of targeted buffer : DP

    auto tr_z_0_y_xx_x = pbuffer.data(idx_op_geom_110_dp + 126);

    auto tr_z_0_y_xx_y = pbuffer.data(idx_op_geom_110_dp + 127);

    auto tr_z_0_y_xx_z = pbuffer.data(idx_op_geom_110_dp + 128);

    #pragma omp simd aligned(tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxz_0, tr_xxz_xy, tr_xxz_yy, tr_xxz_yz, tr_z_0_y_xx_x, tr_z_0_y_xx_y, tr_z_0_y_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xx_x[i] = 4.0 * tr_xxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_y[i] = -2.0 * tr_xxz_0[i] * tbe_0 + 4.0 * tr_xxz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_z[i] = 4.0 * tr_xxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 129-132 components of targeted buffer : DP

    auto tr_z_0_y_xy_x = pbuffer.data(idx_op_geom_110_dp + 129);

    auto tr_z_0_y_xy_y = pbuffer.data(idx_op_geom_110_dp + 130);

    auto tr_z_0_y_xy_z = pbuffer.data(idx_op_geom_110_dp + 131);

    #pragma omp simd aligned(tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyz_0, tr_xyz_xy, tr_xyz_yy, tr_xyz_yz, tr_xz_x, tr_xz_y, tr_xz_z, tr_z_0_y_xy_x, tr_z_0_y_xy_y, tr_z_0_y_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xy_x[i] = -2.0 * tr_xz_x[i] * tbe_0 + 4.0 * tr_xyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_y[i] = -2.0 * tr_xz_y[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_z[i] = -2.0 * tr_xz_z[i] * tbe_0 + 4.0 * tr_xyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 132-135 components of targeted buffer : DP

    auto tr_z_0_y_xz_x = pbuffer.data(idx_op_geom_110_dp + 132);

    auto tr_z_0_y_xz_y = pbuffer.data(idx_op_geom_110_dp + 133);

    auto tr_z_0_y_xz_z = pbuffer.data(idx_op_geom_110_dp + 134);

    #pragma omp simd aligned(tr_x_0, tr_x_xy, tr_x_yy, tr_x_yz, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xzz_0, tr_xzz_xy, tr_xzz_yy, tr_xzz_yz, tr_z_0_y_xz_x, tr_z_0_y_xz_y, tr_z_0_y_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xz_x[i] = -2.0 * tr_x_xy[i] * tke_0 + 4.0 * tr_xzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_x[i] * tbe_0 + 4.0 * tr_xyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_y[i] = tr_x_0[i] - 2.0 * tr_x_yy[i] * tke_0 - 2.0 * tr_xzz_0[i] * tbe_0 + 4.0 * tr_xzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_y[i] * tbe_0 + 4.0 * tr_xyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_z[i] = -2.0 * tr_x_yz[i] * tke_0 + 4.0 * tr_xzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_z[i] * tbe_0 + 4.0 * tr_xyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 135-138 components of targeted buffer : DP

    auto tr_z_0_y_yy_x = pbuffer.data(idx_op_geom_110_dp + 135);

    auto tr_z_0_y_yy_y = pbuffer.data(idx_op_geom_110_dp + 136);

    auto tr_z_0_y_yy_z = pbuffer.data(idx_op_geom_110_dp + 137);

    #pragma omp simd aligned(tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yyz_0, tr_yyz_xy, tr_yyz_yy, tr_yyz_yz, tr_yz_x, tr_yz_y, tr_yz_z, tr_z_0_y_yy_x, tr_z_0_y_yy_y, tr_z_0_y_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yy_x[i] = -4.0 * tr_yz_x[i] * tbe_0 + 4.0 * tr_yyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_y[i] = -4.0 * tr_yz_y[i] * tbe_0 - 2.0 * tr_yyz_0[i] * tbe_0 + 4.0 * tr_yyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_z[i] = -4.0 * tr_yz_z[i] * tbe_0 + 4.0 * tr_yyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 138-141 components of targeted buffer : DP

    auto tr_z_0_y_yz_x = pbuffer.data(idx_op_geom_110_dp + 138);

    auto tr_z_0_y_yz_y = pbuffer.data(idx_op_geom_110_dp + 139);

    auto tr_z_0_y_yz_z = pbuffer.data(idx_op_geom_110_dp + 140);

    #pragma omp simd aligned(tr_0_x, tr_0_y, tr_0_z, tr_y_0, tr_y_xy, tr_y_yy, tr_y_yz, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_yzz_0, tr_yzz_xy, tr_yzz_yy, tr_yzz_yz, tr_z_0_y_yz_x, tr_z_0_y_yz_y, tr_z_0_y_yz_z, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yz_x[i] = tr_0_x[i] - 2.0 * tr_zz_x[i] * tbe_0 - 2.0 * tr_y_xy[i] * tke_0 + 4.0 * tr_yzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_x[i] * tbe_0 + 4.0 * tr_yyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_y[i] = tr_0_y[i] - 2.0 * tr_zz_y[i] * tbe_0 + tr_y_0[i] - 2.0 * tr_y_yy[i] * tke_0 - 2.0 * tr_yzz_0[i] * tbe_0 + 4.0 * tr_yzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_y[i] * tbe_0 + 4.0 * tr_yyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_z[i] = tr_0_z[i] - 2.0 * tr_zz_z[i] * tbe_0 - 2.0 * tr_y_yz[i] * tke_0 + 4.0 * tr_yzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_z[i] * tbe_0 + 4.0 * tr_yyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 141-144 components of targeted buffer : DP

    auto tr_z_0_y_zz_x = pbuffer.data(idx_op_geom_110_dp + 141);

    auto tr_z_0_y_zz_y = pbuffer.data(idx_op_geom_110_dp + 142);

    auto tr_z_0_y_zz_z = pbuffer.data(idx_op_geom_110_dp + 143);

    #pragma omp simd aligned(tr_yz_x, tr_yz_y, tr_yz_z, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_z_0, tr_z_0_y_zz_x, tr_z_0_y_zz_y, tr_z_0_y_zz_z, tr_z_xy, tr_z_yy, tr_z_yz, tr_zzz_0, tr_zzz_xy, tr_zzz_yy, tr_zzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_zz_x[i] = -4.0 * tr_z_xy[i] * tke_0 + 4.0 * tr_zzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_x[i] * tbe_0 + 4.0 * tr_yzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_y[i] = 2.0 * tr_z_0[i] - 4.0 * tr_z_yy[i] * tke_0 - 2.0 * tr_zzz_0[i] * tbe_0 + 4.0 * tr_zzz_yy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_y[i] * tbe_0 + 4.0 * tr_yzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_z[i] = -4.0 * tr_z_yz[i] * tke_0 + 4.0 * tr_zzz_yz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_z[i] * tbe_0 + 4.0 * tr_yzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 144-147 components of targeted buffer : DP

    auto tr_z_0_z_xx_x = pbuffer.data(idx_op_geom_110_dp + 144);

    auto tr_z_0_z_xx_y = pbuffer.data(idx_op_geom_110_dp + 145);

    auto tr_z_0_z_xx_z = pbuffer.data(idx_op_geom_110_dp + 146);

    #pragma omp simd aligned(tr_xx_x, tr_xx_y, tr_xx_z, tr_xxz_0, tr_xxz_xz, tr_xxz_yz, tr_xxz_zz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_z_0_z_xx_x, tr_z_0_z_xx_y, tr_z_0_z_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xx_x[i] = -2.0 * tr_xx_x[i] * tbe_0 + 4.0 * tr_xxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_y[i] = -2.0 * tr_xx_y[i] * tbe_0 + 4.0 * tr_xxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_z[i] = -2.0 * tr_xx_z[i] * tbe_0 - 2.0 * tr_xxz_0[i] * tbe_0 + 4.0 * tr_xxz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 147-150 components of targeted buffer : DP

    auto tr_z_0_z_xy_x = pbuffer.data(idx_op_geom_110_dp + 147);

    auto tr_z_0_z_xy_y = pbuffer.data(idx_op_geom_110_dp + 148);

    auto tr_z_0_z_xy_z = pbuffer.data(idx_op_geom_110_dp + 149);

    #pragma omp simd aligned(tr_xy_x, tr_xy_y, tr_xy_z, tr_xyz_0, tr_xyz_xz, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_z_0_z_xy_x, tr_z_0_z_xy_y, tr_z_0_z_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xy_x[i] = -2.0 * tr_xy_x[i] * tbe_0 + 4.0 * tr_xyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_y[i] = -2.0 * tr_xy_y[i] * tbe_0 + 4.0 * tr_xyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_z[i] = -2.0 * tr_xy_z[i] * tbe_0 - 2.0 * tr_xyz_0[i] * tbe_0 + 4.0 * tr_xyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 150-153 components of targeted buffer : DP

    auto tr_z_0_z_xz_x = pbuffer.data(idx_op_geom_110_dp + 150);

    auto tr_z_0_z_xz_y = pbuffer.data(idx_op_geom_110_dp + 151);

    auto tr_z_0_z_xz_z = pbuffer.data(idx_op_geom_110_dp + 152);

    #pragma omp simd aligned(tr_x_0, tr_x_xz, tr_x_yz, tr_x_zz, tr_xz_x, tr_xz_y, tr_xz_z, tr_xzz_0, tr_xzz_xz, tr_xzz_yz, tr_xzz_zz, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_z_0_z_xz_x, tr_z_0_z_xz_y, tr_z_0_z_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xz_x[i] = -2.0 * tr_x_xz[i] * tke_0 - 6.0 * tr_xz_x[i] * tbe_0 + 4.0 * tr_xzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_y[i] = -2.0 * tr_x_yz[i] * tke_0 - 6.0 * tr_xz_y[i] * tbe_0 + 4.0 * tr_xzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_z[i] = tr_x_0[i] - 2.0 * tr_x_zz[i] * tke_0 - 6.0 * tr_xz_z[i] * tbe_0 - 2.0 * tr_xzz_0[i] * tbe_0 + 4.0 * tr_xzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 153-156 components of targeted buffer : DP

    auto tr_z_0_z_yy_x = pbuffer.data(idx_op_geom_110_dp + 153);

    auto tr_z_0_z_yy_y = pbuffer.data(idx_op_geom_110_dp + 154);

    auto tr_z_0_z_yy_z = pbuffer.data(idx_op_geom_110_dp + 155);

    #pragma omp simd aligned(tr_yy_x, tr_yy_y, tr_yy_z, tr_yyz_0, tr_yyz_xz, tr_yyz_yz, tr_yyz_zz, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_z_0_z_yy_x, tr_z_0_z_yy_y, tr_z_0_z_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yy_x[i] = -2.0 * tr_yy_x[i] * tbe_0 + 4.0 * tr_yyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_y[i] = -2.0 * tr_yy_y[i] * tbe_0 + 4.0 * tr_yyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_z[i] = -2.0 * tr_yy_z[i] * tbe_0 - 2.0 * tr_yyz_0[i] * tbe_0 + 4.0 * tr_yyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 156-159 components of targeted buffer : DP

    auto tr_z_0_z_yz_x = pbuffer.data(idx_op_geom_110_dp + 156);

    auto tr_z_0_z_yz_y = pbuffer.data(idx_op_geom_110_dp + 157);

    auto tr_z_0_z_yz_z = pbuffer.data(idx_op_geom_110_dp + 158);

    #pragma omp simd aligned(tr_y_0, tr_y_xz, tr_y_yz, tr_y_zz, tr_yz_x, tr_yz_y, tr_yz_z, tr_yzz_0, tr_yzz_xz, tr_yzz_yz, tr_yzz_zz, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_z_0_z_yz_x, tr_z_0_z_yz_y, tr_z_0_z_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yz_x[i] = -2.0 * tr_y_xz[i] * tke_0 - 6.0 * tr_yz_x[i] * tbe_0 + 4.0 * tr_yzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_y[i] = -2.0 * tr_y_yz[i] * tke_0 - 6.0 * tr_yz_y[i] * tbe_0 + 4.0 * tr_yzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_z[i] = tr_y_0[i] - 2.0 * tr_y_zz[i] * tke_0 - 6.0 * tr_yz_z[i] * tbe_0 - 2.0 * tr_yzz_0[i] * tbe_0 + 4.0 * tr_yzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 159-162 components of targeted buffer : DP

    auto tr_z_0_z_zz_x = pbuffer.data(idx_op_geom_110_dp + 159);

    auto tr_z_0_z_zz_y = pbuffer.data(idx_op_geom_110_dp + 160);

    auto tr_z_0_z_zz_z = pbuffer.data(idx_op_geom_110_dp + 161);

    #pragma omp simd aligned(tr_0_x, tr_0_y, tr_0_z, tr_z_0, tr_z_0_z_zz_x, tr_z_0_z_zz_y, tr_z_0_z_zz_z, tr_z_xz, tr_z_yz, tr_z_zz, tr_zz_x, tr_zz_y, tr_zz_z, tr_zzz_0, tr_zzz_xz, tr_zzz_yz, tr_zzz_zz, tr_zzzz_x, tr_zzzz_y, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_zz_x[i] = 2.0 * tr_0_x[i] - 4.0 * tr_z_xz[i] * tke_0 - 10.0 * tr_zz_x[i] * tbe_0 + 4.0 * tr_zzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_y[i] = 2.0 * tr_0_y[i] - 4.0 * tr_z_yz[i] * tke_0 - 10.0 * tr_zz_y[i] * tbe_0 + 4.0 * tr_zzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_z[i] = 2.0 * tr_0_z[i] + 2.0 * tr_z_0[i] - 4.0 * tr_z_zz[i] * tke_0 - 10.0 * tr_zz_z[i] * tbe_0 - 2.0 * tr_zzz_0[i] * tbe_0 + 4.0 * tr_zzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_z[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

