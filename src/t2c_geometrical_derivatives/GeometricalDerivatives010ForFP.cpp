#include "GeometricalDerivatives010ForFP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_fp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_fp,
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

    // Set up 0-3 components of targeted buffer : FP

    auto tr_0_0_x_xxx_x = pbuffer.data(idx_op_geom_010_fp);

    auto tr_0_0_x_xxx_y = pbuffer.data(idx_op_geom_010_fp + 1);

    auto tr_0_0_x_xxx_z = pbuffer.data(idx_op_geom_010_fp + 2);

    #pragma omp simd aligned(tr_0_0_x_xxx_x, tr_0_0_x_xxx_y, tr_0_0_x_xxx_z, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxx_0, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxxx_x, tr_xxxx_y, tr_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxx_x[i] = 2.0 * tr_xxxx_x[i] * tbe_0 + 2.0 * tr_xxx_xx[i] * tke_0 - 3.0 * tr_xx_x[i] - tr_xxx_0[i];

        tr_0_0_x_xxx_y[i] = 2.0 * tr_xxxx_y[i] * tbe_0 + 2.0 * tr_xxx_xy[i] * tke_0 - 3.0 * tr_xx_y[i];

        tr_0_0_x_xxx_z[i] = 2.0 * tr_xxxx_z[i] * tbe_0 + 2.0 * tr_xxx_xz[i] * tke_0 - 3.0 * tr_xx_z[i];
    }

    // Set up 3-6 components of targeted buffer : FP

    auto tr_0_0_x_xxy_x = pbuffer.data(idx_op_geom_010_fp + 3);

    auto tr_0_0_x_xxy_y = pbuffer.data(idx_op_geom_010_fp + 4);

    auto tr_0_0_x_xxy_z = pbuffer.data(idx_op_geom_010_fp + 5);

    #pragma omp simd aligned(tr_0_0_x_xxy_x, tr_0_0_x_xxy_y, tr_0_0_x_xxy_z, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, tr_xxy_0, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xy_x, tr_xy_y, tr_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxy_x[i] = 2.0 * tr_xxxy_x[i] * tbe_0 + 2.0 * tr_xxy_xx[i] * tke_0 - 2.0 * tr_xy_x[i] - tr_xxy_0[i];

        tr_0_0_x_xxy_y[i] = 2.0 * tr_xxxy_y[i] * tbe_0 + 2.0 * tr_xxy_xy[i] * tke_0 - 2.0 * tr_xy_y[i];

        tr_0_0_x_xxy_z[i] = 2.0 * tr_xxxy_z[i] * tbe_0 + 2.0 * tr_xxy_xz[i] * tke_0 - 2.0 * tr_xy_z[i];
    }

    // Set up 6-9 components of targeted buffer : FP

    auto tr_0_0_x_xxz_x = pbuffer.data(idx_op_geom_010_fp + 6);

    auto tr_0_0_x_xxz_y = pbuffer.data(idx_op_geom_010_fp + 7);

    auto tr_0_0_x_xxz_z = pbuffer.data(idx_op_geom_010_fp + 8);

    #pragma omp simd aligned(tr_0_0_x_xxz_x, tr_0_0_x_xxz_y, tr_0_0_x_xxz_z, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, tr_xxz_0, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xz_x, tr_xz_y, tr_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxz_x[i] = 2.0 * tr_xxxz_x[i] * tbe_0 + 2.0 * tr_xxz_xx[i] * tke_0 - 2.0 * tr_xz_x[i] - tr_xxz_0[i];

        tr_0_0_x_xxz_y[i] = 2.0 * tr_xxxz_y[i] * tbe_0 + 2.0 * tr_xxz_xy[i] * tke_0 - 2.0 * tr_xz_y[i];

        tr_0_0_x_xxz_z[i] = 2.0 * tr_xxxz_z[i] * tbe_0 + 2.0 * tr_xxz_xz[i] * tke_0 - 2.0 * tr_xz_z[i];
    }

    // Set up 9-12 components of targeted buffer : FP

    auto tr_0_0_x_xyy_x = pbuffer.data(idx_op_geom_010_fp + 9);

    auto tr_0_0_x_xyy_y = pbuffer.data(idx_op_geom_010_fp + 10);

    auto tr_0_0_x_xyy_z = pbuffer.data(idx_op_geom_010_fp + 11);

    #pragma omp simd aligned(tr_0_0_x_xyy_x, tr_0_0_x_xyy_y, tr_0_0_x_xyy_z, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, tr_xyy_0, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_yy_x, tr_yy_y, tr_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyy_x[i] = 2.0 * tr_xxyy_x[i] * tbe_0 + 2.0 * tr_xyy_xx[i] * tke_0 - tr_yy_x[i] - tr_xyy_0[i];

        tr_0_0_x_xyy_y[i] = 2.0 * tr_xxyy_y[i] * tbe_0 + 2.0 * tr_xyy_xy[i] * tke_0 - tr_yy_y[i];

        tr_0_0_x_xyy_z[i] = 2.0 * tr_xxyy_z[i] * tbe_0 + 2.0 * tr_xyy_xz[i] * tke_0 - tr_yy_z[i];
    }

    // Set up 12-15 components of targeted buffer : FP

    auto tr_0_0_x_xyz_x = pbuffer.data(idx_op_geom_010_fp + 12);

    auto tr_0_0_x_xyz_y = pbuffer.data(idx_op_geom_010_fp + 13);

    auto tr_0_0_x_xyz_z = pbuffer.data(idx_op_geom_010_fp + 14);

    #pragma omp simd aligned(tr_0_0_x_xyz_x, tr_0_0_x_xyz_y, tr_0_0_x_xyz_z, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xyz_0, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_yz_x, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyz_x[i] = 2.0 * tr_xxyz_x[i] * tbe_0 + 2.0 * tr_xyz_xx[i] * tke_0 - tr_yz_x[i] - tr_xyz_0[i];

        tr_0_0_x_xyz_y[i] = 2.0 * tr_xxyz_y[i] * tbe_0 + 2.0 * tr_xyz_xy[i] * tke_0 - tr_yz_y[i];

        tr_0_0_x_xyz_z[i] = 2.0 * tr_xxyz_z[i] * tbe_0 + 2.0 * tr_xyz_xz[i] * tke_0 - tr_yz_z[i];
    }

    // Set up 15-18 components of targeted buffer : FP

    auto tr_0_0_x_xzz_x = pbuffer.data(idx_op_geom_010_fp + 15);

    auto tr_0_0_x_xzz_y = pbuffer.data(idx_op_geom_010_fp + 16);

    auto tr_0_0_x_xzz_z = pbuffer.data(idx_op_geom_010_fp + 17);

    #pragma omp simd aligned(tr_0_0_x_xzz_x, tr_0_0_x_xzz_y, tr_0_0_x_xzz_z, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, tr_xzz_0, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xzz_x[i] = 2.0 * tr_xxzz_x[i] * tbe_0 + 2.0 * tr_xzz_xx[i] * tke_0 - tr_zz_x[i] - tr_xzz_0[i];

        tr_0_0_x_xzz_y[i] = 2.0 * tr_xxzz_y[i] * tbe_0 + 2.0 * tr_xzz_xy[i] * tke_0 - tr_zz_y[i];

        tr_0_0_x_xzz_z[i] = 2.0 * tr_xxzz_z[i] * tbe_0 + 2.0 * tr_xzz_xz[i] * tke_0 - tr_zz_z[i];
    }

    // Set up 18-21 components of targeted buffer : FP

    auto tr_0_0_x_yyy_x = pbuffer.data(idx_op_geom_010_fp + 18);

    auto tr_0_0_x_yyy_y = pbuffer.data(idx_op_geom_010_fp + 19);

    auto tr_0_0_x_yyy_z = pbuffer.data(idx_op_geom_010_fp + 20);

    #pragma omp simd aligned(tr_0_0_x_yyy_x, tr_0_0_x_yyy_y, tr_0_0_x_yyy_z, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, tr_yyy_0, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyy_x[i] = 2.0 * tr_xyyy_x[i] * tbe_0 + 2.0 * tr_yyy_xx[i] * tke_0 - tr_yyy_0[i];

        tr_0_0_x_yyy_y[i] = 2.0 * tr_xyyy_y[i] * tbe_0 + 2.0 * tr_yyy_xy[i] * tke_0;

        tr_0_0_x_yyy_z[i] = 2.0 * tr_xyyy_z[i] * tbe_0 + 2.0 * tr_yyy_xz[i] * tke_0;
    }

    // Set up 21-24 components of targeted buffer : FP

    auto tr_0_0_x_yyz_x = pbuffer.data(idx_op_geom_010_fp + 21);

    auto tr_0_0_x_yyz_y = pbuffer.data(idx_op_geom_010_fp + 22);

    auto tr_0_0_x_yyz_z = pbuffer.data(idx_op_geom_010_fp + 23);

    #pragma omp simd aligned(tr_0_0_x_yyz_x, tr_0_0_x_yyz_y, tr_0_0_x_yyz_z, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_yyz_0, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyz_x[i] = 2.0 * tr_xyyz_x[i] * tbe_0 + 2.0 * tr_yyz_xx[i] * tke_0 - tr_yyz_0[i];

        tr_0_0_x_yyz_y[i] = 2.0 * tr_xyyz_y[i] * tbe_0 + 2.0 * tr_yyz_xy[i] * tke_0;

        tr_0_0_x_yyz_z[i] = 2.0 * tr_xyyz_z[i] * tbe_0 + 2.0 * tr_yyz_xz[i] * tke_0;
    }

    // Set up 24-27 components of targeted buffer : FP

    auto tr_0_0_x_yzz_x = pbuffer.data(idx_op_geom_010_fp + 24);

    auto tr_0_0_x_yzz_y = pbuffer.data(idx_op_geom_010_fp + 25);

    auto tr_0_0_x_yzz_z = pbuffer.data(idx_op_geom_010_fp + 26);

    #pragma omp simd aligned(tr_0_0_x_yzz_x, tr_0_0_x_yzz_y, tr_0_0_x_yzz_z, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_yzz_0, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yzz_x[i] = 2.0 * tr_xyzz_x[i] * tbe_0 + 2.0 * tr_yzz_xx[i] * tke_0 - tr_yzz_0[i];

        tr_0_0_x_yzz_y[i] = 2.0 * tr_xyzz_y[i] * tbe_0 + 2.0 * tr_yzz_xy[i] * tke_0;

        tr_0_0_x_yzz_z[i] = 2.0 * tr_xyzz_z[i] * tbe_0 + 2.0 * tr_yzz_xz[i] * tke_0;
    }

    // Set up 27-30 components of targeted buffer : FP

    auto tr_0_0_x_zzz_x = pbuffer.data(idx_op_geom_010_fp + 27);

    auto tr_0_0_x_zzz_y = pbuffer.data(idx_op_geom_010_fp + 28);

    auto tr_0_0_x_zzz_z = pbuffer.data(idx_op_geom_010_fp + 29);

    #pragma omp simd aligned(tr_0_0_x_zzz_x, tr_0_0_x_zzz_y, tr_0_0_x_zzz_z, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, tr_zzz_0, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_zzz_x[i] = 2.0 * tr_xzzz_x[i] * tbe_0 + 2.0 * tr_zzz_xx[i] * tke_0 - tr_zzz_0[i];

        tr_0_0_x_zzz_y[i] = 2.0 * tr_xzzz_y[i] * tbe_0 + 2.0 * tr_zzz_xy[i] * tke_0;

        tr_0_0_x_zzz_z[i] = 2.0 * tr_xzzz_z[i] * tbe_0 + 2.0 * tr_zzz_xz[i] * tke_0;
    }

    // Set up 30-33 components of targeted buffer : FP

    auto tr_0_0_y_xxx_x = pbuffer.data(idx_op_geom_010_fp + 30);

    auto tr_0_0_y_xxx_y = pbuffer.data(idx_op_geom_010_fp + 31);

    auto tr_0_0_y_xxx_z = pbuffer.data(idx_op_geom_010_fp + 32);

    #pragma omp simd aligned(tr_0_0_y_xxx_x, tr_0_0_y_xxx_y, tr_0_0_y_xxx_z, tr_xxx_0, tr_xxx_xy, tr_xxx_yy, tr_xxx_yz, tr_xxxy_x, tr_xxxy_y, tr_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxx_x[i] = 2.0 * tr_xxxy_x[i] * tbe_0 + 2.0 * tr_xxx_xy[i] * tke_0;

        tr_0_0_y_xxx_y[i] = 2.0 * tr_xxxy_y[i] * tbe_0 + 2.0 * tr_xxx_yy[i] * tke_0 - tr_xxx_0[i];

        tr_0_0_y_xxx_z[i] = 2.0 * tr_xxxy_z[i] * tbe_0 + 2.0 * tr_xxx_yz[i] * tke_0;
    }

    // Set up 33-36 components of targeted buffer : FP

    auto tr_0_0_y_xxy_x = pbuffer.data(idx_op_geom_010_fp + 33);

    auto tr_0_0_y_xxy_y = pbuffer.data(idx_op_geom_010_fp + 34);

    auto tr_0_0_y_xxy_z = pbuffer.data(idx_op_geom_010_fp + 35);

    #pragma omp simd aligned(tr_0_0_y_xxy_x, tr_0_0_y_xxy_y, tr_0_0_y_xxy_z, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxy_0, tr_xxy_xy, tr_xxy_yy, tr_xxy_yz, tr_xxyy_x, tr_xxyy_y, tr_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxy_x[i] = 2.0 * tr_xxyy_x[i] * tbe_0 + 2.0 * tr_xxy_xy[i] * tke_0 - tr_xx_x[i];

        tr_0_0_y_xxy_y[i] = 2.0 * tr_xxyy_y[i] * tbe_0 + 2.0 * tr_xxy_yy[i] * tke_0 - tr_xx_y[i] - tr_xxy_0[i];

        tr_0_0_y_xxy_z[i] = 2.0 * tr_xxyy_z[i] * tbe_0 + 2.0 * tr_xxy_yz[i] * tke_0 - tr_xx_z[i];
    }

    // Set up 36-39 components of targeted buffer : FP

    auto tr_0_0_y_xxz_x = pbuffer.data(idx_op_geom_010_fp + 36);

    auto tr_0_0_y_xxz_y = pbuffer.data(idx_op_geom_010_fp + 37);

    auto tr_0_0_y_xxz_z = pbuffer.data(idx_op_geom_010_fp + 38);

    #pragma omp simd aligned(tr_0_0_y_xxz_x, tr_0_0_y_xxz_y, tr_0_0_y_xxz_z, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, tr_xxz_0, tr_xxz_xy, tr_xxz_yy, tr_xxz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxz_x[i] = 2.0 * tr_xxyz_x[i] * tbe_0 + 2.0 * tr_xxz_xy[i] * tke_0;

        tr_0_0_y_xxz_y[i] = 2.0 * tr_xxyz_y[i] * tbe_0 + 2.0 * tr_xxz_yy[i] * tke_0 - tr_xxz_0[i];

        tr_0_0_y_xxz_z[i] = 2.0 * tr_xxyz_z[i] * tbe_0 + 2.0 * tr_xxz_yz[i] * tke_0;
    }

    // Set up 39-42 components of targeted buffer : FP

    auto tr_0_0_y_xyy_x = pbuffer.data(idx_op_geom_010_fp + 39);

    auto tr_0_0_y_xyy_y = pbuffer.data(idx_op_geom_010_fp + 40);

    auto tr_0_0_y_xyy_z = pbuffer.data(idx_op_geom_010_fp + 41);

    #pragma omp simd aligned(tr_0_0_y_xyy_x, tr_0_0_y_xyy_y, tr_0_0_y_xyy_z, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyy_0, tr_xyy_xy, tr_xyy_yy, tr_xyy_yz, tr_xyyy_x, tr_xyyy_y, tr_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyy_x[i] = 2.0 * tr_xyyy_x[i] * tbe_0 + 2.0 * tr_xyy_xy[i] * tke_0 - 2.0 * tr_xy_x[i];

        tr_0_0_y_xyy_y[i] = 2.0 * tr_xyyy_y[i] * tbe_0 + 2.0 * tr_xyy_yy[i] * tke_0 - 2.0 * tr_xy_y[i] - tr_xyy_0[i];

        tr_0_0_y_xyy_z[i] = 2.0 * tr_xyyy_z[i] * tbe_0 + 2.0 * tr_xyy_yz[i] * tke_0 - 2.0 * tr_xy_z[i];
    }

    // Set up 42-45 components of targeted buffer : FP

    auto tr_0_0_y_xyz_x = pbuffer.data(idx_op_geom_010_fp + 42);

    auto tr_0_0_y_xyz_y = pbuffer.data(idx_op_geom_010_fp + 43);

    auto tr_0_0_y_xyz_z = pbuffer.data(idx_op_geom_010_fp + 44);

    #pragma omp simd aligned(tr_0_0_y_xyz_x, tr_0_0_y_xyz_y, tr_0_0_y_xyz_z, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, tr_xyz_0, tr_xyz_xy, tr_xyz_yy, tr_xyz_yz, tr_xz_x, tr_xz_y, tr_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyz_x[i] = 2.0 * tr_xyyz_x[i] * tbe_0 + 2.0 * tr_xyz_xy[i] * tke_0 - tr_xz_x[i];

        tr_0_0_y_xyz_y[i] = 2.0 * tr_xyyz_y[i] * tbe_0 + 2.0 * tr_xyz_yy[i] * tke_0 - tr_xz_y[i] - tr_xyz_0[i];

        tr_0_0_y_xyz_z[i] = 2.0 * tr_xyyz_z[i] * tbe_0 + 2.0 * tr_xyz_yz[i] * tke_0 - tr_xz_z[i];
    }

    // Set up 45-48 components of targeted buffer : FP

    auto tr_0_0_y_xzz_x = pbuffer.data(idx_op_geom_010_fp + 45);

    auto tr_0_0_y_xzz_y = pbuffer.data(idx_op_geom_010_fp + 46);

    auto tr_0_0_y_xzz_z = pbuffer.data(idx_op_geom_010_fp + 47);

    #pragma omp simd aligned(tr_0_0_y_xzz_x, tr_0_0_y_xzz_y, tr_0_0_y_xzz_z, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, tr_xzz_0, tr_xzz_xy, tr_xzz_yy, tr_xzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xzz_x[i] = 2.0 * tr_xyzz_x[i] * tbe_0 + 2.0 * tr_xzz_xy[i] * tke_0;

        tr_0_0_y_xzz_y[i] = 2.0 * tr_xyzz_y[i] * tbe_0 + 2.0 * tr_xzz_yy[i] * tke_0 - tr_xzz_0[i];

        tr_0_0_y_xzz_z[i] = 2.0 * tr_xyzz_z[i] * tbe_0 + 2.0 * tr_xzz_yz[i] * tke_0;
    }

    // Set up 48-51 components of targeted buffer : FP

    auto tr_0_0_y_yyy_x = pbuffer.data(idx_op_geom_010_fp + 48);

    auto tr_0_0_y_yyy_y = pbuffer.data(idx_op_geom_010_fp + 49);

    auto tr_0_0_y_yyy_z = pbuffer.data(idx_op_geom_010_fp + 50);

    #pragma omp simd aligned(tr_0_0_y_yyy_x, tr_0_0_y_yyy_y, tr_0_0_y_yyy_z, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyy_0, tr_yyy_xy, tr_yyy_yy, tr_yyy_yz, tr_yyyy_x, tr_yyyy_y, tr_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyy_x[i] = 2.0 * tr_yyyy_x[i] * tbe_0 + 2.0 * tr_yyy_xy[i] * tke_0 - 3.0 * tr_yy_x[i];

        tr_0_0_y_yyy_y[i] = 2.0 * tr_yyyy_y[i] * tbe_0 + 2.0 * tr_yyy_yy[i] * tke_0 - 3.0 * tr_yy_y[i] - tr_yyy_0[i];

        tr_0_0_y_yyy_z[i] = 2.0 * tr_yyyy_z[i] * tbe_0 + 2.0 * tr_yyy_yz[i] * tke_0 - 3.0 * tr_yy_z[i];
    }

    // Set up 51-54 components of targeted buffer : FP

    auto tr_0_0_y_yyz_x = pbuffer.data(idx_op_geom_010_fp + 51);

    auto tr_0_0_y_yyz_y = pbuffer.data(idx_op_geom_010_fp + 52);

    auto tr_0_0_y_yyz_z = pbuffer.data(idx_op_geom_010_fp + 53);

    #pragma omp simd aligned(tr_0_0_y_yyz_x, tr_0_0_y_yyz_y, tr_0_0_y_yyz_z, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, tr_yyz_0, tr_yyz_xy, tr_yyz_yy, tr_yyz_yz, tr_yz_x, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyz_x[i] = 2.0 * tr_yyyz_x[i] * tbe_0 + 2.0 * tr_yyz_xy[i] * tke_0 - 2.0 * tr_yz_x[i];

        tr_0_0_y_yyz_y[i] = 2.0 * tr_yyyz_y[i] * tbe_0 + 2.0 * tr_yyz_yy[i] * tke_0 - 2.0 * tr_yz_y[i] - tr_yyz_0[i];

        tr_0_0_y_yyz_z[i] = 2.0 * tr_yyyz_z[i] * tbe_0 + 2.0 * tr_yyz_yz[i] * tke_0 - 2.0 * tr_yz_z[i];
    }

    // Set up 54-57 components of targeted buffer : FP

    auto tr_0_0_y_yzz_x = pbuffer.data(idx_op_geom_010_fp + 54);

    auto tr_0_0_y_yzz_y = pbuffer.data(idx_op_geom_010_fp + 55);

    auto tr_0_0_y_yzz_z = pbuffer.data(idx_op_geom_010_fp + 56);

    #pragma omp simd aligned(tr_0_0_y_yzz_x, tr_0_0_y_yzz_y, tr_0_0_y_yzz_z, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, tr_yzz_0, tr_yzz_xy, tr_yzz_yy, tr_yzz_yz, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yzz_x[i] = 2.0 * tr_yyzz_x[i] * tbe_0 + 2.0 * tr_yzz_xy[i] * tke_0 - tr_zz_x[i];

        tr_0_0_y_yzz_y[i] = 2.0 * tr_yyzz_y[i] * tbe_0 + 2.0 * tr_yzz_yy[i] * tke_0 - tr_zz_y[i] - tr_yzz_0[i];

        tr_0_0_y_yzz_z[i] = 2.0 * tr_yyzz_z[i] * tbe_0 + 2.0 * tr_yzz_yz[i] * tke_0 - tr_zz_z[i];
    }

    // Set up 57-60 components of targeted buffer : FP

    auto tr_0_0_y_zzz_x = pbuffer.data(idx_op_geom_010_fp + 57);

    auto tr_0_0_y_zzz_y = pbuffer.data(idx_op_geom_010_fp + 58);

    auto tr_0_0_y_zzz_z = pbuffer.data(idx_op_geom_010_fp + 59);

    #pragma omp simd aligned(tr_0_0_y_zzz_x, tr_0_0_y_zzz_y, tr_0_0_y_zzz_z, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, tr_zzz_0, tr_zzz_xy, tr_zzz_yy, tr_zzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_zzz_x[i] = 2.0 * tr_yzzz_x[i] * tbe_0 + 2.0 * tr_zzz_xy[i] * tke_0;

        tr_0_0_y_zzz_y[i] = 2.0 * tr_yzzz_y[i] * tbe_0 + 2.0 * tr_zzz_yy[i] * tke_0 - tr_zzz_0[i];

        tr_0_0_y_zzz_z[i] = 2.0 * tr_yzzz_z[i] * tbe_0 + 2.0 * tr_zzz_yz[i] * tke_0;
    }

    // Set up 60-63 components of targeted buffer : FP

    auto tr_0_0_z_xxx_x = pbuffer.data(idx_op_geom_010_fp + 60);

    auto tr_0_0_z_xxx_y = pbuffer.data(idx_op_geom_010_fp + 61);

    auto tr_0_0_z_xxx_z = pbuffer.data(idx_op_geom_010_fp + 62);

    #pragma omp simd aligned(tr_0_0_z_xxx_x, tr_0_0_z_xxx_y, tr_0_0_z_xxx_z, tr_xxx_0, tr_xxx_xz, tr_xxx_yz, tr_xxx_zz, tr_xxxz_x, tr_xxxz_y, tr_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxx_x[i] = 2.0 * tr_xxxz_x[i] * tbe_0 + 2.0 * tr_xxx_xz[i] * tke_0;

        tr_0_0_z_xxx_y[i] = 2.0 * tr_xxxz_y[i] * tbe_0 + 2.0 * tr_xxx_yz[i] * tke_0;

        tr_0_0_z_xxx_z[i] = 2.0 * tr_xxxz_z[i] * tbe_0 + 2.0 * tr_xxx_zz[i] * tke_0 - tr_xxx_0[i];
    }

    // Set up 63-66 components of targeted buffer : FP

    auto tr_0_0_z_xxy_x = pbuffer.data(idx_op_geom_010_fp + 63);

    auto tr_0_0_z_xxy_y = pbuffer.data(idx_op_geom_010_fp + 64);

    auto tr_0_0_z_xxy_z = pbuffer.data(idx_op_geom_010_fp + 65);

    #pragma omp simd aligned(tr_0_0_z_xxy_x, tr_0_0_z_xxy_y, tr_0_0_z_xxy_z, tr_xxy_0, tr_xxy_xz, tr_xxy_yz, tr_xxy_zz, tr_xxyz_x, tr_xxyz_y, tr_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxy_x[i] = 2.0 * tr_xxyz_x[i] * tbe_0 + 2.0 * tr_xxy_xz[i] * tke_0;

        tr_0_0_z_xxy_y[i] = 2.0 * tr_xxyz_y[i] * tbe_0 + 2.0 * tr_xxy_yz[i] * tke_0;

        tr_0_0_z_xxy_z[i] = 2.0 * tr_xxyz_z[i] * tbe_0 + 2.0 * tr_xxy_zz[i] * tke_0 - tr_xxy_0[i];
    }

    // Set up 66-69 components of targeted buffer : FP

    auto tr_0_0_z_xxz_x = pbuffer.data(idx_op_geom_010_fp + 66);

    auto tr_0_0_z_xxz_y = pbuffer.data(idx_op_geom_010_fp + 67);

    auto tr_0_0_z_xxz_z = pbuffer.data(idx_op_geom_010_fp + 68);

    #pragma omp simd aligned(tr_0_0_z_xxz_x, tr_0_0_z_xxz_y, tr_0_0_z_xxz_z, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxz_0, tr_xxz_xz, tr_xxz_yz, tr_xxz_zz, tr_xxzz_x, tr_xxzz_y, tr_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxz_x[i] = 2.0 * tr_xxzz_x[i] * tbe_0 + 2.0 * tr_xxz_xz[i] * tke_0 - tr_xx_x[i];

        tr_0_0_z_xxz_y[i] = 2.0 * tr_xxzz_y[i] * tbe_0 + 2.0 * tr_xxz_yz[i] * tke_0 - tr_xx_y[i];

        tr_0_0_z_xxz_z[i] = 2.0 * tr_xxzz_z[i] * tbe_0 + 2.0 * tr_xxz_zz[i] * tke_0 - tr_xx_z[i] - tr_xxz_0[i];
    }

    // Set up 69-72 components of targeted buffer : FP

    auto tr_0_0_z_xyy_x = pbuffer.data(idx_op_geom_010_fp + 69);

    auto tr_0_0_z_xyy_y = pbuffer.data(idx_op_geom_010_fp + 70);

    auto tr_0_0_z_xyy_z = pbuffer.data(idx_op_geom_010_fp + 71);

    #pragma omp simd aligned(tr_0_0_z_xyy_x, tr_0_0_z_xyy_y, tr_0_0_z_xyy_z, tr_xyy_0, tr_xyy_xz, tr_xyy_yz, tr_xyy_zz, tr_xyyz_x, tr_xyyz_y, tr_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyy_x[i] = 2.0 * tr_xyyz_x[i] * tbe_0 + 2.0 * tr_xyy_xz[i] * tke_0;

        tr_0_0_z_xyy_y[i] = 2.0 * tr_xyyz_y[i] * tbe_0 + 2.0 * tr_xyy_yz[i] * tke_0;

        tr_0_0_z_xyy_z[i] = 2.0 * tr_xyyz_z[i] * tbe_0 + 2.0 * tr_xyy_zz[i] * tke_0 - tr_xyy_0[i];
    }

    // Set up 72-75 components of targeted buffer : FP

    auto tr_0_0_z_xyz_x = pbuffer.data(idx_op_geom_010_fp + 72);

    auto tr_0_0_z_xyz_y = pbuffer.data(idx_op_geom_010_fp + 73);

    auto tr_0_0_z_xyz_z = pbuffer.data(idx_op_geom_010_fp + 74);

    #pragma omp simd aligned(tr_0_0_z_xyz_x, tr_0_0_z_xyz_y, tr_0_0_z_xyz_z, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyz_0, tr_xyz_xz, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_y, tr_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyz_x[i] = 2.0 * tr_xyzz_x[i] * tbe_0 + 2.0 * tr_xyz_xz[i] * tke_0 - tr_xy_x[i];

        tr_0_0_z_xyz_y[i] = 2.0 * tr_xyzz_y[i] * tbe_0 + 2.0 * tr_xyz_yz[i] * tke_0 - tr_xy_y[i];

        tr_0_0_z_xyz_z[i] = 2.0 * tr_xyzz_z[i] * tbe_0 + 2.0 * tr_xyz_zz[i] * tke_0 - tr_xy_z[i] - tr_xyz_0[i];
    }

    // Set up 75-78 components of targeted buffer : FP

    auto tr_0_0_z_xzz_x = pbuffer.data(idx_op_geom_010_fp + 75);

    auto tr_0_0_z_xzz_y = pbuffer.data(idx_op_geom_010_fp + 76);

    auto tr_0_0_z_xzz_z = pbuffer.data(idx_op_geom_010_fp + 77);

    #pragma omp simd aligned(tr_0_0_z_xzz_x, tr_0_0_z_xzz_y, tr_0_0_z_xzz_z, tr_xz_x, tr_xz_y, tr_xz_z, tr_xzz_0, tr_xzz_xz, tr_xzz_yz, tr_xzz_zz, tr_xzzz_x, tr_xzzz_y, tr_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xzz_x[i] = 2.0 * tr_xzzz_x[i] * tbe_0 + 2.0 * tr_xzz_xz[i] * tke_0 - 2.0 * tr_xz_x[i];

        tr_0_0_z_xzz_y[i] = 2.0 * tr_xzzz_y[i] * tbe_0 + 2.0 * tr_xzz_yz[i] * tke_0 - 2.0 * tr_xz_y[i];

        tr_0_0_z_xzz_z[i] = 2.0 * tr_xzzz_z[i] * tbe_0 + 2.0 * tr_xzz_zz[i] * tke_0 - 2.0 * tr_xz_z[i] - tr_xzz_0[i];
    }

    // Set up 78-81 components of targeted buffer : FP

    auto tr_0_0_z_yyy_x = pbuffer.data(idx_op_geom_010_fp + 78);

    auto tr_0_0_z_yyy_y = pbuffer.data(idx_op_geom_010_fp + 79);

    auto tr_0_0_z_yyy_z = pbuffer.data(idx_op_geom_010_fp + 80);

    #pragma omp simd aligned(tr_0_0_z_yyy_x, tr_0_0_z_yyy_y, tr_0_0_z_yyy_z, tr_yyy_0, tr_yyy_xz, tr_yyy_yz, tr_yyy_zz, tr_yyyz_x, tr_yyyz_y, tr_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyy_x[i] = 2.0 * tr_yyyz_x[i] * tbe_0 + 2.0 * tr_yyy_xz[i] * tke_0;

        tr_0_0_z_yyy_y[i] = 2.0 * tr_yyyz_y[i] * tbe_0 + 2.0 * tr_yyy_yz[i] * tke_0;

        tr_0_0_z_yyy_z[i] = 2.0 * tr_yyyz_z[i] * tbe_0 + 2.0 * tr_yyy_zz[i] * tke_0 - tr_yyy_0[i];
    }

    // Set up 81-84 components of targeted buffer : FP

    auto tr_0_0_z_yyz_x = pbuffer.data(idx_op_geom_010_fp + 81);

    auto tr_0_0_z_yyz_y = pbuffer.data(idx_op_geom_010_fp + 82);

    auto tr_0_0_z_yyz_z = pbuffer.data(idx_op_geom_010_fp + 83);

    #pragma omp simd aligned(tr_0_0_z_yyz_x, tr_0_0_z_yyz_y, tr_0_0_z_yyz_z, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyz_0, tr_yyz_xz, tr_yyz_yz, tr_yyz_zz, tr_yyzz_x, tr_yyzz_y, tr_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyz_x[i] = 2.0 * tr_yyzz_x[i] * tbe_0 + 2.0 * tr_yyz_xz[i] * tke_0 - tr_yy_x[i];

        tr_0_0_z_yyz_y[i] = 2.0 * tr_yyzz_y[i] * tbe_0 + 2.0 * tr_yyz_yz[i] * tke_0 - tr_yy_y[i];

        tr_0_0_z_yyz_z[i] = 2.0 * tr_yyzz_z[i] * tbe_0 + 2.0 * tr_yyz_zz[i] * tke_0 - tr_yy_z[i] - tr_yyz_0[i];
    }

    // Set up 84-87 components of targeted buffer : FP

    auto tr_0_0_z_yzz_x = pbuffer.data(idx_op_geom_010_fp + 84);

    auto tr_0_0_z_yzz_y = pbuffer.data(idx_op_geom_010_fp + 85);

    auto tr_0_0_z_yzz_z = pbuffer.data(idx_op_geom_010_fp + 86);

    #pragma omp simd aligned(tr_0_0_z_yzz_x, tr_0_0_z_yzz_y, tr_0_0_z_yzz_z, tr_yz_x, tr_yz_y, tr_yz_z, tr_yzz_0, tr_yzz_xz, tr_yzz_yz, tr_yzz_zz, tr_yzzz_x, tr_yzzz_y, tr_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yzz_x[i] = 2.0 * tr_yzzz_x[i] * tbe_0 + 2.0 * tr_yzz_xz[i] * tke_0 - 2.0 * tr_yz_x[i];

        tr_0_0_z_yzz_y[i] = 2.0 * tr_yzzz_y[i] * tbe_0 + 2.0 * tr_yzz_yz[i] * tke_0 - 2.0 * tr_yz_y[i];

        tr_0_0_z_yzz_z[i] = 2.0 * tr_yzzz_z[i] * tbe_0 + 2.0 * tr_yzz_zz[i] * tke_0 - 2.0 * tr_yz_z[i] - tr_yzz_0[i];
    }

    // Set up 87-90 components of targeted buffer : FP

    auto tr_0_0_z_zzz_x = pbuffer.data(idx_op_geom_010_fp + 87);

    auto tr_0_0_z_zzz_y = pbuffer.data(idx_op_geom_010_fp + 88);

    auto tr_0_0_z_zzz_z = pbuffer.data(idx_op_geom_010_fp + 89);

    #pragma omp simd aligned(tr_0_0_z_zzz_x, tr_0_0_z_zzz_y, tr_0_0_z_zzz_z, tr_zz_x, tr_zz_y, tr_zz_z, tr_zzz_0, tr_zzz_xz, tr_zzz_yz, tr_zzz_zz, tr_zzzz_x, tr_zzzz_y, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_zzz_x[i] = 2.0 * tr_zzzz_x[i] * tbe_0 + 2.0 * tr_zzz_xz[i] * tke_0 - 3.0 * tr_zz_x[i];

        tr_0_0_z_zzz_y[i] = 2.0 * tr_zzzz_y[i] * tbe_0 + 2.0 * tr_zzz_yz[i] * tke_0 - 3.0 * tr_zz_y[i];

        tr_0_0_z_zzz_z[i] = 2.0 * tr_zzzz_z[i] * tbe_0 + 2.0 * tr_zzz_zz[i] * tke_0 - 3.0 * tr_zz_z[i] - tr_zzz_0[i];
    }

}

} // t2cgeom namespace

