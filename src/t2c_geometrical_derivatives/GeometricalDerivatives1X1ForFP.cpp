#include "GeometricalDerivatives1X1ForFP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_fp(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_fp,
                        const size_t idx_op_ds,
                        const size_t idx_op_dd,
                        const size_t idx_op_gs,
                        const size_t idx_op_gd,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : DS

        auto to_xx_0 = pbuffer.data(idx_op_ds + i * 6 + 0);

        auto to_xy_0 = pbuffer.data(idx_op_ds + i * 6 + 1);

        auto to_xz_0 = pbuffer.data(idx_op_ds + i * 6 + 2);

        auto to_yy_0 = pbuffer.data(idx_op_ds + i * 6 + 3);

        auto to_yz_0 = pbuffer.data(idx_op_ds + i * 6 + 4);

        auto to_zz_0 = pbuffer.data(idx_op_ds + i * 6 + 5);

        // Set up components of auxiliary buffer : DD

        auto to_xx_xx = pbuffer.data(idx_op_dd + i * 36 + 0);

        auto to_xx_xy = pbuffer.data(idx_op_dd + i * 36 + 1);

        auto to_xx_xz = pbuffer.data(idx_op_dd + i * 36 + 2);

        auto to_xx_yy = pbuffer.data(idx_op_dd + i * 36 + 3);

        auto to_xx_yz = pbuffer.data(idx_op_dd + i * 36 + 4);

        auto to_xx_zz = pbuffer.data(idx_op_dd + i * 36 + 5);

        auto to_xy_xx = pbuffer.data(idx_op_dd + i * 36 + 6);

        auto to_xy_xy = pbuffer.data(idx_op_dd + i * 36 + 7);

        auto to_xy_xz = pbuffer.data(idx_op_dd + i * 36 + 8);

        auto to_xy_yy = pbuffer.data(idx_op_dd + i * 36 + 9);

        auto to_xy_yz = pbuffer.data(idx_op_dd + i * 36 + 10);

        auto to_xy_zz = pbuffer.data(idx_op_dd + i * 36 + 11);

        auto to_xz_xx = pbuffer.data(idx_op_dd + i * 36 + 12);

        auto to_xz_xy = pbuffer.data(idx_op_dd + i * 36 + 13);

        auto to_xz_xz = pbuffer.data(idx_op_dd + i * 36 + 14);

        auto to_xz_yy = pbuffer.data(idx_op_dd + i * 36 + 15);

        auto to_xz_yz = pbuffer.data(idx_op_dd + i * 36 + 16);

        auto to_xz_zz = pbuffer.data(idx_op_dd + i * 36 + 17);

        auto to_yy_xx = pbuffer.data(idx_op_dd + i * 36 + 18);

        auto to_yy_xy = pbuffer.data(idx_op_dd + i * 36 + 19);

        auto to_yy_xz = pbuffer.data(idx_op_dd + i * 36 + 20);

        auto to_yy_yy = pbuffer.data(idx_op_dd + i * 36 + 21);

        auto to_yy_yz = pbuffer.data(idx_op_dd + i * 36 + 22);

        auto to_yy_zz = pbuffer.data(idx_op_dd + i * 36 + 23);

        auto to_yz_xx = pbuffer.data(idx_op_dd + i * 36 + 24);

        auto to_yz_xy = pbuffer.data(idx_op_dd + i * 36 + 25);

        auto to_yz_xz = pbuffer.data(idx_op_dd + i * 36 + 26);

        auto to_yz_yy = pbuffer.data(idx_op_dd + i * 36 + 27);

        auto to_yz_yz = pbuffer.data(idx_op_dd + i * 36 + 28);

        auto to_yz_zz = pbuffer.data(idx_op_dd + i * 36 + 29);

        auto to_zz_xx = pbuffer.data(idx_op_dd + i * 36 + 30);

        auto to_zz_xy = pbuffer.data(idx_op_dd + i * 36 + 31);

        auto to_zz_xz = pbuffer.data(idx_op_dd + i * 36 + 32);

        auto to_zz_yy = pbuffer.data(idx_op_dd + i * 36 + 33);

        auto to_zz_yz = pbuffer.data(idx_op_dd + i * 36 + 34);

        auto to_zz_zz = pbuffer.data(idx_op_dd + i * 36 + 35);

        // Set up components of auxiliary buffer : GS

        auto to_xxxx_0 = pbuffer.data(idx_op_gs + i * 15 + 0);

        auto to_xxxy_0 = pbuffer.data(idx_op_gs + i * 15 + 1);

        auto to_xxxz_0 = pbuffer.data(idx_op_gs + i * 15 + 2);

        auto to_xxyy_0 = pbuffer.data(idx_op_gs + i * 15 + 3);

        auto to_xxyz_0 = pbuffer.data(idx_op_gs + i * 15 + 4);

        auto to_xxzz_0 = pbuffer.data(idx_op_gs + i * 15 + 5);

        auto to_xyyy_0 = pbuffer.data(idx_op_gs + i * 15 + 6);

        auto to_xyyz_0 = pbuffer.data(idx_op_gs + i * 15 + 7);

        auto to_xyzz_0 = pbuffer.data(idx_op_gs + i * 15 + 8);

        auto to_xzzz_0 = pbuffer.data(idx_op_gs + i * 15 + 9);

        auto to_yyyy_0 = pbuffer.data(idx_op_gs + i * 15 + 10);

        auto to_yyyz_0 = pbuffer.data(idx_op_gs + i * 15 + 11);

        auto to_yyzz_0 = pbuffer.data(idx_op_gs + i * 15 + 12);

        auto to_yzzz_0 = pbuffer.data(idx_op_gs + i * 15 + 13);

        auto to_zzzz_0 = pbuffer.data(idx_op_gs + i * 15 + 14);

        // Set up components of auxiliary buffer : GD

        auto to_xxxx_xx = pbuffer.data(idx_op_gd + i * 90 + 0);

        auto to_xxxx_xy = pbuffer.data(idx_op_gd + i * 90 + 1);

        auto to_xxxx_xz = pbuffer.data(idx_op_gd + i * 90 + 2);

        auto to_xxxx_yy = pbuffer.data(idx_op_gd + i * 90 + 3);

        auto to_xxxx_yz = pbuffer.data(idx_op_gd + i * 90 + 4);

        auto to_xxxx_zz = pbuffer.data(idx_op_gd + i * 90 + 5);

        auto to_xxxy_xx = pbuffer.data(idx_op_gd + i * 90 + 6);

        auto to_xxxy_xy = pbuffer.data(idx_op_gd + i * 90 + 7);

        auto to_xxxy_xz = pbuffer.data(idx_op_gd + i * 90 + 8);

        auto to_xxxy_yy = pbuffer.data(idx_op_gd + i * 90 + 9);

        auto to_xxxy_yz = pbuffer.data(idx_op_gd + i * 90 + 10);

        auto to_xxxy_zz = pbuffer.data(idx_op_gd + i * 90 + 11);

        auto to_xxxz_xx = pbuffer.data(idx_op_gd + i * 90 + 12);

        auto to_xxxz_xy = pbuffer.data(idx_op_gd + i * 90 + 13);

        auto to_xxxz_xz = pbuffer.data(idx_op_gd + i * 90 + 14);

        auto to_xxxz_yy = pbuffer.data(idx_op_gd + i * 90 + 15);

        auto to_xxxz_yz = pbuffer.data(idx_op_gd + i * 90 + 16);

        auto to_xxxz_zz = pbuffer.data(idx_op_gd + i * 90 + 17);

        auto to_xxyy_xx = pbuffer.data(idx_op_gd + i * 90 + 18);

        auto to_xxyy_xy = pbuffer.data(idx_op_gd + i * 90 + 19);

        auto to_xxyy_xz = pbuffer.data(idx_op_gd + i * 90 + 20);

        auto to_xxyy_yy = pbuffer.data(idx_op_gd + i * 90 + 21);

        auto to_xxyy_yz = pbuffer.data(idx_op_gd + i * 90 + 22);

        auto to_xxyy_zz = pbuffer.data(idx_op_gd + i * 90 + 23);

        auto to_xxyz_xx = pbuffer.data(idx_op_gd + i * 90 + 24);

        auto to_xxyz_xy = pbuffer.data(idx_op_gd + i * 90 + 25);

        auto to_xxyz_xz = pbuffer.data(idx_op_gd + i * 90 + 26);

        auto to_xxyz_yy = pbuffer.data(idx_op_gd + i * 90 + 27);

        auto to_xxyz_yz = pbuffer.data(idx_op_gd + i * 90 + 28);

        auto to_xxyz_zz = pbuffer.data(idx_op_gd + i * 90 + 29);

        auto to_xxzz_xx = pbuffer.data(idx_op_gd + i * 90 + 30);

        auto to_xxzz_xy = pbuffer.data(idx_op_gd + i * 90 + 31);

        auto to_xxzz_xz = pbuffer.data(idx_op_gd + i * 90 + 32);

        auto to_xxzz_yy = pbuffer.data(idx_op_gd + i * 90 + 33);

        auto to_xxzz_yz = pbuffer.data(idx_op_gd + i * 90 + 34);

        auto to_xxzz_zz = pbuffer.data(idx_op_gd + i * 90 + 35);

        auto to_xyyy_xx = pbuffer.data(idx_op_gd + i * 90 + 36);

        auto to_xyyy_xy = pbuffer.data(idx_op_gd + i * 90 + 37);

        auto to_xyyy_xz = pbuffer.data(idx_op_gd + i * 90 + 38);

        auto to_xyyy_yy = pbuffer.data(idx_op_gd + i * 90 + 39);

        auto to_xyyy_yz = pbuffer.data(idx_op_gd + i * 90 + 40);

        auto to_xyyy_zz = pbuffer.data(idx_op_gd + i * 90 + 41);

        auto to_xyyz_xx = pbuffer.data(idx_op_gd + i * 90 + 42);

        auto to_xyyz_xy = pbuffer.data(idx_op_gd + i * 90 + 43);

        auto to_xyyz_xz = pbuffer.data(idx_op_gd + i * 90 + 44);

        auto to_xyyz_yy = pbuffer.data(idx_op_gd + i * 90 + 45);

        auto to_xyyz_yz = pbuffer.data(idx_op_gd + i * 90 + 46);

        auto to_xyyz_zz = pbuffer.data(idx_op_gd + i * 90 + 47);

        auto to_xyzz_xx = pbuffer.data(idx_op_gd + i * 90 + 48);

        auto to_xyzz_xy = pbuffer.data(idx_op_gd + i * 90 + 49);

        auto to_xyzz_xz = pbuffer.data(idx_op_gd + i * 90 + 50);

        auto to_xyzz_yy = pbuffer.data(idx_op_gd + i * 90 + 51);

        auto to_xyzz_yz = pbuffer.data(idx_op_gd + i * 90 + 52);

        auto to_xyzz_zz = pbuffer.data(idx_op_gd + i * 90 + 53);

        auto to_xzzz_xx = pbuffer.data(idx_op_gd + i * 90 + 54);

        auto to_xzzz_xy = pbuffer.data(idx_op_gd + i * 90 + 55);

        auto to_xzzz_xz = pbuffer.data(idx_op_gd + i * 90 + 56);

        auto to_xzzz_yy = pbuffer.data(idx_op_gd + i * 90 + 57);

        auto to_xzzz_yz = pbuffer.data(idx_op_gd + i * 90 + 58);

        auto to_xzzz_zz = pbuffer.data(idx_op_gd + i * 90 + 59);

        auto to_yyyy_xx = pbuffer.data(idx_op_gd + i * 90 + 60);

        auto to_yyyy_xy = pbuffer.data(idx_op_gd + i * 90 + 61);

        auto to_yyyy_xz = pbuffer.data(idx_op_gd + i * 90 + 62);

        auto to_yyyy_yy = pbuffer.data(idx_op_gd + i * 90 + 63);

        auto to_yyyy_yz = pbuffer.data(idx_op_gd + i * 90 + 64);

        auto to_yyyy_zz = pbuffer.data(idx_op_gd + i * 90 + 65);

        auto to_yyyz_xx = pbuffer.data(idx_op_gd + i * 90 + 66);

        auto to_yyyz_xy = pbuffer.data(idx_op_gd + i * 90 + 67);

        auto to_yyyz_xz = pbuffer.data(idx_op_gd + i * 90 + 68);

        auto to_yyyz_yy = pbuffer.data(idx_op_gd + i * 90 + 69);

        auto to_yyyz_yz = pbuffer.data(idx_op_gd + i * 90 + 70);

        auto to_yyyz_zz = pbuffer.data(idx_op_gd + i * 90 + 71);

        auto to_yyzz_xx = pbuffer.data(idx_op_gd + i * 90 + 72);

        auto to_yyzz_xy = pbuffer.data(idx_op_gd + i * 90 + 73);

        auto to_yyzz_xz = pbuffer.data(idx_op_gd + i * 90 + 74);

        auto to_yyzz_yy = pbuffer.data(idx_op_gd + i * 90 + 75);

        auto to_yyzz_yz = pbuffer.data(idx_op_gd + i * 90 + 76);

        auto to_yyzz_zz = pbuffer.data(idx_op_gd + i * 90 + 77);

        auto to_yzzz_xx = pbuffer.data(idx_op_gd + i * 90 + 78);

        auto to_yzzz_xy = pbuffer.data(idx_op_gd + i * 90 + 79);

        auto to_yzzz_xz = pbuffer.data(idx_op_gd + i * 90 + 80);

        auto to_yzzz_yy = pbuffer.data(idx_op_gd + i * 90 + 81);

        auto to_yzzz_yz = pbuffer.data(idx_op_gd + i * 90 + 82);

        auto to_yzzz_zz = pbuffer.data(idx_op_gd + i * 90 + 83);

        auto to_zzzz_xx = pbuffer.data(idx_op_gd + i * 90 + 84);

        auto to_zzzz_xy = pbuffer.data(idx_op_gd + i * 90 + 85);

        auto to_zzzz_xz = pbuffer.data(idx_op_gd + i * 90 + 86);

        auto to_zzzz_yy = pbuffer.data(idx_op_gd + i * 90 + 87);

        auto to_zzzz_yz = pbuffer.data(idx_op_gd + i * 90 + 88);

        auto to_zzzz_zz = pbuffer.data(idx_op_gd + i * 90 + 89);

        // Set up 0-3 components of targeted buffer : FP

        auto to_x_x_xxx_x = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 0);

        auto to_x_x_xxx_y = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 1);

        auto to_x_x_xxx_z = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_x_x_xxx_x, to_x_x_xxx_y, to_x_x_xxx_z, to_xx_0, to_xx_xx, to_xx_xy, to_xx_xz, to_xxxx_0, to_xxxx_xx, to_xxxx_xy, to_xxxx_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxx_x[k] = 3.0 * to_xx_0[k] - 6.0 * to_xx_xx[k] * tke_0 - 2.0 * to_xxxx_0[k] * tbe_0 + 4.0 * to_xxxx_xx[k] * tbe_0 * tke_0;

            to_x_x_xxx_y[k] = -6.0 * to_xx_xy[k] * tke_0 + 4.0 * to_xxxx_xy[k] * tbe_0 * tke_0;

            to_x_x_xxx_z[k] = -6.0 * to_xx_xz[k] * tke_0 + 4.0 * to_xxxx_xz[k] * tbe_0 * tke_0;
        }

        // Set up 3-6 components of targeted buffer : FP

        auto to_x_x_xxy_x = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 3);

        auto to_x_x_xxy_y = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 4);

        auto to_x_x_xxy_z = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_x_x_xxy_x, to_x_x_xxy_y, to_x_x_xxy_z, to_xxxy_0, to_xxxy_xx, to_xxxy_xy, to_xxxy_xz, to_xy_0, to_xy_xx, to_xy_xy, to_xy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxy_x[k] = 2.0 * to_xy_0[k] - 4.0 * to_xy_xx[k] * tke_0 - 2.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxy_xx[k] * tbe_0 * tke_0;

            to_x_x_xxy_y[k] = -4.0 * to_xy_xy[k] * tke_0 + 4.0 * to_xxxy_xy[k] * tbe_0 * tke_0;

            to_x_x_xxy_z[k] = -4.0 * to_xy_xz[k] * tke_0 + 4.0 * to_xxxy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 6-9 components of targeted buffer : FP

        auto to_x_x_xxz_x = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 6);

        auto to_x_x_xxz_y = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 7);

        auto to_x_x_xxz_z = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_x_x_xxz_x, to_x_x_xxz_y, to_x_x_xxz_z, to_xxxz_0, to_xxxz_xx, to_xxxz_xy, to_xxxz_xz, to_xz_0, to_xz_xx, to_xz_xy, to_xz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxz_x[k] = 2.0 * to_xz_0[k] - 4.0 * to_xz_xx[k] * tke_0 - 2.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxz_xx[k] * tbe_0 * tke_0;

            to_x_x_xxz_y[k] = -4.0 * to_xz_xy[k] * tke_0 + 4.0 * to_xxxz_xy[k] * tbe_0 * tke_0;

            to_x_x_xxz_z[k] = -4.0 * to_xz_xz[k] * tke_0 + 4.0 * to_xxxz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 9-12 components of targeted buffer : FP

        auto to_x_x_xyy_x = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 9);

        auto to_x_x_xyy_y = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 10);

        auto to_x_x_xyy_z = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_x_x_xyy_x, to_x_x_xyy_y, to_x_x_xyy_z, to_xxyy_0, to_xxyy_xx, to_xxyy_xy, to_xxyy_xz, to_yy_0, to_yy_xx, to_yy_xy, to_yy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyy_x[k] = to_yy_0[k] - 2.0 * to_yy_xx[k] * tke_0 - 2.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxyy_xx[k] * tbe_0 * tke_0;

            to_x_x_xyy_y[k] = -2.0 * to_yy_xy[k] * tke_0 + 4.0 * to_xxyy_xy[k] * tbe_0 * tke_0;

            to_x_x_xyy_z[k] = -2.0 * to_yy_xz[k] * tke_0 + 4.0 * to_xxyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 12-15 components of targeted buffer : FP

        auto to_x_x_xyz_x = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 12);

        auto to_x_x_xyz_y = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 13);

        auto to_x_x_xyz_z = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_x_x_xyz_x, to_x_x_xyz_y, to_x_x_xyz_z, to_xxyz_0, to_xxyz_xx, to_xxyz_xy, to_xxyz_xz, to_yz_0, to_yz_xx, to_yz_xy, to_yz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyz_x[k] = to_yz_0[k] - 2.0 * to_yz_xx[k] * tke_0 - 2.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyz_xx[k] * tbe_0 * tke_0;

            to_x_x_xyz_y[k] = -2.0 * to_yz_xy[k] * tke_0 + 4.0 * to_xxyz_xy[k] * tbe_0 * tke_0;

            to_x_x_xyz_z[k] = -2.0 * to_yz_xz[k] * tke_0 + 4.0 * to_xxyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 15-18 components of targeted buffer : FP

        auto to_x_x_xzz_x = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 15);

        auto to_x_x_xzz_y = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 16);

        auto to_x_x_xzz_z = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_x_x_xzz_x, to_x_x_xzz_y, to_x_x_xzz_z, to_xxzz_0, to_xxzz_xx, to_xxzz_xy, to_xxzz_xz, to_zz_0, to_zz_xx, to_zz_xy, to_zz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xzz_x[k] = to_zz_0[k] - 2.0 * to_zz_xx[k] * tke_0 - 2.0 * to_xxzz_0[k] * tbe_0 + 4.0 * to_xxzz_xx[k] * tbe_0 * tke_0;

            to_x_x_xzz_y[k] = -2.0 * to_zz_xy[k] * tke_0 + 4.0 * to_xxzz_xy[k] * tbe_0 * tke_0;

            to_x_x_xzz_z[k] = -2.0 * to_zz_xz[k] * tke_0 + 4.0 * to_xxzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 18-21 components of targeted buffer : FP

        auto to_x_x_yyy_x = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 18);

        auto to_x_x_yyy_y = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 19);

        auto to_x_x_yyy_z = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_x_x_yyy_x, to_x_x_yyy_y, to_x_x_yyy_z, to_xyyy_0, to_xyyy_xx, to_xyyy_xy, to_xyyy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyy_x[k] = -2.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xyyy_xx[k] * tbe_0 * tke_0;

            to_x_x_yyy_y[k] = 4.0 * to_xyyy_xy[k] * tbe_0 * tke_0;

            to_x_x_yyy_z[k] = 4.0 * to_xyyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 21-24 components of targeted buffer : FP

        auto to_x_x_yyz_x = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 21);

        auto to_x_x_yyz_y = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 22);

        auto to_x_x_yyz_z = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_x_x_yyz_x, to_x_x_yyz_y, to_x_x_yyz_z, to_xyyz_0, to_xyyz_xx, to_xyyz_xy, to_xyyz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyz_x[k] = -2.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyz_xx[k] * tbe_0 * tke_0;

            to_x_x_yyz_y[k] = 4.0 * to_xyyz_xy[k] * tbe_0 * tke_0;

            to_x_x_yyz_z[k] = 4.0 * to_xyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 24-27 components of targeted buffer : FP

        auto to_x_x_yzz_x = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 24);

        auto to_x_x_yzz_y = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 25);

        auto to_x_x_yzz_z = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_x_x_yzz_x, to_x_x_yzz_y, to_x_x_yzz_z, to_xyzz_0, to_xyzz_xx, to_xyzz_xy, to_xyzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yzz_x[k] = -2.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzz_xx[k] * tbe_0 * tke_0;

            to_x_x_yzz_y[k] = 4.0 * to_xyzz_xy[k] * tbe_0 * tke_0;

            to_x_x_yzz_z[k] = 4.0 * to_xyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 27-30 components of targeted buffer : FP

        auto to_x_x_zzz_x = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 27);

        auto to_x_x_zzz_y = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 28);

        auto to_x_x_zzz_z = pbuffer.data(idx_op_geom_101_fp + 0 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_x_x_zzz_x, to_x_x_zzz_y, to_x_x_zzz_z, to_xzzz_0, to_xzzz_xx, to_xzzz_xy, to_xzzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zzz_x[k] = -2.0 * to_xzzz_0[k] * tbe_0 + 4.0 * to_xzzz_xx[k] * tbe_0 * tke_0;

            to_x_x_zzz_y[k] = 4.0 * to_xzzz_xy[k] * tbe_0 * tke_0;

            to_x_x_zzz_z[k] = 4.0 * to_xzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 30-33 components of targeted buffer : FP

        auto to_x_y_xxx_x = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 0);

        auto to_x_y_xxx_y = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 1);

        auto to_x_y_xxx_z = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_x_y_xxx_x, to_x_y_xxx_y, to_x_y_xxx_z, to_xx_0, to_xx_xy, to_xx_yy, to_xx_yz, to_xxxx_0, to_xxxx_xy, to_xxxx_yy, to_xxxx_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxx_x[k] = -6.0 * to_xx_xy[k] * tke_0 + 4.0 * to_xxxx_xy[k] * tbe_0 * tke_0;

            to_x_y_xxx_y[k] = 3.0 * to_xx_0[k] - 6.0 * to_xx_yy[k] * tke_0 - 2.0 * to_xxxx_0[k] * tbe_0 + 4.0 * to_xxxx_yy[k] * tbe_0 * tke_0;

            to_x_y_xxx_z[k] = -6.0 * to_xx_yz[k] * tke_0 + 4.0 * to_xxxx_yz[k] * tbe_0 * tke_0;
        }

        // Set up 33-36 components of targeted buffer : FP

        auto to_x_y_xxy_x = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 3);

        auto to_x_y_xxy_y = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 4);

        auto to_x_y_xxy_z = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_x_y_xxy_x, to_x_y_xxy_y, to_x_y_xxy_z, to_xxxy_0, to_xxxy_xy, to_xxxy_yy, to_xxxy_yz, to_xy_0, to_xy_xy, to_xy_yy, to_xy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxy_x[k] = -4.0 * to_xy_xy[k] * tke_0 + 4.0 * to_xxxy_xy[k] * tbe_0 * tke_0;

            to_x_y_xxy_y[k] = 2.0 * to_xy_0[k] - 4.0 * to_xy_yy[k] * tke_0 - 2.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxy_yy[k] * tbe_0 * tke_0;

            to_x_y_xxy_z[k] = -4.0 * to_xy_yz[k] * tke_0 + 4.0 * to_xxxy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 36-39 components of targeted buffer : FP

        auto to_x_y_xxz_x = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 6);

        auto to_x_y_xxz_y = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 7);

        auto to_x_y_xxz_z = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_x_y_xxz_x, to_x_y_xxz_y, to_x_y_xxz_z, to_xxxz_0, to_xxxz_xy, to_xxxz_yy, to_xxxz_yz, to_xz_0, to_xz_xy, to_xz_yy, to_xz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxz_x[k] = -4.0 * to_xz_xy[k] * tke_0 + 4.0 * to_xxxz_xy[k] * tbe_0 * tke_0;

            to_x_y_xxz_y[k] = 2.0 * to_xz_0[k] - 4.0 * to_xz_yy[k] * tke_0 - 2.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxz_yy[k] * tbe_0 * tke_0;

            to_x_y_xxz_z[k] = -4.0 * to_xz_yz[k] * tke_0 + 4.0 * to_xxxz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 39-42 components of targeted buffer : FP

        auto to_x_y_xyy_x = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 9);

        auto to_x_y_xyy_y = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 10);

        auto to_x_y_xyy_z = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_x_y_xyy_x, to_x_y_xyy_y, to_x_y_xyy_z, to_xxyy_0, to_xxyy_xy, to_xxyy_yy, to_xxyy_yz, to_yy_0, to_yy_xy, to_yy_yy, to_yy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyy_x[k] = -2.0 * to_yy_xy[k] * tke_0 + 4.0 * to_xxyy_xy[k] * tbe_0 * tke_0;

            to_x_y_xyy_y[k] = to_yy_0[k] - 2.0 * to_yy_yy[k] * tke_0 - 2.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxyy_yy[k] * tbe_0 * tke_0;

            to_x_y_xyy_z[k] = -2.0 * to_yy_yz[k] * tke_0 + 4.0 * to_xxyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 42-45 components of targeted buffer : FP

        auto to_x_y_xyz_x = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 12);

        auto to_x_y_xyz_y = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 13);

        auto to_x_y_xyz_z = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_x_y_xyz_x, to_x_y_xyz_y, to_x_y_xyz_z, to_xxyz_0, to_xxyz_xy, to_xxyz_yy, to_xxyz_yz, to_yz_0, to_yz_xy, to_yz_yy, to_yz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyz_x[k] = -2.0 * to_yz_xy[k] * tke_0 + 4.0 * to_xxyz_xy[k] * tbe_0 * tke_0;

            to_x_y_xyz_y[k] = to_yz_0[k] - 2.0 * to_yz_yy[k] * tke_0 - 2.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyz_yy[k] * tbe_0 * tke_0;

            to_x_y_xyz_z[k] = -2.0 * to_yz_yz[k] * tke_0 + 4.0 * to_xxyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 45-48 components of targeted buffer : FP

        auto to_x_y_xzz_x = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 15);

        auto to_x_y_xzz_y = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 16);

        auto to_x_y_xzz_z = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_x_y_xzz_x, to_x_y_xzz_y, to_x_y_xzz_z, to_xxzz_0, to_xxzz_xy, to_xxzz_yy, to_xxzz_yz, to_zz_0, to_zz_xy, to_zz_yy, to_zz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xzz_x[k] = -2.0 * to_zz_xy[k] * tke_0 + 4.0 * to_xxzz_xy[k] * tbe_0 * tke_0;

            to_x_y_xzz_y[k] = to_zz_0[k] - 2.0 * to_zz_yy[k] * tke_0 - 2.0 * to_xxzz_0[k] * tbe_0 + 4.0 * to_xxzz_yy[k] * tbe_0 * tke_0;

            to_x_y_xzz_z[k] = -2.0 * to_zz_yz[k] * tke_0 + 4.0 * to_xxzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 48-51 components of targeted buffer : FP

        auto to_x_y_yyy_x = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 18);

        auto to_x_y_yyy_y = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 19);

        auto to_x_y_yyy_z = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_x_y_yyy_x, to_x_y_yyy_y, to_x_y_yyy_z, to_xyyy_0, to_xyyy_xy, to_xyyy_yy, to_xyyy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyy_x[k] = 4.0 * to_xyyy_xy[k] * tbe_0 * tke_0;

            to_x_y_yyy_y[k] = -2.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xyyy_yy[k] * tbe_0 * tke_0;

            to_x_y_yyy_z[k] = 4.0 * to_xyyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 51-54 components of targeted buffer : FP

        auto to_x_y_yyz_x = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 21);

        auto to_x_y_yyz_y = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 22);

        auto to_x_y_yyz_z = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_x_y_yyz_x, to_x_y_yyz_y, to_x_y_yyz_z, to_xyyz_0, to_xyyz_xy, to_xyyz_yy, to_xyyz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyz_x[k] = 4.0 * to_xyyz_xy[k] * tbe_0 * tke_0;

            to_x_y_yyz_y[k] = -2.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyz_yy[k] * tbe_0 * tke_0;

            to_x_y_yyz_z[k] = 4.0 * to_xyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 54-57 components of targeted buffer : FP

        auto to_x_y_yzz_x = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 24);

        auto to_x_y_yzz_y = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 25);

        auto to_x_y_yzz_z = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_x_y_yzz_x, to_x_y_yzz_y, to_x_y_yzz_z, to_xyzz_0, to_xyzz_xy, to_xyzz_yy, to_xyzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yzz_x[k] = 4.0 * to_xyzz_xy[k] * tbe_0 * tke_0;

            to_x_y_yzz_y[k] = -2.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzz_yy[k] * tbe_0 * tke_0;

            to_x_y_yzz_z[k] = 4.0 * to_xyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 57-60 components of targeted buffer : FP

        auto to_x_y_zzz_x = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 27);

        auto to_x_y_zzz_y = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 28);

        auto to_x_y_zzz_z = pbuffer.data(idx_op_geom_101_fp + 1 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_x_y_zzz_x, to_x_y_zzz_y, to_x_y_zzz_z, to_xzzz_0, to_xzzz_xy, to_xzzz_yy, to_xzzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zzz_x[k] = 4.0 * to_xzzz_xy[k] * tbe_0 * tke_0;

            to_x_y_zzz_y[k] = -2.0 * to_xzzz_0[k] * tbe_0 + 4.0 * to_xzzz_yy[k] * tbe_0 * tke_0;

            to_x_y_zzz_z[k] = 4.0 * to_xzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 60-63 components of targeted buffer : FP

        auto to_x_z_xxx_x = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 0);

        auto to_x_z_xxx_y = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 1);

        auto to_x_z_xxx_z = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_x_z_xxx_x, to_x_z_xxx_y, to_x_z_xxx_z, to_xx_0, to_xx_xz, to_xx_yz, to_xx_zz, to_xxxx_0, to_xxxx_xz, to_xxxx_yz, to_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxx_x[k] = -6.0 * to_xx_xz[k] * tke_0 + 4.0 * to_xxxx_xz[k] * tbe_0 * tke_0;

            to_x_z_xxx_y[k] = -6.0 * to_xx_yz[k] * tke_0 + 4.0 * to_xxxx_yz[k] * tbe_0 * tke_0;

            to_x_z_xxx_z[k] = 3.0 * to_xx_0[k] - 6.0 * to_xx_zz[k] * tke_0 - 2.0 * to_xxxx_0[k] * tbe_0 + 4.0 * to_xxxx_zz[k] * tbe_0 * tke_0;
        }

        // Set up 63-66 components of targeted buffer : FP

        auto to_x_z_xxy_x = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 3);

        auto to_x_z_xxy_y = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 4);

        auto to_x_z_xxy_z = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_x_z_xxy_x, to_x_z_xxy_y, to_x_z_xxy_z, to_xxxy_0, to_xxxy_xz, to_xxxy_yz, to_xxxy_zz, to_xy_0, to_xy_xz, to_xy_yz, to_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxy_x[k] = -4.0 * to_xy_xz[k] * tke_0 + 4.0 * to_xxxy_xz[k] * tbe_0 * tke_0;

            to_x_z_xxy_y[k] = -4.0 * to_xy_yz[k] * tke_0 + 4.0 * to_xxxy_yz[k] * tbe_0 * tke_0;

            to_x_z_xxy_z[k] = 2.0 * to_xy_0[k] - 4.0 * to_xy_zz[k] * tke_0 - 2.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 66-69 components of targeted buffer : FP

        auto to_x_z_xxz_x = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 6);

        auto to_x_z_xxz_y = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 7);

        auto to_x_z_xxz_z = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_x_z_xxz_x, to_x_z_xxz_y, to_x_z_xxz_z, to_xxxz_0, to_xxxz_xz, to_xxxz_yz, to_xxxz_zz, to_xz_0, to_xz_xz, to_xz_yz, to_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxz_x[k] = -4.0 * to_xz_xz[k] * tke_0 + 4.0 * to_xxxz_xz[k] * tbe_0 * tke_0;

            to_x_z_xxz_y[k] = -4.0 * to_xz_yz[k] * tke_0 + 4.0 * to_xxxz_yz[k] * tbe_0 * tke_0;

            to_x_z_xxz_z[k] = 2.0 * to_xz_0[k] - 4.0 * to_xz_zz[k] * tke_0 - 2.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 69-72 components of targeted buffer : FP

        auto to_x_z_xyy_x = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 9);

        auto to_x_z_xyy_y = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 10);

        auto to_x_z_xyy_z = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_x_z_xyy_x, to_x_z_xyy_y, to_x_z_xyy_z, to_xxyy_0, to_xxyy_xz, to_xxyy_yz, to_xxyy_zz, to_yy_0, to_yy_xz, to_yy_yz, to_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyy_x[k] = -2.0 * to_yy_xz[k] * tke_0 + 4.0 * to_xxyy_xz[k] * tbe_0 * tke_0;

            to_x_z_xyy_y[k] = -2.0 * to_yy_yz[k] * tke_0 + 4.0 * to_xxyy_yz[k] * tbe_0 * tke_0;

            to_x_z_xyy_z[k] = to_yy_0[k] - 2.0 * to_yy_zz[k] * tke_0 - 2.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 72-75 components of targeted buffer : FP

        auto to_x_z_xyz_x = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 12);

        auto to_x_z_xyz_y = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 13);

        auto to_x_z_xyz_z = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_x_z_xyz_x, to_x_z_xyz_y, to_x_z_xyz_z, to_xxyz_0, to_xxyz_xz, to_xxyz_yz, to_xxyz_zz, to_yz_0, to_yz_xz, to_yz_yz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyz_x[k] = -2.0 * to_yz_xz[k] * tke_0 + 4.0 * to_xxyz_xz[k] * tbe_0 * tke_0;

            to_x_z_xyz_y[k] = -2.0 * to_yz_yz[k] * tke_0 + 4.0 * to_xxyz_yz[k] * tbe_0 * tke_0;

            to_x_z_xyz_z[k] = to_yz_0[k] - 2.0 * to_yz_zz[k] * tke_0 - 2.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 75-78 components of targeted buffer : FP

        auto to_x_z_xzz_x = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 15);

        auto to_x_z_xzz_y = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 16);

        auto to_x_z_xzz_z = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_x_z_xzz_x, to_x_z_xzz_y, to_x_z_xzz_z, to_xxzz_0, to_xxzz_xz, to_xxzz_yz, to_xxzz_zz, to_zz_0, to_zz_xz, to_zz_yz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xzz_x[k] = -2.0 * to_zz_xz[k] * tke_0 + 4.0 * to_xxzz_xz[k] * tbe_0 * tke_0;

            to_x_z_xzz_y[k] = -2.0 * to_zz_yz[k] * tke_0 + 4.0 * to_xxzz_yz[k] * tbe_0 * tke_0;

            to_x_z_xzz_z[k] = to_zz_0[k] - 2.0 * to_zz_zz[k] * tke_0 - 2.0 * to_xxzz_0[k] * tbe_0 + 4.0 * to_xxzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 78-81 components of targeted buffer : FP

        auto to_x_z_yyy_x = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 18);

        auto to_x_z_yyy_y = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 19);

        auto to_x_z_yyy_z = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_x_z_yyy_x, to_x_z_yyy_y, to_x_z_yyy_z, to_xyyy_0, to_xyyy_xz, to_xyyy_yz, to_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyy_x[k] = 4.0 * to_xyyy_xz[k] * tbe_0 * tke_0;

            to_x_z_yyy_y[k] = 4.0 * to_xyyy_yz[k] * tbe_0 * tke_0;

            to_x_z_yyy_z[k] = -2.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xyyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 81-84 components of targeted buffer : FP

        auto to_x_z_yyz_x = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 21);

        auto to_x_z_yyz_y = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 22);

        auto to_x_z_yyz_z = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_x_z_yyz_x, to_x_z_yyz_y, to_x_z_yyz_z, to_xyyz_0, to_xyyz_xz, to_xyyz_yz, to_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyz_x[k] = 4.0 * to_xyyz_xz[k] * tbe_0 * tke_0;

            to_x_z_yyz_y[k] = 4.0 * to_xyyz_yz[k] * tbe_0 * tke_0;

            to_x_z_yyz_z[k] = -2.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 84-87 components of targeted buffer : FP

        auto to_x_z_yzz_x = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 24);

        auto to_x_z_yzz_y = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 25);

        auto to_x_z_yzz_z = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_x_z_yzz_x, to_x_z_yzz_y, to_x_z_yzz_z, to_xyzz_0, to_xyzz_xz, to_xyzz_yz, to_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yzz_x[k] = 4.0 * to_xyzz_xz[k] * tbe_0 * tke_0;

            to_x_z_yzz_y[k] = 4.0 * to_xyzz_yz[k] * tbe_0 * tke_0;

            to_x_z_yzz_z[k] = -2.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 87-90 components of targeted buffer : FP

        auto to_x_z_zzz_x = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 27);

        auto to_x_z_zzz_y = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 28);

        auto to_x_z_zzz_z = pbuffer.data(idx_op_geom_101_fp + 2 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_x_z_zzz_x, to_x_z_zzz_y, to_x_z_zzz_z, to_xzzz_0, to_xzzz_xz, to_xzzz_yz, to_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zzz_x[k] = 4.0 * to_xzzz_xz[k] * tbe_0 * tke_0;

            to_x_z_zzz_y[k] = 4.0 * to_xzzz_yz[k] * tbe_0 * tke_0;

            to_x_z_zzz_z[k] = -2.0 * to_xzzz_0[k] * tbe_0 + 4.0 * to_xzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 90-93 components of targeted buffer : FP

        auto to_y_x_xxx_x = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 0);

        auto to_y_x_xxx_y = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 1);

        auto to_y_x_xxx_z = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_xxxy_0, to_xxxy_xx, to_xxxy_xy, to_xxxy_xz, to_y_x_xxx_x, to_y_x_xxx_y, to_y_x_xxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxx_x[k] = -2.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxy_xx[k] * tbe_0 * tke_0;

            to_y_x_xxx_y[k] = 4.0 * to_xxxy_xy[k] * tbe_0 * tke_0;

            to_y_x_xxx_z[k] = 4.0 * to_xxxy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 93-96 components of targeted buffer : FP

        auto to_y_x_xxy_x = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 3);

        auto to_y_x_xxy_y = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 4);

        auto to_y_x_xxy_z = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_xx_0, to_xx_xx, to_xx_xy, to_xx_xz, to_xxyy_0, to_xxyy_xx, to_xxyy_xy, to_xxyy_xz, to_y_x_xxy_x, to_y_x_xxy_y, to_y_x_xxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxy_x[k] = to_xx_0[k] - 2.0 * to_xx_xx[k] * tke_0 - 2.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxyy_xx[k] * tbe_0 * tke_0;

            to_y_x_xxy_y[k] = -2.0 * to_xx_xy[k] * tke_0 + 4.0 * to_xxyy_xy[k] * tbe_0 * tke_0;

            to_y_x_xxy_z[k] = -2.0 * to_xx_xz[k] * tke_0 + 4.0 * to_xxyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 96-99 components of targeted buffer : FP

        auto to_y_x_xxz_x = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 6);

        auto to_y_x_xxz_y = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 7);

        auto to_y_x_xxz_z = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_xxyz_0, to_xxyz_xx, to_xxyz_xy, to_xxyz_xz, to_y_x_xxz_x, to_y_x_xxz_y, to_y_x_xxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxz_x[k] = -2.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyz_xx[k] * tbe_0 * tke_0;

            to_y_x_xxz_y[k] = 4.0 * to_xxyz_xy[k] * tbe_0 * tke_0;

            to_y_x_xxz_z[k] = 4.0 * to_xxyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 99-102 components of targeted buffer : FP

        auto to_y_x_xyy_x = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 9);

        auto to_y_x_xyy_y = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 10);

        auto to_y_x_xyy_z = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_xy_0, to_xy_xx, to_xy_xy, to_xy_xz, to_xyyy_0, to_xyyy_xx, to_xyyy_xy, to_xyyy_xz, to_y_x_xyy_x, to_y_x_xyy_y, to_y_x_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyy_x[k] = 2.0 * to_xy_0[k] - 4.0 * to_xy_xx[k] * tke_0 - 2.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xyyy_xx[k] * tbe_0 * tke_0;

            to_y_x_xyy_y[k] = -4.0 * to_xy_xy[k] * tke_0 + 4.0 * to_xyyy_xy[k] * tbe_0 * tke_0;

            to_y_x_xyy_z[k] = -4.0 * to_xy_xz[k] * tke_0 + 4.0 * to_xyyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 102-105 components of targeted buffer : FP

        auto to_y_x_xyz_x = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 12);

        auto to_y_x_xyz_y = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 13);

        auto to_y_x_xyz_z = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_xyyz_0, to_xyyz_xx, to_xyyz_xy, to_xyyz_xz, to_xz_0, to_xz_xx, to_xz_xy, to_xz_xz, to_y_x_xyz_x, to_y_x_xyz_y, to_y_x_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyz_x[k] = to_xz_0[k] - 2.0 * to_xz_xx[k] * tke_0 - 2.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyz_xx[k] * tbe_0 * tke_0;

            to_y_x_xyz_y[k] = -2.0 * to_xz_xy[k] * tke_0 + 4.0 * to_xyyz_xy[k] * tbe_0 * tke_0;

            to_y_x_xyz_z[k] = -2.0 * to_xz_xz[k] * tke_0 + 4.0 * to_xyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 105-108 components of targeted buffer : FP

        auto to_y_x_xzz_x = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 15);

        auto to_y_x_xzz_y = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 16);

        auto to_y_x_xzz_z = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_xyzz_0, to_xyzz_xx, to_xyzz_xy, to_xyzz_xz, to_y_x_xzz_x, to_y_x_xzz_y, to_y_x_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xzz_x[k] = -2.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzz_xx[k] * tbe_0 * tke_0;

            to_y_x_xzz_y[k] = 4.0 * to_xyzz_xy[k] * tbe_0 * tke_0;

            to_y_x_xzz_z[k] = 4.0 * to_xyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 108-111 components of targeted buffer : FP

        auto to_y_x_yyy_x = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 18);

        auto to_y_x_yyy_y = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 19);

        auto to_y_x_yyy_z = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_y_x_yyy_x, to_y_x_yyy_y, to_y_x_yyy_z, to_yy_0, to_yy_xx, to_yy_xy, to_yy_xz, to_yyyy_0, to_yyyy_xx, to_yyyy_xy, to_yyyy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyy_x[k] = 3.0 * to_yy_0[k] - 6.0 * to_yy_xx[k] * tke_0 - 2.0 * to_yyyy_0[k] * tbe_0 + 4.0 * to_yyyy_xx[k] * tbe_0 * tke_0;

            to_y_x_yyy_y[k] = -6.0 * to_yy_xy[k] * tke_0 + 4.0 * to_yyyy_xy[k] * tbe_0 * tke_0;

            to_y_x_yyy_z[k] = -6.0 * to_yy_xz[k] * tke_0 + 4.0 * to_yyyy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 111-114 components of targeted buffer : FP

        auto to_y_x_yyz_x = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 21);

        auto to_y_x_yyz_y = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 22);

        auto to_y_x_yyz_z = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_y_x_yyz_x, to_y_x_yyz_y, to_y_x_yyz_z, to_yyyz_0, to_yyyz_xx, to_yyyz_xy, to_yyyz_xz, to_yz_0, to_yz_xx, to_yz_xy, to_yz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyz_x[k] = 2.0 * to_yz_0[k] - 4.0 * to_yz_xx[k] * tke_0 - 2.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_yyyz_xx[k] * tbe_0 * tke_0;

            to_y_x_yyz_y[k] = -4.0 * to_yz_xy[k] * tke_0 + 4.0 * to_yyyz_xy[k] * tbe_0 * tke_0;

            to_y_x_yyz_z[k] = -4.0 * to_yz_xz[k] * tke_0 + 4.0 * to_yyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 114-117 components of targeted buffer : FP

        auto to_y_x_yzz_x = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 24);

        auto to_y_x_yzz_y = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 25);

        auto to_y_x_yzz_z = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_y_x_yzz_x, to_y_x_yzz_y, to_y_x_yzz_z, to_yyzz_0, to_yyzz_xx, to_yyzz_xy, to_yyzz_xz, to_zz_0, to_zz_xx, to_zz_xy, to_zz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yzz_x[k] = to_zz_0[k] - 2.0 * to_zz_xx[k] * tke_0 - 2.0 * to_yyzz_0[k] * tbe_0 + 4.0 * to_yyzz_xx[k] * tbe_0 * tke_0;

            to_y_x_yzz_y[k] = -2.0 * to_zz_xy[k] * tke_0 + 4.0 * to_yyzz_xy[k] * tbe_0 * tke_0;

            to_y_x_yzz_z[k] = -2.0 * to_zz_xz[k] * tke_0 + 4.0 * to_yyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 117-120 components of targeted buffer : FP

        auto to_y_x_zzz_x = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 27);

        auto to_y_x_zzz_y = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 28);

        auto to_y_x_zzz_z = pbuffer.data(idx_op_geom_101_fp + 3 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_y_x_zzz_x, to_y_x_zzz_y, to_y_x_zzz_z, to_yzzz_0, to_yzzz_xx, to_yzzz_xy, to_yzzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zzz_x[k] = -2.0 * to_yzzz_0[k] * tbe_0 + 4.0 * to_yzzz_xx[k] * tbe_0 * tke_0;

            to_y_x_zzz_y[k] = 4.0 * to_yzzz_xy[k] * tbe_0 * tke_0;

            to_y_x_zzz_z[k] = 4.0 * to_yzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 120-123 components of targeted buffer : FP

        auto to_y_y_xxx_x = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 0);

        auto to_y_y_xxx_y = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 1);

        auto to_y_y_xxx_z = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_xxxy_0, to_xxxy_xy, to_xxxy_yy, to_xxxy_yz, to_y_y_xxx_x, to_y_y_xxx_y, to_y_y_xxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxx_x[k] = 4.0 * to_xxxy_xy[k] * tbe_0 * tke_0;

            to_y_y_xxx_y[k] = -2.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxy_yy[k] * tbe_0 * tke_0;

            to_y_y_xxx_z[k] = 4.0 * to_xxxy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 123-126 components of targeted buffer : FP

        auto to_y_y_xxy_x = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 3);

        auto to_y_y_xxy_y = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 4);

        auto to_y_y_xxy_z = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_xx_0, to_xx_xy, to_xx_yy, to_xx_yz, to_xxyy_0, to_xxyy_xy, to_xxyy_yy, to_xxyy_yz, to_y_y_xxy_x, to_y_y_xxy_y, to_y_y_xxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxy_x[k] = -2.0 * to_xx_xy[k] * tke_0 + 4.0 * to_xxyy_xy[k] * tbe_0 * tke_0;

            to_y_y_xxy_y[k] = to_xx_0[k] - 2.0 * to_xx_yy[k] * tke_0 - 2.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxyy_yy[k] * tbe_0 * tke_0;

            to_y_y_xxy_z[k] = -2.0 * to_xx_yz[k] * tke_0 + 4.0 * to_xxyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 126-129 components of targeted buffer : FP

        auto to_y_y_xxz_x = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 6);

        auto to_y_y_xxz_y = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 7);

        auto to_y_y_xxz_z = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_xxyz_0, to_xxyz_xy, to_xxyz_yy, to_xxyz_yz, to_y_y_xxz_x, to_y_y_xxz_y, to_y_y_xxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxz_x[k] = 4.0 * to_xxyz_xy[k] * tbe_0 * tke_0;

            to_y_y_xxz_y[k] = -2.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyz_yy[k] * tbe_0 * tke_0;

            to_y_y_xxz_z[k] = 4.0 * to_xxyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 129-132 components of targeted buffer : FP

        auto to_y_y_xyy_x = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 9);

        auto to_y_y_xyy_y = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 10);

        auto to_y_y_xyy_z = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_xy_0, to_xy_xy, to_xy_yy, to_xy_yz, to_xyyy_0, to_xyyy_xy, to_xyyy_yy, to_xyyy_yz, to_y_y_xyy_x, to_y_y_xyy_y, to_y_y_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyy_x[k] = -4.0 * to_xy_xy[k] * tke_0 + 4.0 * to_xyyy_xy[k] * tbe_0 * tke_0;

            to_y_y_xyy_y[k] = 2.0 * to_xy_0[k] - 4.0 * to_xy_yy[k] * tke_0 - 2.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xyyy_yy[k] * tbe_0 * tke_0;

            to_y_y_xyy_z[k] = -4.0 * to_xy_yz[k] * tke_0 + 4.0 * to_xyyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 132-135 components of targeted buffer : FP

        auto to_y_y_xyz_x = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 12);

        auto to_y_y_xyz_y = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 13);

        auto to_y_y_xyz_z = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_xyyz_0, to_xyyz_xy, to_xyyz_yy, to_xyyz_yz, to_xz_0, to_xz_xy, to_xz_yy, to_xz_yz, to_y_y_xyz_x, to_y_y_xyz_y, to_y_y_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyz_x[k] = -2.0 * to_xz_xy[k] * tke_0 + 4.0 * to_xyyz_xy[k] * tbe_0 * tke_0;

            to_y_y_xyz_y[k] = to_xz_0[k] - 2.0 * to_xz_yy[k] * tke_0 - 2.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyz_yy[k] * tbe_0 * tke_0;

            to_y_y_xyz_z[k] = -2.0 * to_xz_yz[k] * tke_0 + 4.0 * to_xyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 135-138 components of targeted buffer : FP

        auto to_y_y_xzz_x = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 15);

        auto to_y_y_xzz_y = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 16);

        auto to_y_y_xzz_z = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_xyzz_0, to_xyzz_xy, to_xyzz_yy, to_xyzz_yz, to_y_y_xzz_x, to_y_y_xzz_y, to_y_y_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xzz_x[k] = 4.0 * to_xyzz_xy[k] * tbe_0 * tke_0;

            to_y_y_xzz_y[k] = -2.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzz_yy[k] * tbe_0 * tke_0;

            to_y_y_xzz_z[k] = 4.0 * to_xyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 138-141 components of targeted buffer : FP

        auto to_y_y_yyy_x = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 18);

        auto to_y_y_yyy_y = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 19);

        auto to_y_y_yyy_z = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_y_y_yyy_x, to_y_y_yyy_y, to_y_y_yyy_z, to_yy_0, to_yy_xy, to_yy_yy, to_yy_yz, to_yyyy_0, to_yyyy_xy, to_yyyy_yy, to_yyyy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyy_x[k] = -6.0 * to_yy_xy[k] * tke_0 + 4.0 * to_yyyy_xy[k] * tbe_0 * tke_0;

            to_y_y_yyy_y[k] = 3.0 * to_yy_0[k] - 6.0 * to_yy_yy[k] * tke_0 - 2.0 * to_yyyy_0[k] * tbe_0 + 4.0 * to_yyyy_yy[k] * tbe_0 * tke_0;

            to_y_y_yyy_z[k] = -6.0 * to_yy_yz[k] * tke_0 + 4.0 * to_yyyy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 141-144 components of targeted buffer : FP

        auto to_y_y_yyz_x = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 21);

        auto to_y_y_yyz_y = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 22);

        auto to_y_y_yyz_z = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_y_y_yyz_x, to_y_y_yyz_y, to_y_y_yyz_z, to_yyyz_0, to_yyyz_xy, to_yyyz_yy, to_yyyz_yz, to_yz_0, to_yz_xy, to_yz_yy, to_yz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyz_x[k] = -4.0 * to_yz_xy[k] * tke_0 + 4.0 * to_yyyz_xy[k] * tbe_0 * tke_0;

            to_y_y_yyz_y[k] = 2.0 * to_yz_0[k] - 4.0 * to_yz_yy[k] * tke_0 - 2.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_yyyz_yy[k] * tbe_0 * tke_0;

            to_y_y_yyz_z[k] = -4.0 * to_yz_yz[k] * tke_0 + 4.0 * to_yyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 144-147 components of targeted buffer : FP

        auto to_y_y_yzz_x = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 24);

        auto to_y_y_yzz_y = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 25);

        auto to_y_y_yzz_z = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_y_y_yzz_x, to_y_y_yzz_y, to_y_y_yzz_z, to_yyzz_0, to_yyzz_xy, to_yyzz_yy, to_yyzz_yz, to_zz_0, to_zz_xy, to_zz_yy, to_zz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yzz_x[k] = -2.0 * to_zz_xy[k] * tke_0 + 4.0 * to_yyzz_xy[k] * tbe_0 * tke_0;

            to_y_y_yzz_y[k] = to_zz_0[k] - 2.0 * to_zz_yy[k] * tke_0 - 2.0 * to_yyzz_0[k] * tbe_0 + 4.0 * to_yyzz_yy[k] * tbe_0 * tke_0;

            to_y_y_yzz_z[k] = -2.0 * to_zz_yz[k] * tke_0 + 4.0 * to_yyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 147-150 components of targeted buffer : FP

        auto to_y_y_zzz_x = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 27);

        auto to_y_y_zzz_y = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 28);

        auto to_y_y_zzz_z = pbuffer.data(idx_op_geom_101_fp + 4 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_y_y_zzz_x, to_y_y_zzz_y, to_y_y_zzz_z, to_yzzz_0, to_yzzz_xy, to_yzzz_yy, to_yzzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zzz_x[k] = 4.0 * to_yzzz_xy[k] * tbe_0 * tke_0;

            to_y_y_zzz_y[k] = -2.0 * to_yzzz_0[k] * tbe_0 + 4.0 * to_yzzz_yy[k] * tbe_0 * tke_0;

            to_y_y_zzz_z[k] = 4.0 * to_yzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 150-153 components of targeted buffer : FP

        auto to_y_z_xxx_x = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 0);

        auto to_y_z_xxx_y = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 1);

        auto to_y_z_xxx_z = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_xxxy_0, to_xxxy_xz, to_xxxy_yz, to_xxxy_zz, to_y_z_xxx_x, to_y_z_xxx_y, to_y_z_xxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxx_x[k] = 4.0 * to_xxxy_xz[k] * tbe_0 * tke_0;

            to_y_z_xxx_y[k] = 4.0 * to_xxxy_yz[k] * tbe_0 * tke_0;

            to_y_z_xxx_z[k] = -2.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 153-156 components of targeted buffer : FP

        auto to_y_z_xxy_x = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 3);

        auto to_y_z_xxy_y = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 4);

        auto to_y_z_xxy_z = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_xx_0, to_xx_xz, to_xx_yz, to_xx_zz, to_xxyy_0, to_xxyy_xz, to_xxyy_yz, to_xxyy_zz, to_y_z_xxy_x, to_y_z_xxy_y, to_y_z_xxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxy_x[k] = -2.0 * to_xx_xz[k] * tke_0 + 4.0 * to_xxyy_xz[k] * tbe_0 * tke_0;

            to_y_z_xxy_y[k] = -2.0 * to_xx_yz[k] * tke_0 + 4.0 * to_xxyy_yz[k] * tbe_0 * tke_0;

            to_y_z_xxy_z[k] = to_xx_0[k] - 2.0 * to_xx_zz[k] * tke_0 - 2.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 156-159 components of targeted buffer : FP

        auto to_y_z_xxz_x = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 6);

        auto to_y_z_xxz_y = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 7);

        auto to_y_z_xxz_z = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_xxyz_0, to_xxyz_xz, to_xxyz_yz, to_xxyz_zz, to_y_z_xxz_x, to_y_z_xxz_y, to_y_z_xxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxz_x[k] = 4.0 * to_xxyz_xz[k] * tbe_0 * tke_0;

            to_y_z_xxz_y[k] = 4.0 * to_xxyz_yz[k] * tbe_0 * tke_0;

            to_y_z_xxz_z[k] = -2.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 159-162 components of targeted buffer : FP

        auto to_y_z_xyy_x = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 9);

        auto to_y_z_xyy_y = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 10);

        auto to_y_z_xyy_z = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_xy_0, to_xy_xz, to_xy_yz, to_xy_zz, to_xyyy_0, to_xyyy_xz, to_xyyy_yz, to_xyyy_zz, to_y_z_xyy_x, to_y_z_xyy_y, to_y_z_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyy_x[k] = -4.0 * to_xy_xz[k] * tke_0 + 4.0 * to_xyyy_xz[k] * tbe_0 * tke_0;

            to_y_z_xyy_y[k] = -4.0 * to_xy_yz[k] * tke_0 + 4.0 * to_xyyy_yz[k] * tbe_0 * tke_0;

            to_y_z_xyy_z[k] = 2.0 * to_xy_0[k] - 4.0 * to_xy_zz[k] * tke_0 - 2.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xyyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 162-165 components of targeted buffer : FP

        auto to_y_z_xyz_x = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 12);

        auto to_y_z_xyz_y = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 13);

        auto to_y_z_xyz_z = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_xyyz_0, to_xyyz_xz, to_xyyz_yz, to_xyyz_zz, to_xz_0, to_xz_xz, to_xz_yz, to_xz_zz, to_y_z_xyz_x, to_y_z_xyz_y, to_y_z_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyz_x[k] = -2.0 * to_xz_xz[k] * tke_0 + 4.0 * to_xyyz_xz[k] * tbe_0 * tke_0;

            to_y_z_xyz_y[k] = -2.0 * to_xz_yz[k] * tke_0 + 4.0 * to_xyyz_yz[k] * tbe_0 * tke_0;

            to_y_z_xyz_z[k] = to_xz_0[k] - 2.0 * to_xz_zz[k] * tke_0 - 2.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 165-168 components of targeted buffer : FP

        auto to_y_z_xzz_x = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 15);

        auto to_y_z_xzz_y = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 16);

        auto to_y_z_xzz_z = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_xyzz_0, to_xyzz_xz, to_xyzz_yz, to_xyzz_zz, to_y_z_xzz_x, to_y_z_xzz_y, to_y_z_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xzz_x[k] = 4.0 * to_xyzz_xz[k] * tbe_0 * tke_0;

            to_y_z_xzz_y[k] = 4.0 * to_xyzz_yz[k] * tbe_0 * tke_0;

            to_y_z_xzz_z[k] = -2.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 168-171 components of targeted buffer : FP

        auto to_y_z_yyy_x = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 18);

        auto to_y_z_yyy_y = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 19);

        auto to_y_z_yyy_z = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_y_z_yyy_x, to_y_z_yyy_y, to_y_z_yyy_z, to_yy_0, to_yy_xz, to_yy_yz, to_yy_zz, to_yyyy_0, to_yyyy_xz, to_yyyy_yz, to_yyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyy_x[k] = -6.0 * to_yy_xz[k] * tke_0 + 4.0 * to_yyyy_xz[k] * tbe_0 * tke_0;

            to_y_z_yyy_y[k] = -6.0 * to_yy_yz[k] * tke_0 + 4.0 * to_yyyy_yz[k] * tbe_0 * tke_0;

            to_y_z_yyy_z[k] = 3.0 * to_yy_0[k] - 6.0 * to_yy_zz[k] * tke_0 - 2.0 * to_yyyy_0[k] * tbe_0 + 4.0 * to_yyyy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 171-174 components of targeted buffer : FP

        auto to_y_z_yyz_x = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 21);

        auto to_y_z_yyz_y = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 22);

        auto to_y_z_yyz_z = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_y_z_yyz_x, to_y_z_yyz_y, to_y_z_yyz_z, to_yyyz_0, to_yyyz_xz, to_yyyz_yz, to_yyyz_zz, to_yz_0, to_yz_xz, to_yz_yz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyz_x[k] = -4.0 * to_yz_xz[k] * tke_0 + 4.0 * to_yyyz_xz[k] * tbe_0 * tke_0;

            to_y_z_yyz_y[k] = -4.0 * to_yz_yz[k] * tke_0 + 4.0 * to_yyyz_yz[k] * tbe_0 * tke_0;

            to_y_z_yyz_z[k] = 2.0 * to_yz_0[k] - 4.0 * to_yz_zz[k] * tke_0 - 2.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_yyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 174-177 components of targeted buffer : FP

        auto to_y_z_yzz_x = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 24);

        auto to_y_z_yzz_y = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 25);

        auto to_y_z_yzz_z = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_y_z_yzz_x, to_y_z_yzz_y, to_y_z_yzz_z, to_yyzz_0, to_yyzz_xz, to_yyzz_yz, to_yyzz_zz, to_zz_0, to_zz_xz, to_zz_yz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yzz_x[k] = -2.0 * to_zz_xz[k] * tke_0 + 4.0 * to_yyzz_xz[k] * tbe_0 * tke_0;

            to_y_z_yzz_y[k] = -2.0 * to_zz_yz[k] * tke_0 + 4.0 * to_yyzz_yz[k] * tbe_0 * tke_0;

            to_y_z_yzz_z[k] = to_zz_0[k] - 2.0 * to_zz_zz[k] * tke_0 - 2.0 * to_yyzz_0[k] * tbe_0 + 4.0 * to_yyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 177-180 components of targeted buffer : FP

        auto to_y_z_zzz_x = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 27);

        auto to_y_z_zzz_y = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 28);

        auto to_y_z_zzz_z = pbuffer.data(idx_op_geom_101_fp + 5 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_y_z_zzz_x, to_y_z_zzz_y, to_y_z_zzz_z, to_yzzz_0, to_yzzz_xz, to_yzzz_yz, to_yzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zzz_x[k] = 4.0 * to_yzzz_xz[k] * tbe_0 * tke_0;

            to_y_z_zzz_y[k] = 4.0 * to_yzzz_yz[k] * tbe_0 * tke_0;

            to_y_z_zzz_z[k] = -2.0 * to_yzzz_0[k] * tbe_0 + 4.0 * to_yzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 180-183 components of targeted buffer : FP

        auto to_z_x_xxx_x = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 0);

        auto to_z_x_xxx_y = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 1);

        auto to_z_x_xxx_z = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_xxxz_0, to_xxxz_xx, to_xxxz_xy, to_xxxz_xz, to_z_x_xxx_x, to_z_x_xxx_y, to_z_x_xxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxx_x[k] = -2.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxz_xx[k] * tbe_0 * tke_0;

            to_z_x_xxx_y[k] = 4.0 * to_xxxz_xy[k] * tbe_0 * tke_0;

            to_z_x_xxx_z[k] = 4.0 * to_xxxz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 183-186 components of targeted buffer : FP

        auto to_z_x_xxy_x = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 3);

        auto to_z_x_xxy_y = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 4);

        auto to_z_x_xxy_z = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_xxyz_0, to_xxyz_xx, to_xxyz_xy, to_xxyz_xz, to_z_x_xxy_x, to_z_x_xxy_y, to_z_x_xxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxy_x[k] = -2.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyz_xx[k] * tbe_0 * tke_0;

            to_z_x_xxy_y[k] = 4.0 * to_xxyz_xy[k] * tbe_0 * tke_0;

            to_z_x_xxy_z[k] = 4.0 * to_xxyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 186-189 components of targeted buffer : FP

        auto to_z_x_xxz_x = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 6);

        auto to_z_x_xxz_y = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 7);

        auto to_z_x_xxz_z = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_xx_0, to_xx_xx, to_xx_xy, to_xx_xz, to_xxzz_0, to_xxzz_xx, to_xxzz_xy, to_xxzz_xz, to_z_x_xxz_x, to_z_x_xxz_y, to_z_x_xxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxz_x[k] = to_xx_0[k] - 2.0 * to_xx_xx[k] * tke_0 - 2.0 * to_xxzz_0[k] * tbe_0 + 4.0 * to_xxzz_xx[k] * tbe_0 * tke_0;

            to_z_x_xxz_y[k] = -2.0 * to_xx_xy[k] * tke_0 + 4.0 * to_xxzz_xy[k] * tbe_0 * tke_0;

            to_z_x_xxz_z[k] = -2.0 * to_xx_xz[k] * tke_0 + 4.0 * to_xxzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 189-192 components of targeted buffer : FP

        auto to_z_x_xyy_x = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 9);

        auto to_z_x_xyy_y = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 10);

        auto to_z_x_xyy_z = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_xyyz_0, to_xyyz_xx, to_xyyz_xy, to_xyyz_xz, to_z_x_xyy_x, to_z_x_xyy_y, to_z_x_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyy_x[k] = -2.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyz_xx[k] * tbe_0 * tke_0;

            to_z_x_xyy_y[k] = 4.0 * to_xyyz_xy[k] * tbe_0 * tke_0;

            to_z_x_xyy_z[k] = 4.0 * to_xyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 192-195 components of targeted buffer : FP

        auto to_z_x_xyz_x = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 12);

        auto to_z_x_xyz_y = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 13);

        auto to_z_x_xyz_z = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_xy_0, to_xy_xx, to_xy_xy, to_xy_xz, to_xyzz_0, to_xyzz_xx, to_xyzz_xy, to_xyzz_xz, to_z_x_xyz_x, to_z_x_xyz_y, to_z_x_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyz_x[k] = to_xy_0[k] - 2.0 * to_xy_xx[k] * tke_0 - 2.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzz_xx[k] * tbe_0 * tke_0;

            to_z_x_xyz_y[k] = -2.0 * to_xy_xy[k] * tke_0 + 4.0 * to_xyzz_xy[k] * tbe_0 * tke_0;

            to_z_x_xyz_z[k] = -2.0 * to_xy_xz[k] * tke_0 + 4.0 * to_xyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 195-198 components of targeted buffer : FP

        auto to_z_x_xzz_x = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 15);

        auto to_z_x_xzz_y = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 16);

        auto to_z_x_xzz_z = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_xz_0, to_xz_xx, to_xz_xy, to_xz_xz, to_xzzz_0, to_xzzz_xx, to_xzzz_xy, to_xzzz_xz, to_z_x_xzz_x, to_z_x_xzz_y, to_z_x_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xzz_x[k] = 2.0 * to_xz_0[k] - 4.0 * to_xz_xx[k] * tke_0 - 2.0 * to_xzzz_0[k] * tbe_0 + 4.0 * to_xzzz_xx[k] * tbe_0 * tke_0;

            to_z_x_xzz_y[k] = -4.0 * to_xz_xy[k] * tke_0 + 4.0 * to_xzzz_xy[k] * tbe_0 * tke_0;

            to_z_x_xzz_z[k] = -4.0 * to_xz_xz[k] * tke_0 + 4.0 * to_xzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 198-201 components of targeted buffer : FP

        auto to_z_x_yyy_x = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 18);

        auto to_z_x_yyy_y = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 19);

        auto to_z_x_yyy_z = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_yyyz_0, to_yyyz_xx, to_yyyz_xy, to_yyyz_xz, to_z_x_yyy_x, to_z_x_yyy_y, to_z_x_yyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyy_x[k] = -2.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_yyyz_xx[k] * tbe_0 * tke_0;

            to_z_x_yyy_y[k] = 4.0 * to_yyyz_xy[k] * tbe_0 * tke_0;

            to_z_x_yyy_z[k] = 4.0 * to_yyyz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 201-204 components of targeted buffer : FP

        auto to_z_x_yyz_x = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 21);

        auto to_z_x_yyz_y = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 22);

        auto to_z_x_yyz_z = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_yy_0, to_yy_xx, to_yy_xy, to_yy_xz, to_yyzz_0, to_yyzz_xx, to_yyzz_xy, to_yyzz_xz, to_z_x_yyz_x, to_z_x_yyz_y, to_z_x_yyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyz_x[k] = to_yy_0[k] - 2.0 * to_yy_xx[k] * tke_0 - 2.0 * to_yyzz_0[k] * tbe_0 + 4.0 * to_yyzz_xx[k] * tbe_0 * tke_0;

            to_z_x_yyz_y[k] = -2.0 * to_yy_xy[k] * tke_0 + 4.0 * to_yyzz_xy[k] * tbe_0 * tke_0;

            to_z_x_yyz_z[k] = -2.0 * to_yy_xz[k] * tke_0 + 4.0 * to_yyzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 204-207 components of targeted buffer : FP

        auto to_z_x_yzz_x = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 24);

        auto to_z_x_yzz_y = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 25);

        auto to_z_x_yzz_z = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_yz_0, to_yz_xx, to_yz_xy, to_yz_xz, to_yzzz_0, to_yzzz_xx, to_yzzz_xy, to_yzzz_xz, to_z_x_yzz_x, to_z_x_yzz_y, to_z_x_yzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yzz_x[k] = 2.0 * to_yz_0[k] - 4.0 * to_yz_xx[k] * tke_0 - 2.0 * to_yzzz_0[k] * tbe_0 + 4.0 * to_yzzz_xx[k] * tbe_0 * tke_0;

            to_z_x_yzz_y[k] = -4.0 * to_yz_xy[k] * tke_0 + 4.0 * to_yzzz_xy[k] * tbe_0 * tke_0;

            to_z_x_yzz_z[k] = -4.0 * to_yz_xz[k] * tke_0 + 4.0 * to_yzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 207-210 components of targeted buffer : FP

        auto to_z_x_zzz_x = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 27);

        auto to_z_x_zzz_y = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 28);

        auto to_z_x_zzz_z = pbuffer.data(idx_op_geom_101_fp + 6 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_z_x_zzz_x, to_z_x_zzz_y, to_z_x_zzz_z, to_zz_0, to_zz_xx, to_zz_xy, to_zz_xz, to_zzzz_0, to_zzzz_xx, to_zzzz_xy, to_zzzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zzz_x[k] = 3.0 * to_zz_0[k] - 6.0 * to_zz_xx[k] * tke_0 - 2.0 * to_zzzz_0[k] * tbe_0 + 4.0 * to_zzzz_xx[k] * tbe_0 * tke_0;

            to_z_x_zzz_y[k] = -6.0 * to_zz_xy[k] * tke_0 + 4.0 * to_zzzz_xy[k] * tbe_0 * tke_0;

            to_z_x_zzz_z[k] = -6.0 * to_zz_xz[k] * tke_0 + 4.0 * to_zzzz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 210-213 components of targeted buffer : FP

        auto to_z_y_xxx_x = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 0);

        auto to_z_y_xxx_y = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 1);

        auto to_z_y_xxx_z = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_xxxz_0, to_xxxz_xy, to_xxxz_yy, to_xxxz_yz, to_z_y_xxx_x, to_z_y_xxx_y, to_z_y_xxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxx_x[k] = 4.0 * to_xxxz_xy[k] * tbe_0 * tke_0;

            to_z_y_xxx_y[k] = -2.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxz_yy[k] * tbe_0 * tke_0;

            to_z_y_xxx_z[k] = 4.0 * to_xxxz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 213-216 components of targeted buffer : FP

        auto to_z_y_xxy_x = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 3);

        auto to_z_y_xxy_y = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 4);

        auto to_z_y_xxy_z = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_xxyz_0, to_xxyz_xy, to_xxyz_yy, to_xxyz_yz, to_z_y_xxy_x, to_z_y_xxy_y, to_z_y_xxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxy_x[k] = 4.0 * to_xxyz_xy[k] * tbe_0 * tke_0;

            to_z_y_xxy_y[k] = -2.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyz_yy[k] * tbe_0 * tke_0;

            to_z_y_xxy_z[k] = 4.0 * to_xxyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 216-219 components of targeted buffer : FP

        auto to_z_y_xxz_x = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 6);

        auto to_z_y_xxz_y = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 7);

        auto to_z_y_xxz_z = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_xx_0, to_xx_xy, to_xx_yy, to_xx_yz, to_xxzz_0, to_xxzz_xy, to_xxzz_yy, to_xxzz_yz, to_z_y_xxz_x, to_z_y_xxz_y, to_z_y_xxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxz_x[k] = -2.0 * to_xx_xy[k] * tke_0 + 4.0 * to_xxzz_xy[k] * tbe_0 * tke_0;

            to_z_y_xxz_y[k] = to_xx_0[k] - 2.0 * to_xx_yy[k] * tke_0 - 2.0 * to_xxzz_0[k] * tbe_0 + 4.0 * to_xxzz_yy[k] * tbe_0 * tke_0;

            to_z_y_xxz_z[k] = -2.0 * to_xx_yz[k] * tke_0 + 4.0 * to_xxzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 219-222 components of targeted buffer : FP

        auto to_z_y_xyy_x = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 9);

        auto to_z_y_xyy_y = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 10);

        auto to_z_y_xyy_z = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_xyyz_0, to_xyyz_xy, to_xyyz_yy, to_xyyz_yz, to_z_y_xyy_x, to_z_y_xyy_y, to_z_y_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyy_x[k] = 4.0 * to_xyyz_xy[k] * tbe_0 * tke_0;

            to_z_y_xyy_y[k] = -2.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyz_yy[k] * tbe_0 * tke_0;

            to_z_y_xyy_z[k] = 4.0 * to_xyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 222-225 components of targeted buffer : FP

        auto to_z_y_xyz_x = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 12);

        auto to_z_y_xyz_y = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 13);

        auto to_z_y_xyz_z = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_xy_0, to_xy_xy, to_xy_yy, to_xy_yz, to_xyzz_0, to_xyzz_xy, to_xyzz_yy, to_xyzz_yz, to_z_y_xyz_x, to_z_y_xyz_y, to_z_y_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyz_x[k] = -2.0 * to_xy_xy[k] * tke_0 + 4.0 * to_xyzz_xy[k] * tbe_0 * tke_0;

            to_z_y_xyz_y[k] = to_xy_0[k] - 2.0 * to_xy_yy[k] * tke_0 - 2.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzz_yy[k] * tbe_0 * tke_0;

            to_z_y_xyz_z[k] = -2.0 * to_xy_yz[k] * tke_0 + 4.0 * to_xyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 225-228 components of targeted buffer : FP

        auto to_z_y_xzz_x = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 15);

        auto to_z_y_xzz_y = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 16);

        auto to_z_y_xzz_z = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_xz_0, to_xz_xy, to_xz_yy, to_xz_yz, to_xzzz_0, to_xzzz_xy, to_xzzz_yy, to_xzzz_yz, to_z_y_xzz_x, to_z_y_xzz_y, to_z_y_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xzz_x[k] = -4.0 * to_xz_xy[k] * tke_0 + 4.0 * to_xzzz_xy[k] * tbe_0 * tke_0;

            to_z_y_xzz_y[k] = 2.0 * to_xz_0[k] - 4.0 * to_xz_yy[k] * tke_0 - 2.0 * to_xzzz_0[k] * tbe_0 + 4.0 * to_xzzz_yy[k] * tbe_0 * tke_0;

            to_z_y_xzz_z[k] = -4.0 * to_xz_yz[k] * tke_0 + 4.0 * to_xzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 228-231 components of targeted buffer : FP

        auto to_z_y_yyy_x = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 18);

        auto to_z_y_yyy_y = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 19);

        auto to_z_y_yyy_z = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_yyyz_0, to_yyyz_xy, to_yyyz_yy, to_yyyz_yz, to_z_y_yyy_x, to_z_y_yyy_y, to_z_y_yyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyy_x[k] = 4.0 * to_yyyz_xy[k] * tbe_0 * tke_0;

            to_z_y_yyy_y[k] = -2.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_yyyz_yy[k] * tbe_0 * tke_0;

            to_z_y_yyy_z[k] = 4.0 * to_yyyz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 231-234 components of targeted buffer : FP

        auto to_z_y_yyz_x = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 21);

        auto to_z_y_yyz_y = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 22);

        auto to_z_y_yyz_z = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_yy_0, to_yy_xy, to_yy_yy, to_yy_yz, to_yyzz_0, to_yyzz_xy, to_yyzz_yy, to_yyzz_yz, to_z_y_yyz_x, to_z_y_yyz_y, to_z_y_yyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyz_x[k] = -2.0 * to_yy_xy[k] * tke_0 + 4.0 * to_yyzz_xy[k] * tbe_0 * tke_0;

            to_z_y_yyz_y[k] = to_yy_0[k] - 2.0 * to_yy_yy[k] * tke_0 - 2.0 * to_yyzz_0[k] * tbe_0 + 4.0 * to_yyzz_yy[k] * tbe_0 * tke_0;

            to_z_y_yyz_z[k] = -2.0 * to_yy_yz[k] * tke_0 + 4.0 * to_yyzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 234-237 components of targeted buffer : FP

        auto to_z_y_yzz_x = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 24);

        auto to_z_y_yzz_y = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 25);

        auto to_z_y_yzz_z = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_yz_0, to_yz_xy, to_yz_yy, to_yz_yz, to_yzzz_0, to_yzzz_xy, to_yzzz_yy, to_yzzz_yz, to_z_y_yzz_x, to_z_y_yzz_y, to_z_y_yzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yzz_x[k] = -4.0 * to_yz_xy[k] * tke_0 + 4.0 * to_yzzz_xy[k] * tbe_0 * tke_0;

            to_z_y_yzz_y[k] = 2.0 * to_yz_0[k] - 4.0 * to_yz_yy[k] * tke_0 - 2.0 * to_yzzz_0[k] * tbe_0 + 4.0 * to_yzzz_yy[k] * tbe_0 * tke_0;

            to_z_y_yzz_z[k] = -4.0 * to_yz_yz[k] * tke_0 + 4.0 * to_yzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 237-240 components of targeted buffer : FP

        auto to_z_y_zzz_x = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 27);

        auto to_z_y_zzz_y = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 28);

        auto to_z_y_zzz_z = pbuffer.data(idx_op_geom_101_fp + 7 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_z_y_zzz_x, to_z_y_zzz_y, to_z_y_zzz_z, to_zz_0, to_zz_xy, to_zz_yy, to_zz_yz, to_zzzz_0, to_zzzz_xy, to_zzzz_yy, to_zzzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zzz_x[k] = -6.0 * to_zz_xy[k] * tke_0 + 4.0 * to_zzzz_xy[k] * tbe_0 * tke_0;

            to_z_y_zzz_y[k] = 3.0 * to_zz_0[k] - 6.0 * to_zz_yy[k] * tke_0 - 2.0 * to_zzzz_0[k] * tbe_0 + 4.0 * to_zzzz_yy[k] * tbe_0 * tke_0;

            to_z_y_zzz_z[k] = -6.0 * to_zz_yz[k] * tke_0 + 4.0 * to_zzzz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 240-243 components of targeted buffer : FP

        auto to_z_z_xxx_x = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 0);

        auto to_z_z_xxx_y = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 1);

        auto to_z_z_xxx_z = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 2);

        #pragma omp simd aligned(to_xxxz_0, to_xxxz_xz, to_xxxz_yz, to_xxxz_zz, to_z_z_xxx_x, to_z_z_xxx_y, to_z_z_xxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxx_x[k] = 4.0 * to_xxxz_xz[k] * tbe_0 * tke_0;

            to_z_z_xxx_y[k] = 4.0 * to_xxxz_yz[k] * tbe_0 * tke_0;

            to_z_z_xxx_z[k] = -2.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 243-246 components of targeted buffer : FP

        auto to_z_z_xxy_x = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 3);

        auto to_z_z_xxy_y = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 4);

        auto to_z_z_xxy_z = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 5);

        #pragma omp simd aligned(to_xxyz_0, to_xxyz_xz, to_xxyz_yz, to_xxyz_zz, to_z_z_xxy_x, to_z_z_xxy_y, to_z_z_xxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxy_x[k] = 4.0 * to_xxyz_xz[k] * tbe_0 * tke_0;

            to_z_z_xxy_y[k] = 4.0 * to_xxyz_yz[k] * tbe_0 * tke_0;

            to_z_z_xxy_z[k] = -2.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 246-249 components of targeted buffer : FP

        auto to_z_z_xxz_x = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 6);

        auto to_z_z_xxz_y = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 7);

        auto to_z_z_xxz_z = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 8);

        #pragma omp simd aligned(to_xx_0, to_xx_xz, to_xx_yz, to_xx_zz, to_xxzz_0, to_xxzz_xz, to_xxzz_yz, to_xxzz_zz, to_z_z_xxz_x, to_z_z_xxz_y, to_z_z_xxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxz_x[k] = -2.0 * to_xx_xz[k] * tke_0 + 4.0 * to_xxzz_xz[k] * tbe_0 * tke_0;

            to_z_z_xxz_y[k] = -2.0 * to_xx_yz[k] * tke_0 + 4.0 * to_xxzz_yz[k] * tbe_0 * tke_0;

            to_z_z_xxz_z[k] = to_xx_0[k] - 2.0 * to_xx_zz[k] * tke_0 - 2.0 * to_xxzz_0[k] * tbe_0 + 4.0 * to_xxzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 249-252 components of targeted buffer : FP

        auto to_z_z_xyy_x = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 9);

        auto to_z_z_xyy_y = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 10);

        auto to_z_z_xyy_z = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 11);

        #pragma omp simd aligned(to_xyyz_0, to_xyyz_xz, to_xyyz_yz, to_xyyz_zz, to_z_z_xyy_x, to_z_z_xyy_y, to_z_z_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyy_x[k] = 4.0 * to_xyyz_xz[k] * tbe_0 * tke_0;

            to_z_z_xyy_y[k] = 4.0 * to_xyyz_yz[k] * tbe_0 * tke_0;

            to_z_z_xyy_z[k] = -2.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 252-255 components of targeted buffer : FP

        auto to_z_z_xyz_x = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 12);

        auto to_z_z_xyz_y = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 13);

        auto to_z_z_xyz_z = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 14);

        #pragma omp simd aligned(to_xy_0, to_xy_xz, to_xy_yz, to_xy_zz, to_xyzz_0, to_xyzz_xz, to_xyzz_yz, to_xyzz_zz, to_z_z_xyz_x, to_z_z_xyz_y, to_z_z_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyz_x[k] = -2.0 * to_xy_xz[k] * tke_0 + 4.0 * to_xyzz_xz[k] * tbe_0 * tke_0;

            to_z_z_xyz_y[k] = -2.0 * to_xy_yz[k] * tke_0 + 4.0 * to_xyzz_yz[k] * tbe_0 * tke_0;

            to_z_z_xyz_z[k] = to_xy_0[k] - 2.0 * to_xy_zz[k] * tke_0 - 2.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 255-258 components of targeted buffer : FP

        auto to_z_z_xzz_x = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 15);

        auto to_z_z_xzz_y = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 16);

        auto to_z_z_xzz_z = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 17);

        #pragma omp simd aligned(to_xz_0, to_xz_xz, to_xz_yz, to_xz_zz, to_xzzz_0, to_xzzz_xz, to_xzzz_yz, to_xzzz_zz, to_z_z_xzz_x, to_z_z_xzz_y, to_z_z_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xzz_x[k] = -4.0 * to_xz_xz[k] * tke_0 + 4.0 * to_xzzz_xz[k] * tbe_0 * tke_0;

            to_z_z_xzz_y[k] = -4.0 * to_xz_yz[k] * tke_0 + 4.0 * to_xzzz_yz[k] * tbe_0 * tke_0;

            to_z_z_xzz_z[k] = 2.0 * to_xz_0[k] - 4.0 * to_xz_zz[k] * tke_0 - 2.0 * to_xzzz_0[k] * tbe_0 + 4.0 * to_xzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 258-261 components of targeted buffer : FP

        auto to_z_z_yyy_x = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 18);

        auto to_z_z_yyy_y = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 19);

        auto to_z_z_yyy_z = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 20);

        #pragma omp simd aligned(to_yyyz_0, to_yyyz_xz, to_yyyz_yz, to_yyyz_zz, to_z_z_yyy_x, to_z_z_yyy_y, to_z_z_yyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyy_x[k] = 4.0 * to_yyyz_xz[k] * tbe_0 * tke_0;

            to_z_z_yyy_y[k] = 4.0 * to_yyyz_yz[k] * tbe_0 * tke_0;

            to_z_z_yyy_z[k] = -2.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_yyyz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 261-264 components of targeted buffer : FP

        auto to_z_z_yyz_x = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 21);

        auto to_z_z_yyz_y = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 22);

        auto to_z_z_yyz_z = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 23);

        #pragma omp simd aligned(to_yy_0, to_yy_xz, to_yy_yz, to_yy_zz, to_yyzz_0, to_yyzz_xz, to_yyzz_yz, to_yyzz_zz, to_z_z_yyz_x, to_z_z_yyz_y, to_z_z_yyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyz_x[k] = -2.0 * to_yy_xz[k] * tke_0 + 4.0 * to_yyzz_xz[k] * tbe_0 * tke_0;

            to_z_z_yyz_y[k] = -2.0 * to_yy_yz[k] * tke_0 + 4.0 * to_yyzz_yz[k] * tbe_0 * tke_0;

            to_z_z_yyz_z[k] = to_yy_0[k] - 2.0 * to_yy_zz[k] * tke_0 - 2.0 * to_yyzz_0[k] * tbe_0 + 4.0 * to_yyzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 264-267 components of targeted buffer : FP

        auto to_z_z_yzz_x = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 24);

        auto to_z_z_yzz_y = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 25);

        auto to_z_z_yzz_z = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 26);

        #pragma omp simd aligned(to_yz_0, to_yz_xz, to_yz_yz, to_yz_zz, to_yzzz_0, to_yzzz_xz, to_yzzz_yz, to_yzzz_zz, to_z_z_yzz_x, to_z_z_yzz_y, to_z_z_yzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yzz_x[k] = -4.0 * to_yz_xz[k] * tke_0 + 4.0 * to_yzzz_xz[k] * tbe_0 * tke_0;

            to_z_z_yzz_y[k] = -4.0 * to_yz_yz[k] * tke_0 + 4.0 * to_yzzz_yz[k] * tbe_0 * tke_0;

            to_z_z_yzz_z[k] = 2.0 * to_yz_0[k] - 4.0 * to_yz_zz[k] * tke_0 - 2.0 * to_yzzz_0[k] * tbe_0 + 4.0 * to_yzzz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 267-270 components of targeted buffer : FP

        auto to_z_z_zzz_x = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 27);

        auto to_z_z_zzz_y = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 28);

        auto to_z_z_zzz_z = pbuffer.data(idx_op_geom_101_fp + 8 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_z_z_zzz_x, to_z_z_zzz_y, to_z_z_zzz_z, to_zz_0, to_zz_xz, to_zz_yz, to_zz_zz, to_zzzz_0, to_zzzz_xz, to_zzzz_yz, to_zzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zzz_x[k] = -6.0 * to_zz_xz[k] * tke_0 + 4.0 * to_zzzz_xz[k] * tbe_0 * tke_0;

            to_z_z_zzz_y[k] = -6.0 * to_zz_yz[k] * tke_0 + 4.0 * to_zzzz_yz[k] * tbe_0 * tke_0;

            to_z_z_zzz_z[k] = 3.0 * to_zz_0[k] - 6.0 * to_zz_zz[k] * tke_0 - 2.0 * to_zzzz_0[k] * tbe_0 + 4.0 * to_zzzz_zz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

